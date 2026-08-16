"""
phagecore.lysis
==============
Detection harness for lysis-cassette and host-recognition proteins.

Everything here is profile-driven: keyword sets, length priors, markers and the
curation registry all come from the HostProfile. The engine logic (collect ->
classify -> rank -> select -> optional tBLASTn rescue) is host-agnostic.

Phase 2 module registry
-----------------------
Lysis biology differs by host. The current modules ("endolysin", "holin") run
for every Gram-positive host. Gram-negative hosts additionally need "spanin"
(outer-membrane disruption); Mycobacterium needs "lysinB" (mycolic-acid
esterase). Those are registered as stubs below and switched on by
profile.active_lysis_modules — written once in the engine, enabled per host,
never by forking the codebase.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from Bio.SeqRecord import SeqRecord

from . import (CLASS_FREE, CLASS_VAPH, CLASS_NONLYTIC, CLASS_UNCERTAIN,
               CLASS_NONENDO, normalise_class)

log = logging.getLogger("phagecore.lysis")


# ---------------------------------------------------------------------------
# Qualifier helpers
# ---------------------------------------------------------------------------

def qualifier_text(feature) -> str:
    """Concatenate product/gene/note/function qualifiers, lowercased."""
    parts = []
    for key in ("product", "gene", "note", "function"):
        parts.extend(feature.qualifiers.get(key, []))
    return " ".join(parts).lower()


def identity_text(feature) -> str:
    """Product + gene only (the identity fields), lowercased.

    Used for the non-endolysin denylist so that a real endolysin whose *note*
    merely mentions a neighbouring gene (e.g. "upstream of dUTPase") is not
    wrongly excluded — only its own product/gene identity is tested.
    """
    parts = []
    for key in ("product", "gene"):
        parts.extend(feature.qualifiers.get(key, []))
    return " ".join(parts).lower()


def _kw_hit(text: str, keyword: str) -> bool:
    """
    Substring match, except bare alphanumeric tokens <=3 chars (e.g. 'rbp',
    'm15') are matched on word boundaries to avoid false positives inside longer
    words.
    """
    if len(keyword) <= 3 and keyword.isalnum():
        return re.search(rf"\b{re.escape(keyword)}\b", text) is not None
    return keyword in text


def any_kw(text: str, keywords) -> Optional[str]:
    for kw in keywords:
        if _kw_hit(text, kw):
            return kw
    return None


# ===========================================================================
# S2 surface: holin & RBP presence
# ===========================================================================

@dataclass
class HolinRBP:
    holin_present: bool
    holin_evidence: str
    rbp_present: bool
    rbp_evidence: str
    rbp_count: int          # number of distinct CDS matching an RBP keyword


def detect_holin_rbp(record: SeqRecord, profile) -> HolinRBP:
    """
    Detect holin (first hit) and RBP (first hit + total count).

    Improvement over the original: RBP multiplicity is reported (phages often
    carry several tail-fibre/RBP genes — directly relevant to host-range work),
    and bare-token keywords use word-boundary matching.
    """
    holin_ev = rbp_ev = None
    rbp_count = 0
    for f in record.features:
        if f.type != "CDS":
            continue
        text = qualifier_text(f)
        if holin_ev is None:
            kw = any_kw(text, profile.holin_keywords)
            if kw:
                prod = f.qualifiers.get("product", ["(no product)"])[0]
                holin_ev = (kw, prod)
        kw = any_kw(text, profile.rbp_keywords)
        if kw:
            rbp_count += 1
            if rbp_ev is None:
                prod = f.qualifiers.get("product", ["(no product)"])[0]
                rbp_ev = (kw, prod)

    def fmt(ev):
        return f"[{ev[0]}] {ev[1]}" if ev else "—"

    return HolinRBP(
        holin_present=holin_ev is not None, holin_evidence=fmt(holin_ev),
        rbp_present=rbp_ev is not None, rbp_evidence=fmt(rbp_ev),
        rbp_count=rbp_count,
    )


# ===========================================================================
# S4 surface: endolysin candidates
# ===========================================================================

@dataclass
class Candidate:
    protein_id: str
    accession: str
    organism: str
    product: str
    sequence: str
    inferred_domain: str = "unknown"     # host display label (from profile)
    domain_key: str = "unknown"          # host-neutral evidence key (engine)
    classification: str = CLASS_UNCERTAIN
    selected: bool = False
    tblastn_identity: Optional[float] = None
    tblastn_note: str = ""
    runtime_flag: str = ""
    evidence: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def fasta_header(self) -> str:
        org = self.organism.replace(" ", "_")
        prod = self.product.replace("|", "/")
        # accession first => survives as the InterProScan sequence ID for a
        # clean, reproducible join when reconciling InterPro results (--interpro).
        return (f"{self.accession}|{self.protein_id}|{org}|"
                f"status={self.classification}|{prod}")


# ---------------------------------------------------------------------------
# Neutral evidence vocabulary (v3.0)
# ---------------------------------------------------------------------------
# infer_evidence() emits ONE of these HOST-NEUTRAL keys. It performs domain
# DETECTION only; it never names a host-specific class ("LysK-type", "NADAR").
# Naming is a host decision: the profile maps each key -> display label
# (profile.domain_labels), and the profile declares which keys count as a
# free endolysin (profile.lytic_domain_terms). This breaks the v2.x coupling
# where a Staph-flavoured LABEL ("LysK-type") was load-bearing for the
# free-endolysin CLASSIFICATION — the bug that dropped 46 B. subtilis
# endolysins to 'uncertain' when a non-Staph profile (correctly) lacked "lysk".
EV_NONLYTIC = "nonlytic"
EV_TAIL_LYSOZYME = "tail_lysozyme"
EV_CHAP_AMIDASE = "chap_amidase"
EV_CHAP = "chap"
EV_AMIDASE = "amidase"
EV_PEPTIDASE_M15 = "peptidase_m15"
EV_GLUCOSAMINIDASE = "glucosaminidase"
EV_NLPC_P60 = "nlpc_p60"
EV_MURAMIDASE = "muramidase"
EV_GLYCOSIDASE = "glycosidase"
EV_ENDOLYSIN_GENERIC = "endolysin_generic"
EV_PEPTIDOGLYCAN_HYDROLASE = "peptidoglycan_hydrolase"
EV_PEPTIDASE_GENERIC = "peptidase_generic"
EV_HYDROLASE_GENERIC = "hydrolase_generic"
EV_UNKNOWN = "unknown"

EVIDENCE_KEYS = frozenset({
    EV_NONLYTIC, EV_TAIL_LYSOZYME, EV_CHAP_AMIDASE, EV_CHAP, EV_AMIDASE,
    EV_PEPTIDASE_M15, EV_GLUCOSAMINIDASE, EV_NLPC_P60, EV_MURAMIDASE,
    EV_GLYCOSIDASE, EV_ENDOLYSIN_GENERIC, EV_PEPTIDOGLYCAN_HYDROLASE,
    EV_PEPTIDASE_GENERIC, EV_HYDROLASE_GENERIC, EV_UNKNOWN,
})

# Engine fallback labels (host-neutral). A profile SHOULD override these via
# profile.domain_labels; the Staph profile does so to reproduce its exact v2.x
# display strings byte-for-byte. Detection vocabulary (chap/amidase/muramidase
# ...) is the pan-phage peptidoglycan-hydrolase core and stays in the engine;
# only the NAMING and the free/not-free decision are host data.
_DEFAULT_DOMAIN_LABELS = {
    EV_NONLYTIC: "Non-lytic (marker)",
    EV_TAIL_LYSOZYME: "Tail-associated lysozyme",
    EV_CHAP_AMIDASE: "CHAP + Amidase",
    EV_CHAP: "CHAP",
    EV_AMIDASE: "Amidase",
    EV_PEPTIDASE_M15: "Peptidase_M15",
    EV_GLUCOSAMINIDASE: "Glucosaminidase",
    EV_NLPC_P60: "NlpC/P60",
    EV_MURAMIDASE: "Muramidase",
    EV_GLYCOSIDASE: "Glycosidase",
    EV_ENDOLYSIN_GENERIC: "Endolysin (unspecified)",
    EV_PEPTIDOGLYCAN_HYDROLASE: "Peptidoglycan hydrolase",
    EV_PEPTIDASE_GENERIC: "Peptidase (unspecified — verify)",
    EV_HYDROLASE_GENERIC: "Hydrolase (unspecified — verify; may be non-PG)",
    EV_UNKNOWN: "unknown",
}


def infer_evidence(text: str, profile) -> str:
    """Host-neutral domain DETECTION. Returns an EV_* evidence key.

    Decision tree is identical to the v2.x infer_domain() so that, once the
    profile maps keys->labels, the display column and classification reproduce
    v2.x byte-for-byte on the Staph set. The ONLY change is that this returns a
    neutral key instead of a Staph-flavoured display string.
    """
    if any(m in text for m in profile.nonlytic_markers):
        return EV_NONLYTIC
    if "phage_lysozyme2" in text or ("tail" in text and "lysozyme" in text):
        return EV_TAIL_LYSOZYME
    has_chap = "chap" in text
    has_ami = "amidase" in text or "n-acetylmuramoyl" in text
    if has_chap and has_ami:
        return EV_CHAP_AMIDASE
    if has_chap:
        return EV_CHAP
    if has_ami:
        return EV_AMIDASE
    if _kw_hit(text, "m15"):
        return EV_PEPTIDASE_M15
    if "glucosaminidase" in text:
        return EV_GLUCOSAMINIDASE
    if "nlpc" in text or "p60" in text:
        return EV_NLPC_P60
    if "muramidase" in text or "lysozyme" in text:
        return EV_MURAMIDASE
    if "glycosidase" in text:
        return EV_GLYCOSIDASE
    if "lysk" in text or "endolysin" in text or "lysin" in text:
        return EV_ENDOLYSIN_GENERIC
    if "peptidoglycan" in text:
        return EV_PEPTIDOGLYCAN_HYDROLASE
    if "peptidase" in text:
        return EV_PEPTIDASE_GENERIC
    if "hydrolase" in text:
        # bare "hydrolase" is NOT necessarily a peptidoglycan hydrolase
        # (dUTPase = nucleotidohydrolase, RNase Z = metal-dependent hydrolase).
        return EV_HYDROLASE_GENERIC
    return EV_UNKNOWN


def domain_label(evidence_key: str, profile) -> str:
    """Interpret a neutral evidence key into a host display label.

    Host layer owns naming: profile.domain_labels overrides the engine default.
    Falls back to the neutral default for any key a profile does not map.
    """
    labels = getattr(profile, "domain_labels", None) or {}
    return labels.get(evidence_key,
                       _DEFAULT_DOMAIN_LABELS.get(evidence_key, evidence_key))


def classify_candidate(text: str, length: int, domain_key: str,
                       profile, ident_text: str = "") -> str:
    """Per-candidate classification (canonical lowercase vocabulary).

    v3.0: the free-endolysin decision is read from the PROFILE, not a hardcoded
    tuple. A candidate is a free endolysin iff its neutral evidence key (from
    infer_evidence) is listed in profile.lytic_domain_terms. This is the live
    knob the v2.x field only pretended to be: editing lytic_domain_terms now
    changes classification. Order (denylist -> non-lytic -> VAPH -> free) and
    the VAPH/denylist logic are unchanged, so Staph output is byte-identical
    once STAPH_AUREUS.lytic_domain_terms carries the former hardcoded set.

    `ident_text` is the product/gene identity; if it matches the profile's
    non-endolysin denylist the candidate is tagged CLASS_NONENDO and can never
    be selected as the endolysin representative (dUTPase, chaperonin, RNase,
    NTPase, tail spike — broad-keyword false positives confirmed by InterPro).
    """
    probe = ident_text or text
    if profile.non_endolysin_markers and \
            any(m in probe for m in profile.non_endolysin_markers):
        return CLASS_NONENDO
    if domain_key == EV_NONLYTIC:
        return CLASS_NONLYTIC
    if any(m in text for m in profile.vaph_markers) or \
            length > profile.vaph_length_threshold_aa:
        return CLASS_VAPH
    if domain_key in profile.lytic_domain_terms:
        return CLASS_FREE
    return CLASS_UNCERTAIN


def collect_candidates(record: SeqRecord, profile) -> list[Candidate]:
    """Collect & classify every lysis-keyword CDS in a genome."""
    from phagecore.genbank_io import resolve_organism
    organism = resolve_organism(record, getattr(record, "_source_path", None))
    accession = record.id
    acc_base = accession.split(".")[0]
    runtime = ("intron-split" if acc_base in profile.intron_split_accessions else
               "hnh-disrupted" if acc_base in profile.hnh_fragmented_accessions
               else "")
    out: list[Candidate] = []
    for idx, f in enumerate(record.features):
        if f.type != "CDS":
            continue
        text = qualifier_text(f)
        if not any(_kw_hit(text, kw) for kw in profile.lysis_keywords):
            continue
        translation = f.qualifiers.get("translation", [])
        if not translation or not translation[0]:
            continue
        product = f.qualifiers.get("product", ["hypothetical protein"])[0]
        pid = f.qualifiers.get("protein_id", [f"CDS_{idx:04d}"])[0]
        seq = translation[0]
        key = infer_evidence(text, profile)          # host-neutral detection
        label = domain_label(key, profile)           # host names it
        out.append(Candidate(
            protein_id=pid, accession=accession, organism=organism,
            product=product, sequence=seq, inferred_domain=label,
            domain_key=key,
            classification=classify_candidate(text, len(seq), key, profile,
                                              ident_text=identity_text(f)),
            runtime_flag=runtime,
        ))
    return out


def rank_score(c: Candidate, profile) -> int:
    """Higher = stronger free-endolysin; windows come from the profile."""
    text = (c.product + " " + c.inferred_domain).lower()
    score = 0
    if c.classification == CLASS_FREE:
        score += 1000
    if "endolysin" in text or "lysk" in text or "lysin" in text:
        score += 100
    if "chap" in text:
        score += 30
    if "amidase" in text:
        score += 30
    lo_c, hi_c = profile.endolysin_canonical_window_aa
    lo_p, hi_p = profile.endolysin_plausible_window_aa
    if lo_c <= c.length <= hi_c:
        score += 50
    elif lo_p <= c.length <= hi_p:
        score += 20
    if c.classification == CLASS_VAPH:
        score -= 500
    if c.classification == CLASS_NONLYTIC:
        score -= 1000
    return score


# ---------------------------------------------------------------------------
# tBLASTn fallback
# ---------------------------------------------------------------------------

def find_tblastn() -> Optional[str]:
    return shutil.which("tblastn")


def parse_outfmt6(text: str) -> Optional[dict]:
    best = None
    for line in text.strip().splitlines():
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        try:
            hit = dict(pident=float(cols[2]), length=int(cols[3]),
                       sstart=int(cols[4]), send=int(cols[5]),
                       evalue=float(cols[6]), bitscore=float(cols[7]))
        except ValueError:
            continue
        if best is None or hit["bitscore"] > best["bitscore"]:
            best = hit
    return best


def run_tblastn(record: SeqRecord, reference: Path,
                identity_threshold: float) -> tuple[Optional[float], str, str]:
    """tBLASTn LysK reference vs this genome; returns (identity, note, suffix)."""
    exe = find_tblastn()
    if exe is None:
        cmd = (f'tblastn -query "{reference}" -subject <genome.fasta> '
               f'-outfmt "6 qseqid sseqid pident length sstart send evalue bitscore"')
        return None, f"BLAST+ not on PATH — run manually: {cmd}", "tblastn-required"
    with tempfile.NamedTemporaryFile("w", suffix=".fna", delete=False,
                                     encoding="utf-8") as fh:
        fh.write(f">{record.id}\n{str(record.seq)}\n")
        subject = fh.name
    try:
        proc = subprocess.run(
            [exe, "-query", str(reference), "-subject", subject,
             "-outfmt", "6 qseqid sseqid pident length sstart send evalue bitscore",
             "-max_target_seqs", "5"],
            capture_output=True, text=True, timeout=300)
    except (subprocess.SubprocessError, OSError) as exc:
        return None, f"tblastn execution failed: {exc}", "tblastn-error"
    finally:
        Path(subject).unlink(missing_ok=True)
    hit = parse_outfmt6(proc.stdout)
    if hit is None:
        return None, "tblastn returned no hit", "tblastn-no-hit"
    coords = (f"c{hit['sstart']}-{hit['send']}" if hit["sstart"] > hit["send"]
              else f"{hit['sstart']}-{hit['send']}")
    if hit["pident"] >= identity_threshold:
        return hit["pident"], (f"tBLASTn {hit['pident']:.1f}% over {hit['length']} aa "
                               f"at {coords} — free endolysin recovered"), "tblastn-recovered"
    return hit["pident"], (f"tBLASTn only {hit['pident']:.1f}% over {hit['length']} aa "
                           f"at {coords} — divergent/uncertain"), "divergent-uncertain"


# ---------------------------------------------------------------------------
# Per-genome representative selection (curation registry first, then heuristic)
# ---------------------------------------------------------------------------

def select_for_genome(cands: list[Candidate], record: SeqRecord, profile,
                      reference: Optional[Path], run_tblastn_flag: bool,
                      identity_threshold: float
                      ) -> tuple[Optional[Candidate], list[str]]:
    """Choose one representative endolysin; resolve validated special cases."""
    msgs: list[str] = []
    acc = record.id.split(".")[0]

    if acc in profile.known_cases:
        kc = profile.known_cases[acc]
        status = normalise_class(kc["status"])
        pick = None
        if cands:
            if status == CLASS_VAPH:
                pick = max(cands, key=lambda c: c.length)
            elif status in ("intron-split", "divergent-endolysin"):
                pick = max(cands, key=lambda c: rank_score(c, profile))
            # tblastn-recovered / hnh-disrupted: annotation has only decoys/fragments
        if pick is not None:
            pick.classification = status
            pick.tblastn_identity = kc.get("identity")
            pick.tblastn_note = kc.get("note", "")
            pick.evidence = (f"matched '{pick.inferred_domain}'; validated status "
                             f"'{status}' (curation registry)")
            pick.selected = True
        ident_txt = (f"{kc['identity']:.0f}%" if kc.get("identity") is not None
                     else "n/a")
        msgs.append(f"[{status}] {record.id}: {kc.get('note','')} "
                    f"(tBLASTn identity: {ident_txt})")
        return pick, msgs

    # general heuristic
    free = [c for c in cands if c.classification == CLASS_FREE]
    if free:
        pick = max(free, key=lambda c: rank_score(c, profile))
        pick.selected = True
        pick.evidence = (f"matched '{pick.product[:30]}'; inferred "
                         f"{pick.inferred_domain}; free-endolysin "
                         f"(top of {len(cands)} candidate(s))")
        weak = ("nlpc" in pick.inferred_domain.lower()
                or not any(t in pick.inferred_domain.lower()
                           for t in ("chap", "amidase", "lysk")))
        if weak and run_tblastn_flag and reference is not None:
            ident, note, suffix = run_tblastn(record, reference, identity_threshold)
            pick.tblastn_identity, pick.tblastn_note = ident, note
            if suffix == "tblastn-recovered":
                msgs.append(f"[verify] {record.id}: weak keyword pick "
                            f"({pick.inferred_domain}); {note}")
        return pick, msgs

    vaph = [c for c in cands if c.classification == CLASS_VAPH]
    if run_tblastn_flag and reference is not None:
        ident, note, suffix = run_tblastn(record, reference, identity_threshold)
        msgs.append(f"[{suffix}] {record.id}: {note}")
        if vaph:
            v = max(vaph, key=lambda c: c.length)
            v.selected = True
            v.tblastn_identity, v.tblastn_note = ident, note
            v.evidence = f"VAPH ({v.inferred_domain}); no free endolysin; {note}"
            return v, msgs
        return None, msgs

    if vaph:
        v = max(vaph, key=lambda c: c.length)
        v.selected = True
        v.evidence = f"VAPH ({v.inferred_domain}); no free endolysin by keyword"
        msgs.append(f"[vaph] {record.id}: only a virion-associated hydrolase "
                    f"annotated; run --run-tblastn to seek a free endolysin.")
        return v, msgs

    msgs.append(f"[no-endolysin] {record.id}: no free endolysin in annotation; "
                f"run --run-tblastn with a LysK reference.")
    return None, msgs


# ===========================================================================
# Phase 2 module registry — stubs, switched on by profile.active_lysis_modules
# ===========================================================================

def _module_spanin(record: SeqRecord, profile) -> list[Candidate]:
    """Outer-membrane spanin detection for Gram-negative hosts. Phase 2."""
    raise NotImplementedError(
        "spanin module is a Phase-2 capability (Gram-negative hosts). "
        "Implement Rz/Rz1 (i-spanin/o-spanin) detection here and enable it via "
        "profile.active_lysis_modules.")


def _module_lysinB(record: SeqRecord, profile) -> list[Candidate]:
    """Mycolic-acid esterase (LysinB) detection for Mycobacterium. Phase 2."""
    raise NotImplementedError(
        "lysinB module is a Phase-2 capability (Mycobacterium). Implement "
        "mycolyl-arabinogalactan esterase detection here and enable it via "
        "profile.active_lysis_modules.")


# Modules handled by the core surfaces above (always available in Phase 1).
IMPLEMENTED_MODULES = frozenset({"endolysin", "holin"})

# Phase-2 modules: registered (so misconfiguration is distinguishable from a
# typo) but not yet implemented. Calling one, or activating one via a profile,
# fails loudly — a missing lysis mechanism is never silently skipped.
LYSIS_MODULE_REGISTRY: dict[str, Callable] = {
    "spanin": _module_spanin,
    "lysinB": _module_lysinB,
}


def assert_modules_available(profile) -> None:
    """
    Fail loudly if a profile requests a lysis module that is unknown (typo) or
    registered-but-not-yet-implemented (Phase-2 scaffold).
    """
    for m in profile.active_lysis_modules:
        if m in IMPLEMENTED_MODULES:
            continue
        if m not in LYSIS_MODULE_REGISTRY:
            raise ValueError(
                f"Profile '{profile.name}' requests unknown lysis module '{m}'. "
                f"Implemented: {sorted(IMPLEMENTED_MODULES)}; "
                f"registered Phase-2 scaffolds: {sorted(LYSIS_MODULE_REGISTRY)}.")
        raise NotImplementedError(
            f"Profile '{profile.name}' activates lysis module '{m}', which is a "
            f"Phase-2 scaffold not yet implemented. Implement _module_{m} in "
            f"phagecore/lysis.py before enabling it, or remove '{m}' from "
            f"active_lysis_modules.")
