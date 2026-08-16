"""
phagecore.interpro
=================
Reconcile S4's keyword-derived endolysin set against InterPro domain calls.

Keyword scanning has high recall but imperfect precision: broad terms like
"hydrolase" admit non-PG enzymes (dUTPase, RNase, NTPase), and proteins whose
endolysin annotation is only a note can pull in tail spikes or chaperonins.
InterPro is the authoritative arbiter. This module reads an InterPro result,
classifies each sequence's domains, and produces a confirmed-endolysin verdict.

Two input formats are auto-detected:
  1. InterProScan TSV (the reproducible EBI download) — tab-separated; column 1
     is the sequence ID (the FASTA header, accession-first), and a later column
     carries the signature/entry description.
  2. The "header = domains" summary a user may paste from the InterPro web
     entries panel (best effort; joined by accession when present, else by
     normalised organism name).

Verdicts (canonical):
  confirmed-endolysin  domains include a PG-hydrolase / CBD signature
  non-endolysin        domains match the profile denylist (dUTPase, GroEL, ...)
  no-data              InterPro returned nothing (hypothetical — needs structure)
  uncertain            domains present but neither allow- nor deny-listed
"""

from __future__ import annotations

import re
from pathlib import Path

# PG-hydrolase / endolysin domain signatures (allowlist).
_PG_DOMAIN_TERMS = (
    "chap", "amidase", "muramoyl", "sh3", "glucosaminidase", "muramidase",
    "lysozyme", "peptidoglycan binding", "pgbd", "glycosyl hydrolase family 73",
    "mannosyl-glycoprotein", "phage tail lysozyme", "nlpc", "p60",
    "transglycosylase", "peptidase_m15", "cwlq", "lytic transglycosylase",
)

# Non-endolysin domain signatures (denylist) — mirror of the profile denylist
# but expressed as InterPro entry-name fragments.
_NON_ENDO_DOMAIN_TERMS = (
    "dutpase", "nucleotidohydrolase", "pyrophosphat", "chaperonin", "groel",
    "tcp-1", "cpn60", "ribonuclease", "hydroxyacylglutathione", "p-loop",
    "nucleoside triphosphate hydrolase", "tail spike", "tailspike",
    "helicase", "primase", "polymerase", "terminase", "dut-like",
)


def classify_domains(domain_str: str) -> str:
    """Return a verdict for one sequence's InterPro domain string."""
    d = (domain_str or "").lower().strip()
    if not d or "no data available" in d or d in ("-", "none"):
        return "no-data"
    if any(t in d for t in _NON_ENDO_DOMAIN_TERMS):
        return "non-endolysin"
    if any(t in d for t in _PG_DOMAIN_TERMS):
        return "confirmed-endolysin"
    return "uncertain"


def _norm_org(name: str) -> str:
    """Normalise an organism / phage name to a join key."""
    n = name.lower()
    n = re.sub(r"staphylococcus[_ ]phage[_ ]", "", n)
    n = re.sub(r"[^a-z0-9]", "", n)
    return n


_ACC_RE = re.compile(r"\b([A-Z]{1,2}_?\d{5,8}\.\d+|[A-Z]{2}\d{6}\.\d+)\b")


def _extract_accession(text: str) -> str:
    m = _ACC_RE.search(text)
    return m.group(1) if m else ""


def parse_interpro(path: Path) -> dict:
    """
    Parse an InterPro result file into {key: {"domains", "verdict", "raw"}}.

    Keys are accession (preferred) and normalised organism (fallback); both are
    registered when available so reconciliation can match on either.

    `path` may be a single file, a DIRECTORY (every *.tsv/*.txt inside is merged),
    or a comma-separated list of files — InterPro splits a submission of >100
    sequences into several result batches, and all of them must be reconciled in
    one pass. Matches the behaviour of the S5 reconciler.
    """
    if isinstance(path, Path) and path.is_dir():
        files = sorted(f for f in path.iterdir()
                       if f.suffix.lower() in (".tsv", ".txt"))
        if not files:
            raise FileNotFoundError(f"no *.tsv / *.txt inside {path}")
    elif "," in str(path):
        files = [Path(p.strip()) for p in str(path).split(",") if p.strip()]
    else:
        files = [Path(path)]

    text = "\n".join(f.read_text(encoding="utf-8", errors="ignore") for f in files)
    lines = [l.rstrip("\n") for l in text.splitlines() if l.strip()]
    out: dict = {}

    # Detect TSV (InterProScan standard) vs the pasted "header = domains" format.
    is_tsv = any("\t" in l and " = " not in l for l in lines[:5])

    if is_tsv:
        # InterProScan TSV: col0 = seq id; collect signature/entry descriptions.
        agg: dict = {}
        for l in lines:
            cols = l.split("\t")
            if len(cols) < 6:
                continue
            seqid = cols[0]
            # entry description is typically the last informative column;
            # gather columns 5 and 12 (signature desc / InterPro entry name).
            desc_bits = [c for i, c in enumerate(cols)
                         if i in (5, 11, 12) and c and c != "-"]
            agg.setdefault(seqid, []).extend(desc_bits)
        for seqid, bits in agg.items():
            domains = ", ".join(dict.fromkeys(bits))   # dedup, keep order
            rec = {"domains": domains, "verdict": classify_domains(domains),
                   "raw": seqid}
            acc = _extract_accession(seqid)
            if acc:
                out[acc] = rec
                out[acc.split(".")[0]] = rec
            out.setdefault(_norm_org(seqid), rec)
        return out

    # Pasted "header = domains" format.
    for l in lines:
        if " = " not in l or l.startswith("Ada "):
            continue
        head, _, domains = l.partition(" = ")
        verdict = classify_domains(domains)
        rec = {"domains": domains.strip(), "verdict": verdict, "raw": head.strip()}
        acc = _extract_accession(head)
        if acc:
            out[acc] = rec
            out[acc.split(".")[0]] = rec
        m = re.search(r"phage[_ ](.+?)\|status=", head)
        if m:
            out.setdefault(_norm_org("phage_" + m.group(1)), rec)
    return out


def reconcile(selected: list[dict], interpro_map: dict) -> tuple[list[dict], dict]:
    """
    Attach an InterPro verdict + domains to each selected representative.

    `selected` is a list of dicts with at least 'accession' and 'organism'.
    Returns (annotated_selected, summary). Matching tries accession (full then
    base) then normalised organism. Unmatched reps get verdict 'not-in-interpro'.
    """
    summary = {"confirmed-endolysin": 0, "non-endolysin": 0, "no-data": 0,
               "uncertain": 0, "not-in-interpro": 0}
    annotated = []
    for p in selected:
        acc = p.get("accession", "")
        rec = (interpro_map.get(acc)
               or interpro_map.get(acc.split(".")[0])
               or interpro_map.get(_norm_org(p.get("organism", ""))))
        if rec is None:
            verdict, domains = "not-in-interpro", ""
        else:
            verdict, domains = rec["verdict"], rec["domains"]
        q = dict(p)
        q["interpro_verdict"] = verdict
        q["interpro_domains"] = domains
        annotated.append(q)
        summary[verdict] = summary.get(verdict, 0) + 1
    return annotated, summary
