#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phagecore.trna — universal tRNA de-novo analysis engine (S6 engine).

Moved out of the standalone S6 script (v3.0 S5-S7 refactor). The DETECTION and
parsing logic is host-agnostic and lives here ONCE so it cannot drift per host.
Only two host-specific values are read from the profile, never hardcoded:

  * profile.trna_ground_truth_acc      primary calibration anchor
  * profile.trna_secondary_anchor_acc  secondary anchor for concordance
  * profile.trna_canonical_isotypes    DERIVED signature (pattern classification)

Universal invariants that stay here unchanged:
  * CAT ambiguity: tRNAscan-SE cannot separate initiator/elongator Met from
    Ile-CAT; CAT rows are tagged cat_ambiguous so S7 handles Met/Ile2, not guesses.
  * pseudo / Undet(NNN) handling (excluded from the functional count).
The Ile2-CAT -> ATA (lysidine) decoding itself lives in phagecore.codon (S7),
consumed downstream; detection here only preserves the anticodon truthfully.

Two modes mirror the standalone script: PREPARE (nucleotide FASTA for the
tRNAscan-SE web server) and PARSE+ANALYSE (join the result table back per phage).
"""

from __future__ import annotations

import csv
import logging
import re
import sys
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

LOG = logging.getLogger("phagecore.trna")

# --------------------------------------------------------------------------- #
# Universal constants (NOT host-specific)
# --------------------------------------------------------------------------- #
COVE_DEFAULT_CUTOFF: float = 20.0
PSEUDO_TAG: str = "pseudo"
MANUAL_CHECK_TAGS: tuple[str, ...] = ("isotype mismatch", "unexpected anticodon")
CAT_AMBIGUOUS_ANTICODON: str = "CAT"          # Met / Ile-CAT not distinguished
GENBANK_SUFFIXES: tuple[str, ...] = (".gb", ".gbk", ".gbff", ".genbank", ".gbf")


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #
def clean_accession(record_id: str) -> str:
    """Normalise a record id to accession.version; unwrap 'gi|..|gb|X|'."""
    rid = (record_id or "").strip()
    if "|" not in rid:
        return rid
    parts = rid.split("|")
    for tag in ("gb", "ref", "emb", "dbj", "tpg", "tpe", "tpd"):
        if tag in parts:
            i = parts.index(tag)
            if i + 1 < len(parts) and parts[i + 1]:
                return parts[i + 1]
    toks = [p for p in parts if p]
    return toks[-1] if toks else rid


def safe_seqname(acc: str) -> str:
    """tRNAscan-SE web server accepts only [A-Za-z0-9_:-]; dot -> '_'."""
    return re.sub(r"[^A-Za-z0-9_:-]", "_", acc)


def _acc_base(acc: str) -> str:
    return acc.split(".")[0]


def iter_genbank_files(input_dir: Path):
    for path in sorted(input_dir.iterdir()):
        if path.suffix.lower() in GENBANK_SUFFIXES:
            yield path


def load_records(input_dir: Path) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    for path in iter_genbank_files(input_dir):
        try:
            for rec in SeqIO.parse(str(path), "genbank"):
                records.append(rec)
        except Exception as exc:
            LOG.warning("Could not parse %s: %s", path.name, exc)
    return records


# --------------------------------------------------------------------------- #
# Mode 1 — prepare nucleotide FASTA
# --------------------------------------------------------------------------- #
def _write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for acc, seq in entries:
            fh.write(f">{safe_seqname(acc)}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def write_nucleotide_fasta(input_dir: Path, out_path: Path,
                           chunk_mbp: Optional[float]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    records = load_records(input_dir)
    if not records:
        LOG.error("No GenBank files found in '%s'.", input_dir)
        sys.exit(1)

    entries: list[tuple[str, str]] = []
    skipped = 0
    for rec in records:
        seq = str(rec.seq)
        if not seq or set(seq.upper()) <= {"N"}:
            skipped += 1
            LOG.warning("  %s: no usable nucleotide sequence — skipped.",
                        clean_accession(rec.id))
            continue
        entries.append((clean_accession(rec.id), seq))

    total_bp = sum(len(s) for _, s in entries)
    LOG.info("Genomes with sequence: %d (skipped %d). Total %.2f Mbp.",
             len(entries), skipped, total_bp / 1e6)

    if chunk_mbp is None:
        _write_fasta(out_path, entries)
        LOG.info("Wrote single nucleotide FASTA: '%s' (%d genomes).",
                 out_path, len(entries))
        LOG.info("  Submit to https://trna.ucsc.edu/tRNAscan-SE/ "
                 "(Sequence source: Bacterial), save TABULAR output.")
        if total_bp / 1e6 > 50:
            LOG.warning("  %.1f Mbp is large; if rejected, re-run with --chunk-mbp 5.",
                        total_bp / 1e6)
        return

    limit = chunk_mbp * 1e6
    chunks: list[list[tuple[str, str]]] = [[]]
    running = 0.0
    for acc, seq in entries:
        if running + len(seq) > limit and chunks[-1]:
            chunks.append([])
            running = 0.0
        chunks[-1].append((acc, seq))
        running += len(seq)
    stem, suf = out_path.with_suffix(""), out_path.suffix or ".fasta"
    for i, chunk in enumerate(chunks, 1):
        part = Path(f"{stem}_{i}{suf}")
        _write_fasta(part, chunk)
        LOG.info("  chunk %d/%d: '%s' (%d genomes, %.2f Mbp)",
                 i, len(chunks), part.name, len(chunk),
                 sum(len(s) for _, s in chunk) / 1e6)
    LOG.info("Submit each chunk, then concatenate result tables before Mode 2.")


# --------------------------------------------------------------------------- #
# Mode 2 — parse tRNAscan-SE tabular output
# --------------------------------------------------------------------------- #
class TrnaHit:
    __slots__ = ("seqname", "index", "begin", "end", "aa", "anticodon",
                 "intron_begin", "intron_end", "score", "note")

    def __init__(self, seqname, index, begin, end, aa, anticodon,
                 intron_begin, intron_end, score, note):
        self.seqname = seqname
        self.index = index
        self.begin = begin
        self.end = end
        self.aa = aa
        self.anticodon = anticodon
        self.intron_begin = intron_begin
        self.intron_end = intron_end
        self.score = score
        self.note = note

    @property
    def strand(self) -> str:
        return "-" if self.begin > self.end else "+"

    @property
    def start_bp(self) -> int:
        return min(self.begin, self.end)

    @property
    def end_bp(self) -> int:
        return max(self.begin, self.end)

    @property
    def is_pseudo(self) -> bool:
        return PSEUDO_TAG in self.note.lower()

    @property
    def cat_ambiguous(self) -> bool:
        return self.anticodon.upper() == CAT_AMBIGUOUS_ANTICODON

    @property
    def manual_check(self) -> str:
        low = self.note.lower()
        hits = [t for t in MANUAL_CHECK_TAGS if t in low]
        return ";".join(hits)


_DATA_ROW = re.compile(r"^\S+\s+\d+\s+\d+\s+\d+\s+\S+")


def parse_trnascan(path: Path) -> dict[str, list[TrnaHit]]:
    """Parse a tRNAscan-SE tabular result file -> {accession: [TrnaHit, ...]}."""
    by_acc: dict[str, list[TrnaHit]] = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or not _DATA_ROW.match(line):
                continue
            cols = line.split("\t") if "\t" in line else line.split()
            cols = [c.strip() for c in cols if c.strip() != ""]
            if len(cols) < 8:
                continue
            try:
                seqname = cols[0].strip()
                index = int(cols[1])
                begin = int(cols[2])
                end = int(cols[3])
                aa = cols[4]
                anticodon = cols[5].upper()
                intron_begin = int(cols[6])
                intron_end = int(cols[7])
                score = float(cols[8]) if len(cols) > 8 else float("nan")
                note = " ".join(cols[9:]) if len(cols) > 9 else ""
            except (ValueError, IndexError):
                continue
            by_acc.setdefault(seqname, []).append(
                TrnaHit(seqname, index, begin, end, aa, anticodon,
                        intron_begin, intron_end, score, note))
    return by_acc


# --------------------------------------------------------------------------- #
# Metadata + pattern + ground-truth (host anchors/signature from the profile)
# --------------------------------------------------------------------------- #
def genome_metadata(records: list[SeqRecord]) -> dict[str, dict]:
    meta: dict[str, dict] = {}
    for rec in records:
        acc = clean_accession(rec.id)
        from phagecore.genbank_io import resolve_organism
        org = resolve_organism(rec, getattr(rec, "_source_path", None))
        meta[acc] = {"organism": org, "size": len(rec.seq) or 0}
    return meta


def classify_pattern(aas: set[str], canonical: frozenset[str]) -> str:
    """Pattern relative to the host's DERIVED canonical isotype set (from profile)."""
    if not aas:
        return "none"
    if not canonical:
        return "uncalibrated"          # profile did not supply a signature
    if aas == canonical:
        return "canonical"
    if canonical.issubset(aas):
        return "canonical+"
    if aas & canonical:
        return "partial-canonical"
    return "noncanonical"


def annotated_trna_from_genbank(records: list[SeqRecord], acc: str) -> list[tuple[str, str]]:
    """Extract (aa, anticodon) from GenBank tRNA features. anticodon may be ''
    when the depositor annotated isotype only (common for B. subtilis phages)."""
    out: list[tuple[str, str]] = []
    for rec in records:
        if _acc_base(clean_accession(rec.id)) != _acc_base(acc):
            continue
        for feat in rec.features:
            if feat.type != "tRNA":
                continue
            product = " ".join(feat.qualifiers.get("product", [""]))
            m_aa = re.search(r"tRNA-([A-Za-z]{3})", product)
            aa = m_aa.group(1) if m_aa else "Undet"
            anti = ""
            quals = (feat.qualifiers.get("anticodon", [])
                     + feat.qualifiers.get("note", []))
            for q in quals:
                m_aa2 = re.search(r"aa:([A-Za-z]{3})", q)
                if m_aa2:
                    aa = m_aa2.group(1)
                m = re.search(r"seq:\(?([ACGTUacgtu]{3})\)?", q)
                if not m:
                    m = re.search(r"anticodon[:= ]*\(?([ACGTUacgtu]{3})\)?", q)
                if m:
                    anti = m.group(1).upper().replace("U", "T")
                    break
            out.append((aa, anti))
    return out


def _denovo_isotypes(by_acc, records, acc) -> tuple[Optional[str], list, list]:
    """Resolve the full accession, its GenBank (aa,anti), and de-novo hits."""
    full = next((clean_accession(r.id) for r in records
                 if _acc_base(clean_accession(r.id)) == _acc_base(acc)), None)
    if full is None:
        return None, [], []
    by_safe: dict[str, list] = {}
    for k, v in by_acc.items():
        by_safe.setdefault(safe_seqname(k), []).extend(v)
    dn = [h for h in by_safe.get(safe_seqname(full), []) if not h.is_pseudo]
    ann = annotated_trna_from_genbank(records, full)
    return full, ann, dn


def ground_truth_check(records, by_acc, ground_truth_acc: str,
                       secondary_acc: str = "") -> Optional[dict]:
    """
    Calibrate the de-novo pipeline against the profile's anchor(s).

    ANTICODON-level when GenBank supplies anticodons (e.g. Staph Sb-1). ISOTYPE-
    level fallback when GenBank supplies isotype only (B. subtilis phages, which
    lack GenBank anticodons) — Decision D2. A CAT anticodon shared counts as
    agreement (Met vs Ile2-CAT is the same gene by design). Reports the honest
    recovery fraction and a PASS/REVIEW verdict; adds primary-vs-secondary de-novo
    isotype concordance when a secondary anchor is configured.

    Returns a dict {line, verdict, level, recovered, denom, missing} or None. The
    caller uses `recovered < denom` to mark Pattern PROVISIONAL: classifying every
    phage against a canonical set while the anchor itself under-recovers de-novo is
    a comparison artefact, not biology (e.g. Met is invisible to tRNAscan-SE here).
    """
    if not ground_truth_acc:
        return None
    full, ann, dn = _denovo_isotypes(by_acc, records, ground_truth_acc)
    if full is None:
        return None

    ann_aa = sorted(a for a, _ in ann)
    dn_typed = {h.aa for h in dn if not h.aa.lower().startswith("undet")
                and h.anticodon.upper() != "NNN"}
    dn_undet = sum(1 for h in dn if h.aa.lower().startswith("undet")
                   or h.anticodon.upper() == "NNN")
    ann_anti = {c for _, c in ann if c}

    if ann_anti:  # ---- anticodon-level (GenBank has anticodons) ----
        dn_anti = {h.anticodon for h in dn}
        shared = ann_anti & dn_anti
        missing = ann_anti - dn_anti
        verdict = "PASS" if shared == ann_anti else "REVIEW"
        level, denom = "anticodon", len(ann_anti)
        recov = len(shared)
        cat_note = " [CAT gene present; Met/Ile2 differ by design]" if "CAT" in shared else ""
        detail = (f"anticodon agreement {recov}/{denom}"
                  + (f"; MISSING de-novo {sorted(missing)}" if missing else "")
                  + cat_note)
    else:         # ---- isotype-level fallback (Decision D2) ----
        ann_iso = set(ann_aa)
        shared = ann_iso & dn_typed
        missing = ann_iso - dn_typed
        verdict = "PASS" if shared == ann_iso else "REVIEW"
        level, denom = "isotype", len(ann_iso)
        recov = len(shared)
        undet_note = f"; {dn_undet} de-novo Undet(NNN)" if dn_undet else ""
        detail = (f"isotype recovery {recov}/{denom}"
                  + (f"; MISSING de-novo {sorted(missing)}" if missing else "")
                  + undet_note
                  + " [GenBank has no anticodons -> ARAGORN concordance pending]")

    line = (f"[{verdict}] primary {full} ({level}-level): "
            f"GenBank {len(ann)} tRNA {ann_aa}; de-novo {len(dn)} tRNA "
            f"{sorted(h.aa for h in dn)}; {detail}")

    if secondary_acc:
        s_full, _s_ann, s_dn = _denovo_isotypes(by_acc, records, secondary_acc)
        if s_full:
            s_typed = {h.aa for h in s_dn if not h.aa.lower().startswith("undet")
                       and h.anticodon.upper() != "NNN"}
            concord = sorted(dn_typed & s_typed)
            line += (f" | secondary {s_full} de-novo isotypes {sorted(s_typed)}; "
                     f"primary∩secondary (de-novo) = {concord}")
        else:
            line += f" | secondary {secondary_acc} not in set — concordance skipped"
    return {"line": line, "verdict": verdict, "level": level,
            "recovered": recov, "denom": denom, "missing": sorted(missing)}


# --------------------------------------------------------------------------- #
# Mode 2 driver
# --------------------------------------------------------------------------- #
DETAIL_FIELDS = ["Phage", "Accession", "Genome_Size_bp", "tRNA_Index", "tRNA_Type",
                 "Anticodon_Seq", "Cove_Score", "Start_bp", "End_bp", "Strand",
                 "Relative_Position_pct", "Is_Pseudo", "CAT_Ambiguous",
                 "Manual_Check", "Note"]

SUMMARY_BASE_FIELDS = ["Phage", "Accession", "Genome_Size_bp", "tRNA_Count_functional",
                       "tRNA_Count_pseudo", "tRNA_Types", "Anticodons", "Pattern"]
SUMMARY_TAIL_FIELDS = ["Pattern_Status", "CAT_Ambiguous_n", "Flags"]


def analyse(input_dir: Path, results_path: Path,
            detail_out: Path, summary_out: Path, profile) -> None:
    """Parse tRNAscan-SE output, join per phage, write detail + summary CSVs.

    Host layer: canonical signature and anchors come from `profile`. Detection is
    universal. The summary carries one dynamic Has_<isotype> column per canonical
    isotype declared by the profile (Staph: Met/Trp/Phe/Asp; B. subtilis: the
    BSP38-derived set) instead of a hardcoded Has_Met/Trp/Phe/Asp block.
    """
    canonical = frozenset(profile.trna_canonical_isotypes)
    has_cols = [f"Has_{iso}" for iso in profile.trna_canonical_isotypes]
    summary_fields = SUMMARY_BASE_FIELDS + has_cols + SUMMARY_TAIL_FIELDS

    records = load_records(input_dir)
    if not records:
        LOG.error("No GenBank files found in '%s' (needed for metadata).", input_dir)
        sys.exit(1)
    meta = genome_metadata(records)
    by_acc = parse_trnascan(results_path)
    by_safe: dict[str, list] = {}
    for k, v in by_acc.items():
        by_safe.setdefault(safe_seqname(k), []).extend(v)
    LOG.info("Parsed tRNAscan-SE: %d tRNA across %d sequences.",
             sum(len(v) for v in by_acc.values()), len(by_acc))

    # Anchor recovery decides whether Pattern is trustworthy. If the ground-truth
    # anchor itself under-recovers de-novo (e.g. Met invisible to tRNAscan-SE),
    # classifying every phage against a GenBank-derived canonical set is a
    # comparison artefact, not biology -> Pattern is PROVISIONAL, and Has_<iso> +
    # tRNA_Types (raw de-novo, artefact-free) lead the analysis until ARAGORN
    # concordance lifts recovery.
    gt = ground_truth_check(records, by_acc, profile.trna_ground_truth_acc,
                            profile.trna_secondary_anchor_acc)
    if gt and gt["recovered"] < gt["denom"]:
        pattern_status = (f"provisional: de-novo under-recovers anchor "
                          f"{gt['denom'] - gt['recovered']}/{gt['denom']} "
                          f"(missing {gt['missing']}); pending ARAGORN concordance")
        LOG.warning("Pattern is PROVISIONAL — anchor recovery %d/%d de-novo "
                    "(missing %s). Lead with Has_<iso> + tRNA_Types; do not headline "
                    "Pattern until ARAGORN lifts recovery.",
                    gt["recovered"], gt["denom"], gt["missing"])
    elif gt:
        pattern_status = "validated (anchor fully recovered de-novo)"
    else:
        pattern_status = "uncalibrated (anchor not in set)"

    detail_out.parent.mkdir(parents=True, exist_ok=True)
    detail_rows: list[dict] = []
    summary_rows: list[dict] = []

    for acc in sorted(meta):
        m = meta[acc]
        org, size = m["organism"], m["size"]
        hits = by_safe.get(safe_seqname(acc), [])
        pseudo = [h for h in hits if h.is_pseudo]
        undet = [h for h in hits if not h.is_pseudo and
                 (h.aa.lower().startswith("undet") or h.anticodon.upper() == "NNN")]
        functional = [h for h in hits if not h.is_pseudo and h not in undet]
        for h in sorted(hits, key=lambda x: x.index):
            rel = round(100 * h.start_bp / size, 1) if size else ""
            detail_rows.append({
                "Phage": org, "Accession": acc, "Genome_Size_bp": size,
                "tRNA_Index": h.index, "tRNA_Type": h.aa,
                "Anticodon_Seq": h.anticodon, "Cove_Score": h.score,
                "Start_bp": h.start_bp, "End_bp": h.end_bp, "Strand": h.strand,
                "Relative_Position_pct": rel,
                "Is_Pseudo": "Yes" if h.is_pseudo else "No",
                "CAT_Ambiguous": "Yes" if h.cat_ambiguous else "No",
                "Manual_Check": h.manual_check, "Note": h.note})
        aas = {h.aa for h in functional}
        flags = []
        if any(h.manual_check for h in hits):
            flags.append("manual_check_rows")
        if pseudo:
            flags.append(f"{len(pseudo)}_pseudo")
        if undet:
            flags.append(f"{len(undet)}_undetermined_isotype")
        row = {
            "Phage": org, "Accession": acc, "Genome_Size_bp": size,
            "tRNA_Count_functional": len(functional),
            "tRNA_Count_pseudo": len(pseudo),
            "tRNA_Types": ";".join(sorted(aas)) if aas else "—",
            "Anticodons": ";".join(sorted({h.anticodon for h in functional})) or "—",
            "Pattern": classify_pattern(aas, canonical),
            "Pattern_Status": pattern_status,
            "CAT_Ambiguous_n": sum(1 for h in functional if h.cat_ambiguous),
            "Flags": ";".join(flags)}
        for iso in profile.trna_canonical_isotypes:
            row[f"Has_{iso}"] = "Yes" if iso in aas else "No"
        summary_rows.append(row)

    with open(detail_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=DETAIL_FIELDS)
        w.writeheader()
        w.writerows(detail_rows)
    with open(summary_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)

    n_with = sum(1 for r in summary_rows if r["tRNA_Count_functional"] > 0)
    LOG.info("Detail CSV : '%s' (%d tRNA rows)", detail_out, len(detail_rows))
    LOG.info("Summary CSV: '%s' (%d phages, %d carry >=1 functional tRNA)",
             summary_out, len(summary_rows), n_with)
    if gt:
        LOG.info("CALIBRATION: %s", gt["line"])
    else:
        LOG.info("CALIBRATION: ground-truth anchor '%s' not in set — check skipped.",
                 profile.trna_ground_truth_acc or "(none configured)")
