"""
phagecore.qc
============
Quality-control gates, outlier flagging, and sequence-level deduplication.

This is the layer that was MISSING from the original S1/S2/S4. Its job is to
turn silent corruption into explicit, auditable flags:

  * CDS=0 on a genome large enough to require genes  -> qc_status FAIL_NO_CDS
  * genome size / GC outside profile bounds          -> qc_flags size/gc outlier
  * identical nucleotide sequences across accessions -> duplicate group + one
    representative chosen by a documented, deterministic rule (provenance kept,
    nothing deleted)

Nothing here removes data. It annotates. The researcher decides what to exclude,
which matches a conservative, defensible workflow (and lets the same run serve
both "report prevalence over all genomes" and "analyse unique sequences only").
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from .genbank_io import ParsedGenome

log = logging.getLogger("phagecore.qc")


@dataclass
class QCResult:
    """QC verdict for one genome."""
    qc_status: str = "PASS"          # PASS | FAIL_NO_CDS | WARN
    qc_flags: list[str] = field(default_factory=list)
    annotation_present: bool = True  # False when only 'source' (no gene/CDS) exists
    duplicate_group: str = ""        # md5-derived group id when in a dup cluster
    is_representative: bool = True    # one True per duplicate group
    note: str = ""


def evaluate_qc(pg: ParsedGenome, seq_len: int, n_cds: int, gc: float,
                profile) -> QCResult:
    """
    Per-genome QC against the host profile's plausibility bounds.

    `profile` supplies: min_cds_genome_size_bp, size_bounds_bp, gc_bounds_pct.
    """
    r = QCResult()
    census = pg.feature_census
    non_source = {k: v for k, v in census.items() if k != "source"}
    r.annotation_present = bool(non_source)

    # --- biologically impossible: large genome, zero CDS ---
    if n_cds == 0 and seq_len >= profile.min_cds_genome_size_bp:
        r.qc_status = "FAIL_NO_CDS"
        if not r.annotation_present:
            r.note = ("annotation absent (only 'source' feature) — re-download a "
                      "RefSeq/annotated version of this accession")
        else:
            present = ", ".join(f"{k}={v}" for k, v in sorted(non_source.items()))
            r.note = (f"0 CDS but other features present ({present}) — CDS likely "
                      f"stored under a non-standard structure; re-download RefSeq "
                      f"or re-annotate (e.g. Pharokka)")
        r.qc_flags.append("no_cds")

    # --- size / GC outlier flags (warn, do not fail) ---
    lo_bp, hi_bp = profile.size_bounds_bp
    if seq_len < lo_bp:
        r.qc_flags.append("size_below_expected")
    elif seq_len > hi_bp:
        r.qc_flags.append("size_above_expected")    # e.g. jumbo Twort-like phages

    lo_gc, hi_gc = profile.gc_bounds_pct
    if gc < lo_gc:
        r.qc_flags.append("gc_below_expected")
    elif gc > hi_gc:
        r.qc_flags.append("gc_above_expected")      # e.g. FV3 at 43.5% — verify host

    if r.qc_status == "PASS" and r.qc_flags:
        r.qc_status = "WARN"
    return r


def assign_duplicates(genomes: list[ParsedGenome],
                      prefer_refseq: bool = True) -> dict[str, QCResult]:
    """
    Group genomes with identical nucleotide sequence (seq_md5) and choose ONE
    representative per group.

    Returns {accession: QCResult-with-dup-fields-filled}. Caller merges these
    into the per-genome QCResult from evaluate_qc.

    Representative rule (deterministic, documented in SUPPLEMENTARY):
      prefer_refseq=True  -> keep the RefSeq (NC_) copy; uniform NCBI annotation
                             maximises keyword recall in S2/S4. This matches the
                             22-phage validated set, which is RefSeq-dominant.
      prefer_refseq=False -> keep the primary INSDC deposit (original, citable).
      Ties broken by accession string (stable, reproducible).
    """
    by_md5: dict[str, list[ParsedGenome]] = {}
    for g in genomes:
        by_md5.setdefault(g.seq_md5, []).append(g)

    out: dict[str, QCResult] = {}
    for md5, members in by_md5.items():
        if len(members) == 1:
            out[members[0].accession] = QCResult(duplicate_group="",
                                                 is_representative=True)
            continue
        group_id = md5[:10]

        def sort_key(g: ParsedGenome):
            # primary key: provenance preference; secondary: accession string
            pref = (0 if g.is_refseq == prefer_refseq else 1)
            return (pref, g.accession)

        ranked = sorted(members, key=sort_key)
        rep = ranked[0]
        others = ", ".join(g.accession for g in ranked[1:])
        for g in members:
            is_rep = (g.accession == rep.accession)
            out[g.accession] = QCResult(
                qc_flags=["duplicate_sequence"] if not is_rep else [],
                duplicate_group=group_id,
                is_representative=is_rep,
                note=(f"identical genome; representative={rep.accession}"
                      if not is_rep else
                      f"representative of duplicate group {group_id} "
                      f"(identical: {others})"),
            )
    return out


def merge_dup_into_qc(qc: QCResult, dup: QCResult) -> QCResult:
    """Fold duplicate-group fields from assign_duplicates into a per-genome QC."""
    qc.duplicate_group = dup.duplicate_group
    qc.is_representative = dup.is_representative
    for fl in dup.qc_flags:
        if fl not in qc.qc_flags:
            qc.qc_flags.append(fl)
    if dup.note:
        qc.note = (qc.note + " | " + dup.note).strip(" |")
    if not dup.is_representative and qc.qc_status == "PASS":
        qc.qc_status = "WARN"
    return qc
