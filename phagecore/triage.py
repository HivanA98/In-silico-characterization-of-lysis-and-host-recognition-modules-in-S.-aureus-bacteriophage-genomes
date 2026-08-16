#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phagecore.triage — Pharokka re-annotation triage (S4 output).

Flags source genomes whose GenBank annotation is a GAP rather than biology, so the
operator re-annotates them with Pharokka BEFORE trusting downstream feature-table
extraction (S2/S4/S5/S6). Encodes the operator's rule:

  A genome qualifies for Pharokka re-annotation if >=1 TRIGGER fires and no
  exclusion applies.

TRIGGERS
  T1  CDS = 0 on a genome >= 10 kb  -> biologically impossible; annotation absent
      (this is the S1 FAIL_NO_CDS case).
  T2  Genome has CDS + a structural module but ZERO lysis-keyword CDS -> tailed
      phages always carry endolysin/holin; absence is an annotation gap, not
      biology (the PBS1 / MW354668 cases).
  T3  Only accessory fragments annotated (CDS present, but NO structural module
      AND no lysis module - e.g. only dUTPase) -> partial annotation (OL580764).
  T4  No CDS features at all on a <10 kb record -> downloaded as FASTA / novel
      isolate with no annotation to preserve.

EXCLUSIONS (no re-annotation even if a trigger fires)
  X1  Already Pharokka-annotated (source qualifier / structured_comment mentions
      Pharokka) -> re-annotation would be redundant.
  (Add taxon-based exclusions here if a non-tailed lineage legitimately lacks a
  lysis module; not applied by default for Caudoviricetes.)

Genotype != phenotype is unaffected: this only judges ANNOTATION completeness, not
biological capacity.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from Bio import SeqIO

log = logging.getLogger("phagecore.triage")

# Universal phage structural-module vocabulary (host-agnostic).
STRUCTURAL_TERMS = (
    "tail", "capsid", "portal", "terminase", "baseplate", "base plate", "head",
    "sheath", "tail tube", "major capsid", "virion", "tail fiber", "tail fibre",
    "tail spike", "tailspike", "neck", "collar", "prohead", "scaffold",
    "connector", "spike", "receptor",
)

MIN_LARGE_GENOME_BP = 10_000
TRIAGE_FIELDS = ["Accession", "Organism", "Genome_Size_bp", "CDS_Count",
                 "Has_Structural", "Has_Lysis", "Needs_Pharokka",
                 "Trigger", "Detail"]


def _org(record):
    from phagecore.genbank_io import resolve_organism
    return resolve_organism(record, getattr(record, "_source_path", None))


def _pharokka_already(record) -> bool:
    """Best-effort: was this record already annotated by Pharokka?"""
    sc = " ".join(str(v) for v in record.annotations.get("structured_comment", {}).values()) \
         if isinstance(record.annotations.get("structured_comment"), dict) else ""
    blob = (sc + " " + " ".join(record.annotations.get("comment", "") if isinstance(
        record.annotations.get("comment"), str) else [])).lower()
    for feat in record.features:
        if feat.type == "source":
            blob += " " + " ".join(feat.qualifiers.get("note", [])).lower()
    return "pharokka" in blob


def triage_record(record, profile) -> dict:
    """Evaluate one record against the Pharokka triggers/exclusions."""
    size = len(record.seq)
    cds = [f for f in record.features if f.type == "CDS"]
    n_cds = len(cds)
    lysis_terms = tuple(t.lower() for t in profile.lysis_keywords)

    def _prod(f):
        return " ".join(f.qualifiers.get("product", [])).lower()

    has_lysis = any(any(k in _prod(f) for k in lysis_terms) for f in cds)
    has_structural = any(any(k in _prod(f) for k in STRUCTURAL_TERMS) for f in cds)

    triggers: list[str] = []
    detail: list[str] = []

    if size >= MIN_LARGE_GENOME_BP and n_cds == 0:
        triggers.append("T1_no_CDS_large_genome")
        detail.append(f"{size:,} bp but 0 CDS (biologically impossible)")
    if n_cds == 0 and size < MIN_LARGE_GENOME_BP:
        triggers.append("T4_no_annotation")
        detail.append("no CDS features (FASTA / novel isolate, nothing to preserve)")
    if n_cds > 0 and has_structural and not has_lysis:
        triggers.append("T2_no_lysis_gene")
        detail.append("structural module present but 0 lysis-keyword CDS (annotation gap)")
    if n_cds > 0 and not has_structural and not has_lysis:
        triggers.append("T3_only_accessory")
        detail.append("CDS present but no structural or lysis module (accessory-only)")

    excluded = _pharokka_already(record)
    needs = bool(triggers) and not excluded
    if excluded and triggers:
        detail.append("EXCLUDED: already Pharokka-annotated")

    return {
        "Accession": record.id,
        "Organism": _org(record),
        "Genome_Size_bp": size,
        "CDS_Count": n_cds,
        "Has_Structural": "Yes" if has_structural else "No",
        "Has_Lysis": "Yes" if has_lysis else "No",
        "Needs_Pharokka": "Yes" if needs else "No",
        "Trigger": ";".join(triggers) if triggers else "—",
        "Detail": "; ".join(detail) if detail else "annotation looks complete",
    }


def run_triage(input_dir: Path, profile, out_csv: Path) -> list[dict]:
    """Triage every GenBank file in input_dir; write a report; return the rows."""
    exts = (".gb", ".gbk", ".gbff", ".genbank", ".gbf")
    files = sorted(f for f in input_dir.iterdir()
                   if f.is_file() and f.suffix.lower() in exts)
    rows: list[dict] = []
    for f in files:
        rec = next(SeqIO.parse(str(f), "genbank"), None)  # first record = the phage
        if rec is None:
            rows.append({"Accession": f.stem, "Organism": "?", "Genome_Size_bp": 0,
                         "CDS_Count": 0, "Has_Structural": "No", "Has_Lysis": "No",
                         "Needs_Pharokka": "Yes", "Trigger": "T4_no_annotation",
                         "Detail": "file unparaseable / empty"})
            continue
        rec._source_path = str(f)
        rows.append(triage_record(rec, profile))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRIAGE_FIELDS)
        w.writeheader()
        w.writerows(rows)

    need = [r for r in rows if r["Needs_Pharokka"] == "Yes"]
    log.info("Pharokka triage: %d/%d genome(s) need re-annotation -> '%s'",
             len(need), len(rows), out_csv)
    for r in need:
        log.info("  NEEDS PHAROKKA  %s (%s): %s",
                 r["Accession"], r["Trigger"], r["Detail"])
    return rows


# =====================================================================
# S5 RBP-AXIS TRIAGE  (v3.3)
# =====================================================================
"""
Second triage axis: flag genomes whose RECEPTOR-BINDING / DEPOLYMERASE annotation
is a GAP, so they are Pharokka-re-annotated BEFORE any GPU is spent on structure
prediction. Runs at S5 PASS-2, after InterPro reconciliation, because the precise
triggers need the InterPro verdict as evidence.

TRIGGERS (>=1, and no exclusion)
  R1  The genome carries a LARGE RBP (>= large_carrier_aa) but has NO enzyme-class
      product term AND no depolymerase verdict -> if a depolymerase domain is
      present it is invisible to both the keyword scan and the whole-protein
      InterPro verdict. This is the architecture in which the S. aureus B6
      depolymerase domains were hiding.
  R2  EVERY RBP candidate in the genome is 'no_domain' -> neither the product
      annotation nor InterPro provides any evidence at all.

EXCLUSION
  X1  Already Pharokka-annotated.

Threshold rationale (literature-grounded, v3.3):
  Characterised standalone capsule depolymerases cluster at ~576-630 aa
  (KP32gp38 576 aa; GBH038_054 577 aa; KP34gp57 630 aa), and curated
  depolymerase sets span ~150-1267 aa (DePP, Magill & Skvortsov 2023), while
  phage RBPs overall run 200-2000+ aa. A protein LARGER than a complete
  standalone depolymerase is therefore a multi-domain CARRIER, which is exactly
  where a depolymerase domain can hide without dominating the annotation. The
  default is set just above the characterised standalone range.
"""

DEFAULT_LARGE_CARRIER_AA = 700

ENZYME_CLASS_TERMS = (
    "lyase", "hydrolase", "glycosid", "depolymer", "esterase", "sialid",
    "dextran", "levan", "rhamnos", "pectate", "pectin", "alginate", "xylan",
    "mannuron", "endorhamnos", "hyaluron", "deacetyl", "amidase", "peptidase",
)
DEPOL_VERDICTS = ("matrix_depolymerase", "ambiguous_depolymerase")

RBP_TRIAGE_FIELDS = ["Accession", "Organism", "N_RBP_Candidates",
                     "Largest_Carrier_aa", "Has_Enzyme_Term", "Has_Depol_Verdict",
                     "Needs_Pharokka", "Trigger", "Detail"]


def run_rbp_triage(candidates_by_genome: dict, out_csv, 
                   large_carrier_aa: int = DEFAULT_LARGE_CARRIER_AA) -> list[dict]:
    """Triage genomes on the RBP/depolymerase axis. `candidates_by_genome` maps an
    accession -> list of RBPCandidate. Returns the rows and writes the report."""
    rows: list[dict] = []
    for acc, cands in sorted(candidates_by_genome.items()):
        if not cands:
            rows.append({
                "Accession": acc, "Organism": "?", "N_RBP_Candidates": 0,
                "Largest_Carrier_aa": 0, "Has_Enzyme_Term": "No",
                "Has_Depol_Verdict": "No", "Needs_Pharokka": "Yes",
                "Trigger": "R2_no_domain_evidence",
                "Detail": "no RBP candidate at all",
            })
            continue
        org = getattr(cands[0], "organism", "?")
        largest = max(c.length for c in cands)
        has_enz = any(any(t in (c.product or "").lower() for t in ENZYME_CLASS_TERMS)
                      for c in cands)
        has_dep = any(c.interpro_verdict in DEPOL_VERDICTS for c in cands)
        all_nodom = all(c.interpro_verdict in ("", "no_domain") for c in cands)

        triggers, detail = [], []
        if largest >= large_carrier_aa and not has_enz and not has_dep:
            triggers.append("R1_large_carrier_no_evidence")
            detail.append(f"{largest} aa carrier but no enzyme term and no "
                          "depolymerase verdict (domain may be hidden)")
        if all_nodom:
            triggers.append("R2_no_domain_evidence")
            detail.append("every RBP candidate is no_domain in InterPro")

        rows.append({
            "Accession": acc, "Organism": org, "N_RBP_Candidates": len(cands),
            "Largest_Carrier_aa": largest,
            "Has_Enzyme_Term": "Yes" if has_enz else "No",
            "Has_Depol_Verdict": "Yes" if has_dep else "No",
            "Needs_Pharokka": "Yes" if triggers else "No",
            "Trigger": ";".join(triggers) if triggers else "—",
            "Detail": "; ".join(detail) if detail else "RBP annotation looks usable",
        })

    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RBP_TRIAGE_FIELDS)
        w.writeheader()
        w.writerows(rows)
    need = [r for r in rows if r["Needs_Pharokka"] == "Yes"]
    log.info("RBP triage: %d/%d genome(s) need Pharokka BEFORE folding -> '%s'",
             len(need), len(rows), out_csv)
    return rows
