#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
S5_rbp_extractor.py  (v3.0 — thin CLI over phagecore.rbp)
=========================================================
Extract + classify Receptor-Binding Protein (RBP) / tail-fibre sequences AND the
separate depolymerase enzymatic axis (the MSc thesis variable) from phage GenBank
records. Three-mechanism detection + four QC layers + optional InterPro
reconciliation, then (downstream) ColabFold -> Dali/Foldseek structural validation.

v3.0 change (S5-S7 refactor):
  * Universal DETECTION engine moved into phagecore.rbp (one copy, cannot drift).
  * Host CONTENT (Tier-3) arrives from the PROFILE via phagecore.rbp.configure():
      profile.matrix_depolymerase_terms   thesis-axis enzyme signatures
      profile.ambiguous_matrix_terms       host-specific ambiguity (Staph: GH73)
      profile.s2_annotation_gap_accessions per-host S2="No" rescue list
      profile.depolymerase_known_cases     curated Dali verdicts (source-tagged)
  * --profile REQUIRED (no silent host default; consistent with S1-S4/S6/S7).

Genotype != phenotype: S5 assigns GENOTYPE/capacity. Biofilm-effect DIRECTION and
esterase SUBSTRATE are wet-lab readouts, out of scope.

Windows CMD:
  python S5_rbp_extractor.py -i "GenBank_Bsub" -o results\rbp_candidates.faa ^
      --csv results\rbp_audit.csv --min-length 80 ^
      --profile phage_characterization\profiles\bacillus_subtilis.yaml
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from phagecore import rbp
from phagecore.profiles import load_profile, available_builtin_profiles

log = logging.getLogger("S5")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="S5_rbp_extractor.py",
        description="S5 — RBP/tail-fibre + depolymerase extraction (thin CLI over "
                    "phagecore.rbp; host content from --profile).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input_dir", "-i", type=Path, required=True, metavar="DIR",
                   help="Directory of GenBank files (.gb/.gbk/.gbff).")
    p.add_argument("--output", "-o", type=Path, default=Path("rbp_candidates.faa"),
                   metavar="FILE", help="Combined FASTA (respects --min-length).")
    p.add_argument("--csv", type=Path, default=None, metavar="FILE",
                   help="Audit CSV (all candidates). Fill InterPro_Domain after InterPro.")
    p.add_argument("--min-length", type=int, default=80, metavar="AA",
                   help="Minimum aa length for the FASTA (0 = no filter). Default 80.")
    p.add_argument("--include-length-candidates", action="store_true",
                   help="Enable Mechanism 3 (>=600 aa tail/baseplate context).")
    p.add_argument("--deduplicate", action="store_true",
                   help="Prefer RefSeq (NC_/NZ_) representative per organism.")
    p.add_argument("--interpro", type=Path, default=None, metavar="FILE",
                   help="InterPro TSV to reconcile (fills InterPro_Domain/Verdict).")
    p.add_argument("--deposcope", type=Path, default=None, metavar="FILE|DIR",
                   help="DepoScope output CSV, or a DIRECTORY of them (one per genome). Adds an ORTHOGONAL ESM-2 depolymerase call next to the InterPro verdict: columns DepoScope_Score / DepoScope_Call / Method_Agreement. Neither method overrides the other. DepoScope-only positives are added as high-value ColabFold targets.")
    p.add_argument("--deposcope-threshold", type=float, default=0.5, metavar="P",
                   help="DepoScope probability cut-off (default 0.5, the script's own). FIX THIS BEFORE looking at your candidates — choosing it afterwards is post-hoc selection.")
    p.add_argument("--rbp-triage", type=Path, default=None, metavar="FILE",
                   help="Write an RBP-axis Pharokka triage report: which genomes have a GAP in their receptor-binding/depolymerase annotation and should be Pharokka-re-annotated BEFORE folding. Use at PASS 2 (needs the InterPro verdict).")
    p.add_argument("--rbp-triage-carrier-aa", type=int, default=700, metavar="AA",
                   help="R1 large-carrier threshold (default 700: just above the characterised standalone capsule depolymerase range 576-630 aa, so a larger protein is a multi-domain CARRIER that can hide a depolymerase domain).")
    p.add_argument("--fold-large-rbp", type=int, default=0, metavar="AA",
                   help="RESCUE: also emit receptor-binding carriers >= AA aa as ColabFold targets even without a depolymerase verdict. Use when a host yields no depolymerase call: depolymerase DOMAINS hide inside giant carriers (the S. aureus B6 precedent, 2706-3084 aa). Suggested: 900.")
    p.add_argument("--colabfold-targets", type=Path, default=None, metavar="FILE",
                   help="After --interpro reconciliation, write a clean-header, sequence-deduplicated FASTA of the depolymerase fold targets (<=1200 aa), directly usable in ColabFold/AlphaFold2 (no conversion). Oversized carriers go to a .oversized.txt for domain-splitting.")
    p.add_argument("--interpro-ready", type=Path, default=None, metavar="FILE",
                   help="Write a filtered FASTA (confirmed/high, normal length, "
                        "+ multidomain carriers) for InterPro (100-seq limit).")
    p.add_argument("--profile", required=True, metavar="NAME|FILE",
                   help="REQUIRED. Built-in name (%s) or path to a .yaml profile. "
                        "No default (v3.0)." % " | ".join(available_builtin_profiles()))
    return p


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  [%(levelname)-7s]  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    args = build_parser().parse_args()
    profile = load_profile(args.profile)
    rbp.configure(profile)          # inject host CONTENT before detection
    log.info("Profile: %s | matrix=%s | matrix_terms=%d ambiguous=%d s2gap=%d known_cases=%d",
             profile.name, profile.biofilm_matrix or "(none)",
             len(profile.matrix_depolymerase_terms), len(profile.ambiguous_matrix_terms),
             len(profile.s2_annotation_gap_accessions), len(profile.depolymerase_known_cases))

    per_genome, rescued, confirmed_gap = rbp.run(
        args.input_dir, args.output, args.csv,
        args.min_length, args.include_length_candidates,
        args.deduplicate, args.interpro, args.interpro_ready,
        args.colabfold_targets, args.fold_large_rbp,
        args.rbp_triage, args.rbp_triage_carrier_aa,
        args.deposcope, args.deposcope_threshold,
    )

    n_cand = sum(len(v) for v in per_genome.values())
    log.info("S5 done: %d candidate(s) across %d genome(s); "
             "rescued from S2-No: %d; annotation-gap confirmed: %d",
             n_cand, len(per_genome), len(rescued), len(confirmed_gap))


if __name__ == "__main__":
    main()
