#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
S6_tRNA_analyzer.py  (v3.0 — thin CLI over phagecore.trna)
==========================================================
De-novo tRNA analysis orchestrator (tRNAscan-SE). Two modes:

Mode 1 — PREPARE (default): write ONE nucleotide multi-FASTA of all genomes to
  submit to the tRNAscan-SE web server (https://trna.ucsc.edu/tRNAscan-SE/,
  Sequence source: Bacterial). Save the TABULAR output for Mode 2.

Mode 2 — PARSE + ANALYSE: parse the tRNAscan-SE table, join per phage, write a
  per-tRNA detail CSV (with REAL anticodons) and a per-phage summary. GenBank is
  read only for metadata + the ground-truth check.

v3.0 change (S5-S7 refactor):
  * Detection/parsing engine moved VERBATIM into phagecore.trna (universal: CAT
    ambiguity, pseudo/Undet handling live there ONCE and cannot drift per host).
  * Host anchors + signature come from the PROFILE, not hardcoded constants:
      profile.trna_ground_truth_acc      primary calibration anchor
      profile.trna_secondary_anchor_acc  secondary anchor (concordance)
      profile.trna_canonical_isotypes    DERIVED signature (pattern + Has_<iso>)
  * --profile is REQUIRED (no silent host default; consistent with S1-S4/S7 so
    S8/S9 inherit no backdoor).

Ground-truth gate is anticodon-level when GenBank supplies anticodons (Staph
Sb-1) and ISOTYPE-level when it does not (B. subtilis phages lack GenBank
anticodons) — reported honestly with a PASS/REVIEW verdict; ARAGORN anticodon
concordance is the pending strengthening for the isotype-level hosts.

Usage (Windows Command Prompt)
------------------------------
  REM Mode 1 (prepare FASTA):
  python S6_tRNA_analyzer.py -i "GenBank_Bsub" -o results\trnascan_input.fasta ^
      --profile phage_characterization\profiles\bacillus_subtilis.yaml

  REM Mode 2 (parse tRNAscan-SE results):
  python S6_tRNA_analyzer.py -i "GenBank_Bsub" --results results\trnascan_results.txt ^
      -d results\tRNA_detailed.csv -s results\tRNA_summary.csv ^
      --profile phage_characterization\profiles\bacillus_subtilis.yaml
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from phagecore import trna
from phagecore.profiles import load_profile, available_builtin_profiles

LOG = logging.getLogger("S6")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="S6_tRNA_analyzer.py",
        description="S6 — de-novo tRNA analysis via tRNAscan-SE (thin CLI over "
                    "phagecore.trna; host anchors/signature from --profile).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, type=Path,
                   help="Folder of GenBank genomes.")
    p.add_argument("-o", "--output", type=Path,
                   help="Mode 1: nucleotide FASTA to submit to tRNAscan-SE.")
    p.add_argument("--chunk-mbp", type=float, default=None,
                   help="Mode 1: split the FASTA into <=N-Mbp parts (e.g. 5).")
    p.add_argument("--results", type=Path,
                   help="Mode 2: tRNAscan-SE tabular result file to parse.")
    p.add_argument("-d", "--detail", type=Path,
                   default=Path("results/tRNA_detailed.csv"),
                   help="Mode 2: per-tRNA detail CSV out.")
    p.add_argument("-s", "--summary", type=Path,
                   default=Path("results/tRNA_summary.csv"),
                   help="Mode 2: per-phage summary CSV out.")
    p.add_argument("--profile", required=True, metavar="NAME|FILE",
                   help="REQUIRED. Built-in name (%s) or path to a .yaml profile. "
                        "No default (v3.0): prevents silent wrong-host runs."
                        % " | ".join(available_builtin_profiles()))
    return p


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  [%(levelname)-7s]  %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    args = build_parser().parse_args(argv)
    profile = load_profile(args.profile)
    LOG.info("Profile: %s | tRNA anchor=%s secondary=%s canonical=%s",
             profile.name, profile.trna_ground_truth_acc or "(none)",
             profile.trna_secondary_anchor_acc or "(none)",
             ",".join(profile.trna_canonical_isotypes) or "(none)")
    if args.results:                      # Mode 2
        trna.analyse(args.input, args.results, args.detail, args.summary, profile)
    elif args.output:                     # Mode 1
        trna.write_nucleotide_fasta(args.input, args.output, args.chunk_mbp)
    else:
        LOG.error("Choose a mode: -o FASTA (prepare) OR --results FILE (analyse).")
        sys.exit(2)


if __name__ == "__main__":
    main()
