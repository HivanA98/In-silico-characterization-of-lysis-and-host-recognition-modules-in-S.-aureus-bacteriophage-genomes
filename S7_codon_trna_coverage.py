#!/usr/bin/env python3
r"""
S7_codon_trna_coverage.py  (v3.0 — thin CLI over phagecore.codon)
=================================================================
Codon usage (RSCU) + phage-tRNA coverage analysis.

Research question (host-agnostic form):
  "Do phage-encoded tRNAs decode codons that are RARE in the bacterial host?"
  A phage tRNA whose de-novo isotype/anticodon decodes a host-rare codon (low
  host RSCU) is a candidate translational supplement / anti-defense element.

v3.0 change (S5-S7 refactor):
  * The engine moved VERBATIM into phagecore.codon (universal; the Ile2-CAT->ATA
    lysidine rule and the RSCU family-sum invariant live there ONCE and cannot
    drift per host). This file is now a thin CLI.
  * --profile is REQUIRED (no silent host default; consistent with S1-S4 and so
    S8/S9 inherit no backdoor). The host GENOME is still supplied by --host; the
    profile carries host_reference_acc only to GUARD against a host/profile
    mismatch (e.g. running B. subtilis phages against the S. aureus host genome).

Host is fully externalised:
  * S. aureus  : --host <NC_007795.gb>  --profile staphylococcus_aureus
  * B. subtilis: --host <NC_000964.gb>  --profile phage_characterization\profiles\bacillus_subtilis.yaml
                 (NC_000964.3 = B. subtilis subsp. subtilis str. 168, the standard
                  reference; download from NCBI, save as host\NC_000964.gb)

Outputs (in --output_dir)
-------
  host_codon_usage.csv    - 61 sense codons: count, freq/1000, RSCU, Rare
  rare_codons.csv         - codons with RSCU < threshold
  phage_trna_coverage.csv - per (phage tRNA x decoded codon): host RSCU, rare flag
  coverage_summary.csv    - one row per phage with interpretation

Usage (Windows Command Prompt)
------------------------------
  python S7_codon_trna_coverage.py --host host\NC_000964.gb ^
      --trna results\tRNA_detailed.csv -o results\codon_analysis ^
      --profile phage_characterization\profiles\bacillus_subtilis.yaml
"""

import argparse
import logging
import sys
from pathlib import Path

from Bio import SeqIO

from phagecore import codon
from phagecore.profiles import load_profile, available_builtin_profiles

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_RARE_THRESHOLD: float = 0.5


def _host_accession(host_path: Path) -> str:
    """Read the accession/id of the --host genome (best-effort)."""
    try:
        return SeqIO.read(str(host_path), "genbank").id
    except Exception:
        return ""


def guard_host_matches_profile(host_path: Path, profile) -> None:
    """
    Warn if the --host genome does not match the profile's host_reference_acc.
    A guard against silently running the wrong host genome; it does NOT hard-fail
    (a user may legitimately use an updated assembly) but makes a mismatch loud.
    """
    expected = (profile.host_reference_acc or "").split(".")[0]
    if not expected:
        log.info(f"  profile '{profile.name}' declares no host_reference_acc — "
                 f"host guard skipped.")
        return
    got = _host_accession(host_path).split(".")[0]
    if got and got != expected:
        log.warning(
            f"  HOST/PROFILE MISMATCH: --host is {got} but profile "
            f"'{profile.name}' expects {expected} ({profile.host_reference_acc}). "
            f"Check you passed the right host genome for this host.")
    else:
        log.info(f"  host guard OK: {got or '(unread)'} matches profile "
                 f"'{profile.name}' ({profile.host_reference_acc}).")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="S7_codon_trna_coverage.py",
        description="Codon usage (RSCU) + phage tRNA coverage. Host-agnostic engine "
                    "(phagecore.codon); host supplied by --host, guarded by --profile.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=r"""
Host genome (download as GenBank, then pass with --host):
  S. aureus  : NC_007795.1 (NCTC 8325)      -> host\NC_007795.gb
  B. subtilis: NC_000964.3 (str. 168)       -> host\NC_000964.gb

Example (Windows CMD):
  python S7_codon_trna_coverage.py --host host\NC_000964.gb ^
      --trna results\tRNA_detailed.csv -o results\codon_analysis ^
      --profile phage_characterization\profiles\bacillus_subtilis.yaml
""",
    )
    p.add_argument("--host", type=Path, required=True, metavar="FILE|DIR",
                   help="Reference host genome (GenBank). A DIRECTORY aggregates codon counts across all strains inside it (--host-dir mode) — codon usage is a species-level property, so summing a few RefSeq-complete strains is fine.")
    p.add_argument("--trna", "-i", type=Path, required=True, metavar="FILE",
                   help="S6 tRNA_detailed.csv (de-novo tRNAscan-SE output).")
    p.add_argument("--output_dir", "-o", type=Path,
                   default=Path("codon_analysis"), metavar="DIR",
                   help="Output directory. Default: codon_analysis")
    p.add_argument("--rare-threshold", type=float, default=DEFAULT_RARE_THRESHOLD,
                   metavar="RSCU",
                   help=f"RSCU below which a codon is 'rare'. Default {DEFAULT_RARE_THRESHOLD}.")
    p.add_argument("--profile", required=True, metavar="NAME|FILE",
                   help="REQUIRED. Host profile: built-in name "
                        f"({' | '.join(available_builtin_profiles())}) or path to a "
                        ".yaml. No default (v3.0): prevents silent wrong-host runs.")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()

    try:
        profile = load_profile(args.profile)
    except (ValueError, FileNotFoundError) as exc:
        log.error(str(exc))
        sys.exit(2)

    log.info(f"Host profile: '{profile.name}'")
    if args.host.is_file():
        guard_host_matches_profile(args.host, profile)

    codon.run(args.host, args.trna, args.output_dir, args.rare_threshold)
