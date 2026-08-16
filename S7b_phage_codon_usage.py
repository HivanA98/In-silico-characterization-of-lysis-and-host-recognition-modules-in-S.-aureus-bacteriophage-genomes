#!/usr/bin/env python3
r"""
S7b_phage_codon_usage.py  (v1.0 — companion to S7_codon_trna_coverage.py)
=========================================================================
Phage-vs-host codon usage comparison.

Research question (host-agnostic form):
  "Is the codon usage of these phages similar to that of the bacterial host?"

S7_codon_trna_coverage.py answers a different question: whether phage-encoded
tRNAs decode codons that are RARE in the host. It computes HOST RSCU only. This
file adds the missing half — PHAGE RSCU per genome, and its similarity to the
host — and is a thin CLI over the same engine.

Design constraints (identical to S7 v3.0):
  * The RSCU definition is NOT reimplemented here. count_codons_from_cds() and
    calculate_rscu() are imported verbatim from phagecore.codon, so the family-
    sum invariant cannot drift between S7 and S7b.
  * --profile is REQUIRED (no silent host default), and the same
    host/profile guard is applied.

Outputs (in --output_dir):
  phage_codon_usage.csv        RSCU per codon per phage (long format)
  phage_host_rscu_matrix.csv   codon x phage RSCU matrix + host column
  phage_host_similarity.csv    per-phage Pearson r, Spearman rho, GC3, n_codons

Usage (Windows Command Prompt):
  python S7b_phage_codon_usage.py --phage-dir GenBank\AureusPhage ^
      --host host\Staphylococcus_aureus -o result_aureus_04\codon_analysis ^
      --profile phage_characterization\profiles\AureusPhage.yaml
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phagecore.codon import (  # noqa: E402
    CODON_TO_AA,
    GENETIC_CODE,
    aggregate_host_counts,
    calculate_rscu,
    count_codons_from_cds,
)
from phagecore.profiles import (  # noqa: E402
    available_builtin_profiles,
    load_profile,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

STOPS = set(GENETIC_CODE["Stop"])


def _host_accession(host_path: Path) -> str:
    """Read the accession of the --host genome (best-effort). Same as S7."""
    try:
        return SeqIO.read(str(host_path), "genbank").id
    except Exception:
        return ""


def guard_host_matches_profile(host_path: Path, profile) -> None:
    """Warn (not fail) if --host does not match the profile's host_reference_acc.
    Copied from S7 so both stages behave identically."""
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


def gc3(counts: dict[str, int]) -> float:
    """GC at third codon position, stop codons excluded."""
    gc = tot = 0
    for cod, n in counts.items():
        if cod in STOPS or cod not in CODON_TO_AA:
            continue
        tot += n
        if cod[2] in ("G", "C"):
            gc += n
    return round(100.0 * gc / tot, 2) if tot else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phage-dir", type=Path, required=True,
                    help="Directory of phage GenBank files (one genome per file)")
    ap.add_argument("--host", type=Path, required=True,
                    help="Host genome file or directory (same input as S7 --host)")
    ap.add_argument("--output_dir", "-o", type=Path, required=True)
    ap.add_argument("--glob", default="*.gb*", help="Phage file pattern")
    ap.add_argument("--profile", required=True, metavar="NAME|FILE",
                    help="REQUIRED. Host profile: built-in name "
                         f"({' | '.join(available_builtin_profiles())}) or path to "
                         "a .yaml. No default: prevents silent wrong-host runs.")
    args = ap.parse_args()

    try:
        profile = load_profile(args.profile)
    except (ValueError, FileNotFoundError) as exc:
        log.error(str(exc))
        return 2
    log.info(f"Host profile: '{profile.name}'")
    if args.host.is_file():
        guard_host_matches_profile(args.host, profile)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- host reference RSCU (identical routine to S7) ---
    host_counts, host_files = aggregate_host_counts(args.host)
    host_rscu = calculate_rscu(host_counts)
    print(f"host strains aggregated : {len(host_files)}")
    print(f"host codons counted     : {sum(host_counts.values()):,}")
    print(f"host GC3                : {gc3(host_counts)}%")

    files = sorted(p for p in args.phage_dir.glob(args.glob) if p.is_file())
    if not files:
        print(f"ERROR: no files matching {args.glob} in {args.phage_dir}", file=sys.stderr)
        return 1

    long_rows, sim_rows, matrix = [], [], {}
    for fp in files:
        recs = list(SeqIO.parse(fp, "genbank"))
        if not recs:
            print(f"  SKIP {fp.name}: not parseable as GenBank", file=sys.stderr)
            continue
        if len(recs) > 1:
            print(f"  SKIP {fp.name}: {len(recs)} records — one genome per file required",
                  file=sys.stderr)
            continue
        rec = recs[0]
        acc = rec.id
        counts = count_codons_from_cds(rec)
        n_cod = sum(counts.values())
        if n_cod == 0:
            print(f"  SKIP {acc}: no CDS translations found", file=sys.stderr)
            continue
        rscu = calculate_rscu(counts)
        matrix[acc] = rscu
        for cod, val in rscu.items():
            long_rows.append({"Accession": acc, "Codon": cod,
                              "Amino_Acid": CODON_TO_AA.get(cod, "Stop"),
                              "Count": counts.get(cod, 0), "RSCU": round(val, 4)})

        # correlate only over degenerate families; Met/Trp are invariant at 1.0
        cods = [c for c in rscu
                if c not in STOPS and len(GENETIC_CODE[CODON_TO_AA[c]]) > 1]
        s_ph = pd.Series([rscu[c] for c in cods], index=cods)
        s_ho = pd.Series([host_rscu.get(c, 0.0) for c in cods], index=cods)
        sim_rows.append({
            "Accession": acc,
            "N_codons": n_cod,
            "GC3_phage_pct": gc3(counts),
            "GC3_host_pct": gc3(host_counts),
            "Pearson_r": round(s_ph.corr(s_ho, method="pearson"), 4),
            "Spearman_rho": round(s_ph.corr(s_ho, method="spearman"), 4),
            "Mean_abs_dRSCU": round((s_ph - s_ho).abs().mean(), 4),
        })
        print(f"  {acc:14s} codons={n_cod:7,d}  GC3={gc3(counts):5.2f}%  "
              f"r={sim_rows[-1]['Pearson_r']:.3f}")

    if not sim_rows:
        print("ERROR: no phage genome processed", file=sys.stderr)
        return 1

    pd.DataFrame(long_rows).to_csv(args.output_dir / "phage_codon_usage.csv", index=False)

    mat = pd.DataFrame(matrix).round(4)
    mat.insert(0, "Host_RSCU", pd.Series(host_rscu).round(4))
    mat.insert(0, "Amino_Acid", pd.Series({c: CODON_TO_AA.get(c, "Stop") for c in mat.index}))
    mat.index.name = "Codon"
    mat.to_csv(args.output_dir / "phage_host_rscu_matrix.csv")

    sim = pd.DataFrame(sim_rows).sort_values("Pearson_r", ascending=False)
    sim.to_csv(args.output_dir / "phage_host_similarity.csv", index=False)

    print(f"\ngenomes processed : {len(sim)}")
    print(f"Pearson r  range  : {sim.Pearson_r.min():.3f} – {sim.Pearson_r.max():.3f}"
          f"  (median {sim.Pearson_r.median():.3f})")
    print(f"Spearman   range  : {sim.Spearman_rho.min():.3f} – {sim.Spearman_rho.max():.3f}")
    print(f"written to        : {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
