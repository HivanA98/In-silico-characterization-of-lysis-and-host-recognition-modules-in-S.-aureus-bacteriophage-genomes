#!/usr/bin/env python3
r"""
S4b_endolysin_domain_identity.py  (v1.0 — companion to S4_endolysin_extractor.py)
=================================================================================
Domain-partitioned pairwise identity of endolysins.

Research question (Reviewer 2, Major 3):
  '"Highly conserved" and "markedly variable" carry no measurement. Supply
   pairwise identity across the canonical endolysins resolved by domain — if
   SH3b varies while the catalytic domains do not, recognition is not confined
   to the tail module.'

That question is not answerable from a full-length identity matrix: a modular
enzyme can be near-identical in its catalytic half and divergent in its binding
domain, and the full-length figure averages the two away. This script therefore
partitions each endolysin at its InterPro-resolved domain boundaries, aligns
each domain set SEPARATELY, and reports one pairwise percent-identity matrix per
domain alongside the full-length matrix.

Precedent: staphylococcal endolysin CBDs are known to vary enough that
LysIPLA5 defines a Pfam family (SH3b_T, PF24246) distinct from the SH3_5
(PF08460) family of LysRODI and LysC1C (Vazquez et al., Probiotics Antimicrob
Proteins 2024), and SH3b deletion alters the lytic range of the Kayvirus
endolysin LysF1 (Benesik et al., Virus Genes 2018). The reviewer's hypothesis
is therefore live, and this script is built to be able to confirm OR refute it.

Design constraints
------------------
  * Domain boundaries are NOT re-predicted. They are read from the InterPro TSV
    that S4 already produces, so the domain definition in this script cannot
    drift from the one reported in Table II.
  * Alignment follows the S3 convention: MAFFT L-INS-i, run locally if the
    binary is on PATH, otherwise the per-domain FASTA files are exported for the
    MAFFT web server and read back on a second pass. The manuscript must cite
    whichever route was used.
  * --profile is REQUIRED, as in S1-S7. It is used to stamp provenance and to
    guard against running Staphylococcus domain sets over another host's
    endolysins. NOTE: DOMAIN_SETS below is currently hardcoded for
    staphylococcal endolysins; moving it into the profile YAML is the correct
    place for it once a second host is analysed.

Outputs (in --output_dir)
-------------------------
  domain_regions_<label>.csv        per protein x domain: signature, db, coords
  identity_matrix_<domain>.csv      pairwise % identity matrix (one per domain)
  identity_summary_<label>.csv      n, pairs, min/Q1/median/Q3/max per domain
  align/<domain>.faa                unaligned domain slices (submit these)
  align/<domain>.aln.faa            alignment (written locally, or placed here)

Usage (Windows Command Prompt)
------------------------------
  python S4b_endolysin_domain_identity.py ^
      --faa result_aureus_04\endolysin_unique_AureusPhage_confirmed.faa ^
      --interpro result_aureus_04\interpro_endolysin\endolysin_unique_AureusPhage.tsv ^
      -o result_aureus_04\endolysin_identity ^
      --profile phage_characterization\profiles\AureusPhage.yaml

  If MAFFT is not on PATH the run stops after writing align\*.faa. Align each
  file on https://mafft.cbrc.jp/alignment/server/ with L-INS-i, save the result
  as align\<domain>.aln.faa, and re-run the same command with --align read.
"""

import argparse
import itertools
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
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

# InterPro TSV has no header; these are the first 15 columns of the standard format.
IPR_COLS = ["pid", "md5", "len", "db", "sig", "desc", "start", "end",
            "evalue", "status", "date", "ipr", "iprdesc", "go", "pathway"]

# Signature accessions accepted for each domain class, in preference order.
# Every accession below was observed in the S. aureus endolysin InterPro output;
# nothing here is guessed. Several databases hit the same domain, so the widest
# span among the accepted signatures is taken and the source is recorded.
DOMAIN_SETS: dict[str, list[str]] = {
    "CHAP": ["PF05257", "PS50911"],
    "Amidase": ["PF01510", "SM00644", "cd06583",
                "G3DSA:3.40.80.10", "SSF55846"],
    "SH3b": ["PF08460", "SM00287", "PS51781", "G3DSA:2.30.30.40"],
}


def short_label(pid: str) -> str:
    """Collapse the pipe-delimited FASTA header to 'ACCESSION|PROTEIN'."""
    parts = pid.split("|")
    return f"{parts[0]}|{parts[1]}" if len(parts) > 1 else pid


def load_domains(ipr_path: Path, lengths: dict[str, int]) -> pd.DataFrame:
    """One row per (protein, domain class): the widest accepted hit."""
    df = pd.read_csv(ipr_path, sep="\t", header=None,
                     names=IPR_COLS, usecols=range(15))
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])

    rows = []
    for domain, accessions in DOMAIN_SETS.items():
        hits = df[df["sig"].isin(accessions)].copy()
        if hits.empty:
            log.warning("  no InterPro hit for domain '%s' in any protein", domain)
            continue
        hits["span"] = hits["end"] - hits["start"] + 1
        for pid, grp in hits.groupby("pid"):
            best = grp.sort_values("span", ascending=False).iloc[0]
            rows.append({
                "Protein": short_label(pid),
                "Full_ID": pid,
                "Domain": domain,
                "Signature": best["sig"],
                "Database": best["db"],
                "Start": int(best["start"]),
                "End": int(best["end"]),
                "Domain_Length_aa": int(best["span"]),
                "Protein_Length_aa": lengths.get(pid, -1),
            })
    return pd.DataFrame(rows)


def pct_identity(a: str, b: str) -> tuple[float, int]:
    """
    Percent identity over aligned residues: columns where BOTH sequences carry a
    residue. Columns with a gap in either sequence are excluded from numerator
    and denominator alike. This is the conservative convention — it does not let
    a short, well-matching fragment inflate identity against a long partner, and
    it is the one that must be stated in Methods.
    """
    ident = aligned = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        aligned += 1
        if x.upper() == y.upper():
            ident += 1
    return (100.0 * ident / aligned if aligned else float("nan")), aligned


def run_mafft(infile: Path, outfile: Path) -> bool:
    """MAFFT L-INS-i, the same strategy S3 specifies for TerL."""
    exe = shutil.which("mafft")
    if not exe:
        return False
    cmd = [exe, "--localpair", "--maxiterate", "1000", "--amino", "--quiet",
           str(infile)]
    log.info("  mafft L-INS-i -> %s", outfile.name)
    with open(outfile, "w") as fh:
        res = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        log.error("  MAFFT failed: %s", res.stderr.strip()[:300])
        outfile.unlink(missing_ok=True)
        return False
    return True


def matrix_for(aln_path: Path) -> tuple[pd.DataFrame, list[float]]:
    recs = list(SeqIO.parse(aln_path, "fasta"))
    names = [short_label(r.id) for r in recs]
    seqs = [str(r.seq) for r in recs]
    mat = pd.DataFrame(index=names, columns=names, dtype=float)
    vals = []
    for i in range(len(recs)):
        mat.iloc[i, i] = 100.0
        for j in range(i + 1, len(recs)):
            pid_val, _ = pct_identity(seqs[i], seqs[j])
            mat.iloc[i, j] = mat.iloc[j, i] = round(pid_val, 2)
            vals.append(pid_val)
    return mat, vals


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="S4b_endolysin_domain_identity.py",
        description="Domain-partitioned pairwise identity of endolysins "
                    "(Reviewer 2, Major 3).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--faa", type=Path, required=True, metavar="FILE",
                   help="S4 unique endolysin FASTA "
                        "(endolysin_unique_<Host>_confirmed.faa).")
    p.add_argument("--interpro", type=Path, required=True, metavar="FILE",
                   help="InterPro TSV for the SAME FASTA "
                        "(interpro_endolysin/endolysin_unique_<Host>.tsv).")
    p.add_argument("--output_dir", "-o", type=Path, required=True, metavar="DIR")
    p.add_argument("--label", default=None, metavar="NAME",
                   help="Output filename stem. Default: profile name.")
    p.add_argument("--align", choices=["auto", "local", "export", "read"],
                   default="auto",
                   help="auto (default): local if mafft is on PATH, else export. "
                        "export: write FASTA only. read: use alignments already "
                        "present in <output_dir>/align/<domain>.aln.faa.")
    p.add_argument("--status", default=None, metavar="LIST",
                   help="Comma-separated status values to keep, matched against "
                        "the 'status=' field of the S4 FASTA header (e.g. "
                        "'free-endolysin' to restrict to the canonical set the "
                        "reviewer asks about). Default: keep all.")
    p.add_argument("--min-seqs", type=int, default=3, metavar="N",
                   help="Skip a domain with fewer than N sequences. Default 3.")
    p.add_argument("--profile", required=True, metavar="NAME|FILE",
                   help="REQUIRED. Built-in name "
                        f"({' | '.join(available_builtin_profiles())}) or a "
                        ".yaml path. No default, as in S1-S7.")
    return p


def main() -> int:
    args = build_parser().parse_args()

    try:
        profile = load_profile(args.profile)
    except (ValueError, FileNotFoundError) as exc:
        log.error(str(exc))
        return 2
    label = args.label or profile.name
    log.info("Host profile: '%s'", profile.name)

    if not args.faa.is_file():
        log.error("FASTA not found: %s", args.faa)
        return 2
    if not args.interpro.is_file():
        log.error("InterPro TSV not found: %s", args.interpro)
        return 2

    aln_dir = args.output_dir / "align"
    aln_dir.mkdir(parents=True, exist_ok=True)

    records = {r.id: str(r.seq) for r in SeqIO.parse(args.faa, "fasta")}
    if args.status:
        keep = {s.strip() for s in args.status.split(",") if s.strip()}
        before = len(records)
        records = {k: v for k, v in records.items()
                   if any(f"status={s}" in k for s in keep)}
        log.info("Status filter %s: kept %d of %d sequences",
                 sorted(keep), len(records), before)
    lengths = {k: len(v) for k, v in records.items()}
    log.info("Endolysin sequences: %d", len(records))
    if not records:
        log.error("No sequences in %s", args.faa)
        return 1

    regions = load_domains(args.interpro, lengths)
    if not regions.empty:
        regions = regions[regions.Full_ID.isin(records)].reset_index(drop=True)
    if regions.empty:
        log.error("No accepted domain signatures found in %s. Check that the "
                  "TSV corresponds to this FASTA.", args.interpro)
        return 1
    regions.to_csv(args.output_dir / f"domain_regions_{label}.csv", index=False)

    log.info("Domain coverage (proteins with a resolved boundary):")
    for dom in DOMAIN_SETS:
        sub = regions[regions.Domain == dom]
        log.info("  %-8s %2d/%d", dom, len(sub), len(records))

    # ---- build one FASTA per partition -------------------------------------
    partitions: dict[str, Path] = {}

    full = aln_dir / "full_length.faa"
    with open(full, "w") as fh:
        for pid, seq in records.items():
            fh.write(f">{pid}\n{seq}\n")
    partitions["full_length"] = full

    for dom in DOMAIN_SETS:
        sub = regions[regions.Domain == dom]
        if len(sub) < args.min_seqs:
            log.warning("  domain '%s' has %d sequence(s) — skipped "
                        "(--min-seqs %d)", dom, len(sub), args.min_seqs)
            continue
        path = aln_dir / f"{dom}.faa"
        with open(path, "w") as fh:
            for _, r in sub.iterrows():
                seq = records[r.Full_ID][int(r.Start) - 1:int(r.End)]
                fh.write(f">{r.Full_ID}\n{seq}\n")
        partitions[dom] = path

    # ---- align --------------------------------------------------------------
    mode = args.align
    if mode == "auto":
        mode = "local" if shutil.which("mafft") else "export"
        log.info("Alignment mode resolved to '%s'", mode)

    alignments: dict[str, Path] = {}
    for name, path in partitions.items():
        out = aln_dir / f"{name}.aln.faa"
        if mode == "local":
            if run_mafft(path, out):
                alignments[name] = out
        elif mode == "read":
            if out.is_file():
                alignments[name] = out
            else:
                log.warning("  missing alignment: %s", out)

    if mode == "export" or (mode == "read" and not alignments):
        log.info("")
        log.info("MAFFT is not being run here. Align each file below on")
        log.info("  https://mafft.cbrc.jp/alignment/server/  (L-INS-i),")
        log.info("save each result as <name>.aln.faa in the same folder, then")
        log.info("re-run this command with --align read.")
        for name, path in partitions.items():
            log.info("  %s", path)
        return 0

    if not alignments:
        log.error("No alignments produced.")
        return 1

    # ---- identity matrices --------------------------------------------------
    summary = []
    for name, path in alignments.items():
        mat, vals = matrix_for(path)
        mat.to_csv(args.output_dir / f"identity_matrix_{name}.csv")
        s = pd.Series(vals)
        summary.append({
            "Partition": name,
            "N_sequences": len(mat),
            "N_pairs": len(vals),
            "Min_pct": round(s.min(), 2),
            "Q1_pct": round(s.quantile(0.25), 2),
            "Median_pct": round(s.median(), 2),
            "Q3_pct": round(s.quantile(0.75), 2),
            "Max_pct": round(s.max(), 2),
            "Range_pct": round(s.max() - s.min(), 2),
        })

    order = ["full_length"] + list(DOMAIN_SETS)
    sm = pd.DataFrame(summary)
    sm["_o"] = sm.Partition.apply(lambda x: order.index(x) if x in order else 99)
    sm = sm.sort_values("_o").drop(columns="_o")
    sm.to_csv(args.output_dir / f"identity_summary_{label}.csv", index=False)

    print()
    print(sm.to_string(index=False))
    print()
    print("Percent identity is computed over columns where both sequences carry")
    print("a residue; gap-containing columns are excluded from both numerator")
    print("and denominator. State this convention in Methods.")
    print(f"Written to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
