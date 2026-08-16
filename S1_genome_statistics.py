#!/usr/bin/env python3
"""
S1_genome_statistics.py  (v2, phagecore engine)
===============================================
TABLE 1 — General characteristics of phage genomes, with QC and provenance.

This is a THIN entry point over the shared `phagecore` engine. It adds, versus
the v1 script:
  * a QC column set: qc_status, qc_flags, annotation_present, duplicate_group,
    is_representative, taxonomy_flag  — so a 40 kb genome with 0 CDS is FLAGGED
    (FAIL_NO_CDS) instead of silently written as CDS=0;
  * sequence-level deduplication with a documented representative rule;
  * a run manifest (tool versions + per-file status + checksums);
  * optional multiprocessing (--jobs) and recursive discovery (--recursive) for
    massive batches;
  * a --profile switch (REQUIRED in v3.0; no default — prevents wrong-host runs).

Output columns (Table 1 core, unchanged + QC appended):
  Phage Accession Class Family SubFamily Genome_Size_bp GC_Percent
  CDS_Count tRNA_Count NCBI_Status
  | qc_status qc_flags annotation_present duplicate_group is_representative
    taxonomy_flag seq_md5 source_file

Usage (Windows Command Prompt)
------------------------------
  python S1_genome_statistics.py -i GenBank -o results\\Table1.csv
  python S1_genome_statistics.py -i GenBank -o results\\Table1.csv --jobs 8 --recursive
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # find phagecore from any cwd

import pandas as pd
from Bio import SeqIO

from phagecore import __version__
from phagecore.genbank_io import (discover_genbank_files, sequence_md5,
                                   is_refseq_accession, census_features,
                                   FileOutcome, write_manifest)
from phagecore.qc import (QCResult, evaluate_qc, assign_duplicates,
                          merge_dup_into_qc)
from phagecore.taxonomy import resolve_taxonomy, infer_ncbi_status
from phagecore.profiles import load_profile
from phagecore.genbank_io import resolve_organism

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("S1")


def gc_percent(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.00
    return round((s.count("G") + s.count("C")) / len(s) * 100, 2)


# A lightweight, picklable per-genome record (no SeqRecord retained).
def extract_one(path_str: str, profile_spec: str) -> dict:
    profile = load_profile(profile_spec)
    path = Path(path_str)
    try:
        record = next(SeqIO.parse(str(path), "genbank"), None)
        if record is None:
            return {"_status": "empty", "_file": path_str}
        seq = str(record.seq)
        census = census_features(record)
        n_cds = census.get("CDS", 0)
        gc = gc_percent(seq)
        klass, family, subfam, tax_flag = resolve_taxonomy(record, profile)
        return {
            "_status": "ok", "_file": path_str,
            "Phage": resolve_organism(record, locals().get("path") or locals().get("gb")),
            "Accession": record.id,
            "_acc_base": record.id.split(".")[0],
            "Class": klass, "Family": family, "SubFamily": subfam,
            "Genome_Size_bp": len(seq), "GC_Percent": gc,
            "CDS_Count": n_cds,
            "tRNA_Count": census.get("tRNA", 0),
            "NCBI_Status": infer_ncbi_status(record),
            "taxonomy_flag": tax_flag,
            "seq_md5": sequence_md5(seq),
            "_is_refseq": is_refseq_accession(record.id),
            "_census": census, "_seq_len": len(seq), "_gc": gc, "_n_cds": n_cds,
        }
    except Exception as exc:                              # noqa: BLE001
        return {"_status": "parse_error", "_file": path_str, "_msg": str(exc)}


class _ShimGenome:
    """Minimal stand-in so qc.evaluate_qc/assign_duplicates can run post-parse."""
    def __init__(self, d):
        self.accession = d["Accession"]
        self.feature_census = d["_census"]
        self.seq_md5 = d["seq_md5"]
        self.is_refseq = d["_is_refseq"]


def run(input_dir: Path, output_path: Path, profile_spec: str,
        jobs: int, recursive: bool, manifest_path: Path | None) -> pd.DataFrame:
    profile = load_profile(profile_spec)
    files = discover_genbank_files(input_dir, recursive=recursive)
    if not files:
        log.error("No GenBank files in '%s'.", input_dir)
        sys.exit(1)
    log.info("Found %d GenBank file(s); profile='%s'", len(files), profile.name)

    raw: list[dict] = []
    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            raw = list(ex.map(extract_one, [str(f) for f in files],
                              [profile_spec] * len(files), chunksize=8))
    else:
        for i, f in enumerate(files, 1):
            if i % 200 == 0:
                log.info("  ...processed %d/%d", i, len(files))
            raw.append(extract_one(str(f), profile_spec))

    outcomes: list[FileOutcome] = []
    ok = [d for d in raw if d.get("_status") == "ok"]
    for d in raw:
        if d.get("_status") == "ok":
            outcomes.append(FileOutcome(d["_file"], "ok", d["Accession"],
                                        d["Phage"], d["_seq_len"], d["_n_cds"],
                                        d["seq_md5"]))
        else:
            outcomes.append(FileOutcome(d.get("_file", "?"),
                                        d.get("_status", "parse_error"),
                                        message=d.get("_msg", "")))
    if not ok:
        log.error("No records extracted.")
        sys.exit(1)

    # --- QC + dedup (engine) ---
    shims = [_ShimGenome(d) for d in ok]
    dup = assign_duplicates(shims, prefer_refseq=profile.prefer_refseq_representative)
    rows = []
    for d in ok:
        g = _ShimGenome(d)
        qc = evaluate_qc(g, d["_seq_len"], d["_n_cds"], d["_gc"], profile)
        qc = merge_dup_into_qc(qc, dup[d["Accession"]])
        rows.append({
            "Phage": d["Phage"], "Accession": d["Accession"],
            "Class": d["Class"], "Family": d["Family"], "SubFamily": d["SubFamily"],
            "Genome_Size_bp": d["Genome_Size_bp"], "GC_Percent": d["GC_Percent"],
            "CDS_Count": d["CDS_Count"], "tRNA_Count": d["tRNA_Count"],
            "NCBI_Status": d["NCBI_Status"],
            "qc_status": qc.qc_status,
            "qc_flags": ";".join(qc.qc_flags),
            "annotation_present": qc.annotation_present,
            "duplicate_group": qc.duplicate_group,
            "is_representative": qc.is_representative,
            "taxonomy_flag": d["taxonomy_flag"],
            "seq_md5": d["seq_md5"], "source_file": d["_file"],
        })

    df = pd.DataFrame(rows).sort_values("Accession").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8")
    log.info("Table 1 written: '%s' (%d records)", output_path, len(df))

    if manifest_path:
        write_manifest(manifest_path, outcomes, profile.name,
                       f"S1_genome_statistics.py v{__version__}")
    return df


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="S1_genome_statistics.py",
        description="Table 1 genome statistics with QC, dedup and provenance "
                    "(phagecore engine).")
    p.add_argument("-i", "--input_dir", type=Path, required=True, metavar="DIR")
    p.add_argument("-o", "--output", type=Path, default=Path("Table1.csv"),
                   metavar="FILE")
    p.add_argument("--profile", required=True,
                   help="REQUIRED. Built-in name (staphylococcus_aureus | gram_negative_generic | mycobacterium) or path to a .yaml profile. No default (v3.0): prevents silent wrong-host runs.")
    p.add_argument("--jobs", type=int, default=1,
                   help="Parallel workers for massive batches (default 1).")
    p.add_argument("--recursive", action="store_true",
                   help="Recurse into sub-directories of the input.")
    p.add_argument("--manifest", type=Path, default=None,
                   help="Optional run-manifest CSV (versions + per-file status).")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    df = run(args.input_dir, args.output, args.profile, args.jobs,
             args.recursive, args.manifest)

    sep = "=" * 80
    print(f"\n{sep}\nTABLE 1 — GENERAL CHARACTERISTICS (with QC)\n{sep}")
    n_fail = int((df.qc_status == "FAIL_NO_CDS").sum())
    n_warn = int((df.qc_status == "WARN").sum())
    n_dup_nonrep = int((~df.is_representative).sum())
    n_taxflag = int((df.taxonomy_flag != "").sum())
    print(f"  Total genomes           : {len(df)}")
    print(f"  QC PASS                 : {int((df.qc_status=='PASS').sum())}")
    print(f"  QC WARN                 : {n_warn}")
    print(f"  QC FAIL_NO_CDS          : {n_fail}   <-- inspect/re-download these")
    print(f"  Non-representative dups  : {n_dup_nonrep}   "
          f"(unique genomes: {int(df.is_representative.sum())})")
    print(f"  Family needs ICTV check : {n_taxflag}")
    print(f"  Genome size range       : {df.Genome_Size_bp.min():,} – "
          f"{df.Genome_Size_bp.max():,} bp")
    print(sep)
    if n_fail:
        print("  FAIL_NO_CDS accessions:")
        for _, r in df[df.qc_status == "FAIL_NO_CDS"].iterrows():
            print(f"    • {r.Accession}  ({r.Genome_Size_bp:,} bp)  {r.qc_flags}")
        print(sep)
