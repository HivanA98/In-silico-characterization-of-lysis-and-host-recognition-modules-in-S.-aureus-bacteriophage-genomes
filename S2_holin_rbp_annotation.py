#!/usr/bin/env python3
"""
S2_holin_rbp_annotation.py  (v2, phagecore engine)
=================================================
TABLE 2 (partial) — Holin presence and Tail-fibre/RBP presence + multiplicity.

Thin entry point over `phagecore`. Improvements versus v1:
  * RBP multiplicity (RBP_Count) is reported — phages often carry several
    tail-fibre/RBP genes, which matters for host-range / RBP-focused work;
  * bare-token keywords ("rbp") are matched on word boundaries;
  * a qc_status column propagates the S1 QC verdict so a "No/No" that is really
    a parsing artefact (CDS=0 genome) is visible, not mistaken for biology;
  * --profile / --jobs / --recursive / --manifest as in S1.

Output columns:
  Phage Accession Holin_Present Holin_Evidence
  Tail_Fiber_RBP_Present RBP_Evidence RBP_Count
  annotation_qc Annotation_Note

Usage:
  python S2_holin_rbp_annotation.py -i GenBank -o results\\Table2_holin_rbp.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
from Bio import SeqIO

from phagecore import __version__
from phagecore.genbank_io import (discover_genbank_files, census_features,
                                   FileOutcome, write_manifest)
from phagecore.lysis import detect_holin_rbp
from phagecore.profiles import load_profile
from phagecore.genbank_io import resolve_organism

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("S2")

NOTE = ("Detection by GenBank text annotation only. 'No' may reflect incomplete "
        "annotation, not biological absence.")


def detect_one(path_str: str, profile_spec: str) -> dict:
    profile = load_profile(profile_spec)
    path = Path(path_str)
    try:
        record = next(SeqIO.parse(str(path), "genbank"), None)
        if record is None:
            return {"_status": "empty", "_file": path_str}
        census = census_features(record)
        n_cds = census.get("CDS", 0)
        seq_len = len(record.seq)
        hr = detect_holin_rbp(record, profile)
        # propagate the same FAIL_NO_CDS verdict S1 uses, so artefactual No/No is visible
        qc = ("FAIL_NO_CDS" if n_cds == 0 and seq_len >= profile.min_cds_genome_size_bp
              else "PASS")
        return {
            "_status": "ok", "_file": path_str, "_seq_len": seq_len, "_n_cds": n_cds,
            "Phage": resolve_organism(record, locals().get("path") or locals().get("gb")),
            "Accession": record.id,
            "Holin_Present": "Yes" if hr.holin_present else "No",
            "Holin_Evidence": hr.holin_evidence,
            "Tail_Fiber_RBP_Present": "Yes" if hr.rbp_present else "No",
            "RBP_Evidence": hr.rbp_evidence,
            "RBP_Count": hr.rbp_count,
            "annotation_qc": qc,
            "Annotation_Note": NOTE,
        }
    except Exception as exc:                              # noqa: BLE001
        return {"_status": "parse_error", "_file": path_str, "_msg": str(exc)}


def run(input_dir: Path, output_path: Path, profile_spec: str,
        jobs: int, recursive: bool, manifest_path: Path | None) -> pd.DataFrame:
    profile = load_profile(profile_spec)
    files = discover_genbank_files(input_dir, recursive=recursive)
    if not files:
        log.error("No GenBank files in '%s'.", input_dir)
        sys.exit(1)
    log.info("Found %d GenBank file(s); profile='%s'", len(files), profile.name)

    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            raw = list(ex.map(detect_one, [str(f) for f in files],
                              [profile_spec] * len(files), chunksize=8))
    else:
        raw = []
        for i, f in enumerate(files, 1):
            if i % 200 == 0:
                log.info("  ...processed %d/%d", i, len(files))
            raw.append(detect_one(str(f), profile_spec))

    outcomes, rows = [], []
    for d in raw:
        if d.get("_status") == "ok":
            outcomes.append(FileOutcome(d["_file"], "ok", d["Accession"],
                                        d["Phage"], d["_seq_len"], d["_n_cds"]))
            rows.append({k: v for k, v in d.items() if not k.startswith("_")})
        else:
            outcomes.append(FileOutcome(d.get("_file", "?"),
                                        d.get("_status", "parse_error"),
                                        message=d.get("_msg", "")))
    if not rows:
        log.error("No records extracted.")
        sys.exit(1)

    df = pd.DataFrame(rows).sort_values("Accession").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8")
    log.info("Table 2 (partial) written: '%s' (%d records)", output_path, len(df))

    if manifest_path:
        write_manifest(manifest_path, outcomes, profile.name,
                       f"S2_holin_rbp_annotation.py v{__version__}")
    return df


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="S2_holin_rbp_annotation.py",
        description="Holin & Tail-fibre/RBP detection (phagecore engine).")
    p.add_argument("-i", "--input_dir", type=Path, required=True, metavar="DIR")
    p.add_argument("-o", "--output", type=Path,
                   default=Path("Table2_holin_rbp.csv"), metavar="FILE")
    p.add_argument("--profile", required=True,
                   help="REQUIRED. Built-in name (staphylococcus_aureus | gram_negative_generic | mycobacterium) or path to a .yaml profile. No default (v3.0): prevents silent wrong-host runs.")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--recursive", action="store_true")
    p.add_argument("--manifest", type=Path, default=None)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    df = run(args.input_dir, args.output, args.profile, args.jobs,
             args.recursive, args.manifest)
    sep = "=" * 70
    print(f"\n{sep}\nTABLE 2 (PARTIAL) — HOLIN & TAIL-FIBRE/RBP\n{sep}")
    print(f"  Genomes analysed        : {len(df)}")
    print(f"  Holin detected (Yes)    : {int((df.Holin_Present=='Yes').sum())}")
    print(f"  RBP detected (Yes)      : {int((df.Tail_Fiber_RBP_Present=='Yes').sum())}")
    print(f"  Genomes with >1 RBP     : {int((df.RBP_Count>1).sum())}")
    art = int((df.annotation_qc == 'FAIL_NO_CDS').sum())
    if art:
        print(f"  No/No that are artefacts: {art}   (CDS=0 genomes — not biology)")
    print(sep)
    print("Run S4_endolysin_extractor.py + InterPro for the endolysin columns.")
