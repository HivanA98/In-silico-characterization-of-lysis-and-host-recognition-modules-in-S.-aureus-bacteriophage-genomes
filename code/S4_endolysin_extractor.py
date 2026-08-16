#!/usr/bin/env python3
"""
S4_endolysin_extractor.py  (v2, phagecore engine)
================================================
TABLE 2 (lysis enzymes) — endolysin classifier + extractor for InterPro.

Thin entry point over `phagecore`. Improvements versus v1:
  * canonical lowercase classification vocabulary (no more VAPH/vaph split);
  * PROTEIN-LEVEL deduplication for the InterPro FASTA: identical endolysins are
    collapsed to one representative sequence (the 105-set's 92 sequences are only
    59 unique), with a unique-only FASTA produced alongside the per-genome FASTA —
    so InterPro processes the unique set, while prevalence is still per genome;
  * curation registry + length priors come from the profile (Phase-2 ready);
  * audit CSV gains a `shared_by` column listing genomes that carry an identical
    endolysin (the conservation signal made explicit);
  * --profile / --jobs / --recursive / --manifest as in S1/S2.

Outputs
-------
  --output FASTA          : one representative per genome (per-genome view)
  --output-unique FASTA   : deduplicated unique sequences (submit THIS to InterPro)
  --csv audit CSV         : one row per candidate CDS, fully auditable

Usage:
  python S4_endolysin_extractor.py -i GenBank -o results\\endolysin_candidates.faa ^
      --output-unique results\\endolysin_unique.faa --csv results\\endolysin_audit.csv
  # add --run-tblastn (with BLAST+ on PATH) to recover keyword-missed endolysins
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from Bio import SeqIO

from phagecore import __version__
from phagecore.genbank_io import (discover_genbank_files, census_features,
                                   FileOutcome, write_manifest)
from phagecore.lysis import (collect_candidates, select_for_genome, rank_score,
                             find_tblastn, assert_modules_available, Candidate)
from phagecore.interpro import parse_interpro, reconcile
from phagecore.profiles import load_profile
from phagecore.triage import run_triage
from phagecore.genbank_io import resolve_organism

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("S4")

INTERPRO_MAX = 100


def build_reference(input_dir: Path, override: Path | None, profile,
                    recursive: bool) -> Path | None:
    """Resolve the tBLASTn LysK reference FASTA from --reference or the profile."""
    if override is not None and override.is_file():
        log.info("  tBLASTn reference: %s (user-provided)", override)
        return override
    ref_acc = profile.reference_accession
    if not ref_acc:
        return None
    for gb in discover_genbank_files(input_dir, recursive=recursive):
        if ref_acc not in gb.stem:
            continue
        try:
            record = next(SeqIO.parse(str(gb), "genbank"), None)
        except Exception:                                 # noqa: BLE001
            continue
        if record is None:
            continue
        cands = [c for c in collect_candidates(record, profile)
                 if c.classification == "free-endolysin"] or \
                collect_candidates(record, profile)
        if not cands:
            return None
        ref = max(cands, key=lambda c: rank_score(c, profile))
        tmp = Path(tempfile.gettempdir()) / "phagecore_lysK_reference.faa"
        tmp.write_text(f">LysK_reference|{ref.accession}\n{ref.sequence}\n",
                       encoding="utf-8")
        log.info("  tBLASTn reference: auto-extracted %s (%d aa) from %s",
                 ref_acc, ref.length, gb.name)
        return tmp
    log.info("  tBLASTn reference: not found (no %s in input, no --reference).",
             ref_acc)
    return None


def process_one(path_str: str, profile_spec: str, reference_str: str | None,
                run_tblastn_flag: bool, identity_threshold: float) -> dict:
    """Worker: collect candidates + select representative for one genome."""
    profile = load_profile(profile_spec)
    reference = Path(reference_str) if reference_str else None
    path = Path(path_str)
    try:
        record = next(SeqIO.parse(str(path), "genbank"), None)
        if record is None:
            return {"_status": "empty", "_file": path_str}
        census = census_features(record)
        cands = collect_candidates(record, profile)
        pick, msgs = select_for_genome(cands, record, profile, reference,
                                       run_tblastn_flag, identity_threshold)
        # serialise candidates (small: protein sequences only)
        cser = [dict(protein_id=c.protein_id, accession=c.accession,
                     organism=c.organism, product=c.product, length=c.length,
                     inferred_domain=c.inferred_domain,
                     classification=c.classification,
                     selected=c.selected, tblastn_identity=c.tblastn_identity,
                     tblastn_note=c.tblastn_note, runtime_flag=c.runtime_flag,
                     evidence=c.evidence, sequence=c.sequence) for c in cands]
        pser = None
        if pick is not None:
            pser = dict(protein_id=pick.protein_id, accession=pick.accession,
                        organism=pick.organism, product=pick.product,
                        classification=pick.classification,
                        sequence=pick.sequence, fasta_header=pick.fasta_header,
                        length=pick.length)
        return {"_status": "ok", "_file": path_str, "accession": record.id,
                "organism": resolve_organism(record, locals().get("path") or locals().get("gb")),
                "seq_len": len(record.seq), "n_cds": census.get("CDS", 0),
                "candidates": cser, "pick": pser, "msgs": msgs}
    except Exception as exc:                              # noqa: BLE001
        return {"_status": "parse_error", "_file": path_str, "_msg": str(exc)}


def run(args) -> None:
    profile = load_profile(args.profile)
    assert_modules_available(profile)
    files = discover_genbank_files(args.input_dir, recursive=args.recursive)
    if not files:
        log.error("No GenBank files in '%s'.", args.input_dir)
        sys.exit(1)
    log.info("Found %d GenBank file(s); profile='%s'", len(files), profile.name)

    reference = None
    if args.run_tblastn:
        reference = build_reference(args.input_dir, args.reference, profile,
                                    args.recursive)
        if find_tblastn() is None:
            log.warning("  --run-tblastn set but BLAST+ (tblastn) not on PATH; "
                        "commands will be printed instead of executed.")
    ref_str = str(reference) if reference else None

    if args.jobs and args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            raw = list(ex.map(process_one, [str(f) for f in files],
                              [args.profile] * len(files),
                              [ref_str] * len(files),
                              [args.run_tblastn] * len(files),
                              [args.identity_threshold] * len(files),
                              chunksize=8))
    else:
        raw = []
        for i, f in enumerate(files, 1):
            if i % 200 == 0:
                log.info("  ...processed %d/%d", i, len(files))
            raw.append(process_one(str(f), args.profile, ref_str,
                                   args.run_tblastn, args.identity_threshold))

    outcomes, all_candidates, selected, no_pick = [], [], [], []
    for d in raw:
        if d.get("_status") != "ok":
            outcomes.append(FileOutcome(d.get("_file", "?"),
                                        d.get("_status", "parse_error"),
                                        message=d.get("_msg", "")))
            continue
        outcomes.append(FileOutcome(d["_file"], "ok", d["accession"],
                                    d["organism"], d["seq_len"], d["n_cds"]))
        for m in d["msgs"]:
            log.info("  %s", m)
        all_candidates.extend(d["candidates"])
        if d["pick"] is not None:
            selected.append(d["pick"])
        else:
            no_pick.append(d["accession"])

    # --- protein-level dedup (the 92 -> 59 fix) ---
    by_seq: dict[str, list[dict]] = defaultdict(list)
    for p in selected:
        by_seq[p["sequence"]].append(p)
    # representative per unique sequence = lexicographically first accession (stable)
    unique_reps = {}
    shared_by = {}
    for seq, members in by_seq.items():
        members_sorted = sorted(members, key=lambda x: x["accession"])
        rep = members_sorted[0]
        unique_reps[seq] = rep
        shared_by[seq] = [m["accession"] for m in members_sorted]

    # --- optional InterPro reconciliation (authoritative precision filter) ---
    interpro_by_acc: dict = {}
    interpro_summary: dict = {}
    if args.interpro:
        ipath = Path(args.interpro)
        if not ipath.exists():
            log.error("--interpro file not found: %s", ipath)
            sys.exit(1)
        imap = parse_interpro(ipath)
        annotated, interpro_summary = reconcile(list(unique_reps.values()), imap)
        for a in annotated:
            interpro_by_acc[a["accession"]] = (a["interpro_verdict"],
                                               a["interpro_domains"])
        log.info("InterPro reconciliation: %s",
                 ", ".join(f"{k}={v}" for k, v in interpro_summary.items()))

    # --- per-genome FASTA (all selected, one per genome) ---
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as fh:
        for p in sorted(selected, key=lambda x: x["accession"]):
            fh.write(f">{p['fasta_header']}\n{p['sequence']}\n")
    log.info("Per-genome FASTA written: '%s' (%d sequences)",
             args.output, len(selected))

    # --- unique FASTA (submit to InterPro) ---
    if args.output_unique:
        args.output_unique.parent.mkdir(parents=True, exist_ok=True)
        ordered_uni = sorted(unique_reps.items(), key=lambda kv: kv[1]["accession"])
        with open(args.output_unique, "w", encoding="utf-8") as fh:
            for seq, rep in ordered_uni:
                n = len(shared_by[seq])
                fh.write(f">{rep['fasta_header']}|shared_by={n}\n{seq}\n")
        nuni = len(unique_reps)
        log.info("Unique FASTA written:     '%s' (%d unique of %d)",
                 args.output_unique, nuni, len(selected))
        if nuni <= INTERPRO_MAX:
            log.info("  OK %d <= %d: paste the unique file into InterPro at once.",
                     nuni, INTERPRO_MAX)
        else:
            # auto-chunk into <=100-seq files (InterPro limit), like S5 rbp_submit
            n_chunks = (nuni + INTERPRO_MAX - 1) // INTERPRO_MAX
            stem = args.output_unique.with_suffix("")
            suf  = args.output_unique.suffix or ".faa"
            for i in range(n_chunks):
                chunk = ordered_uni[i * INTERPRO_MAX:(i + 1) * INTERPRO_MAX]
                part  = Path(f"{stem}_{i + 1}{suf}")
                with open(part, "w", encoding="utf-8") as fh:
                    for seq, rep in chunk:
                        n = len(shared_by[seq])
                        fh.write(f">{rep['fasta_header']}|shared_by={n}\n{seq}\n")
                log.info("  InterPro chunk %d/%d: '%s' (%d seqs)",
                         i + 1, n_chunks, part.name, len(chunk))
            log.warning("  %d > %d: split into %d files of <=%d. Submit each, "
                        "concatenate the TSVs, then run --interpro on the merged TSV.",
                        nuni, INTERPRO_MAX, n_chunks, INTERPRO_MAX)

        # confirmed-only FASTA (InterPro-verified PG hydrolases)
        if args.interpro:
            confirmed_path = args.output_unique.with_name(
                args.output_unique.stem + "_confirmed" + args.output_unique.suffix)
            n_conf = 0
            with open(confirmed_path, "w", encoding="utf-8") as fh:
                for seq, rep in sorted(unique_reps.items(),
                                       key=lambda kv: kv[1]["accession"]):
                    verdict = interpro_by_acc.get(rep["accession"],
                                                  ("not-in-interpro", ""))[0]
                    if verdict == "confirmed-endolysin":
                        fh.write(f">{rep['fasta_header']}|shared_by="
                                 f"{len(shared_by[seq])}\n{seq}\n")
                        n_conf += 1
            log.info("Confirmed FASTA written:  '%s' (%d InterPro-verified "
                     "endolysins)", confirmed_path, n_conf)

    # --- audit CSV ---
    if args.csv:
        seq_to_shared = {s: v for s, v in shared_by.items()}
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        cols = ["Phage", "Accession", "Protein_ID", "Product", "Length_aa",
                "Inferred_Domain", "InterPro_Domain", "InterPro_Verdict",
                "Classification", "Selected", "tBLASTn_Identity_pct",
                "tBLASTn_Note", "Runtime_Flag", "shared_by", "Evidence"]
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for c in sorted(all_candidates,
                            key=lambda x: (x["accession"], not x["selected"])):
                shared = ""
                if c["selected"] and c["sequence"] in seq_to_shared:
                    others = [a for a in seq_to_shared[c["sequence"]]
                              if a != c["accession"]]
                    shared = ";".join(others)
                ip_dom, ip_verd = "", ""
                if c["selected"]:
                    ip_verd, ip_dom = interpro_by_acc.get(
                        c["accession"], ("", ""))
                w.writerow([
                    c["organism"], c["accession"], c["protein_id"], c["product"],
                    c["length"], c["inferred_domain"], ip_dom, ip_verd,
                    c["classification"], "Yes" if c["selected"] else "No",
                    "" if c["tblastn_identity"] is None
                    else f"{c['tblastn_identity']:.1f}",
                    c["tblastn_note"], c["runtime_flag"], shared, c["evidence"],
                ])
        log.info("Audit CSV written:        '%s' (%d candidate rows)",
                 args.csv, len(all_candidates))

    if getattr(args, "pharokka_triage", None):
        run_triage(args.input_dir, profile, args.pharokka_triage)

    if args.manifest:
        write_manifest(args.manifest, outcomes, profile.name,
                       f"S4_endolysin_extractor.py v{__version__}",
                       extra={"selected": len(selected),
                              "unique_sequences": len(unique_reps),
                              "genomes_without_endolysin": len(no_pick)})

    # --- summary ---
    sep = "=" * 90
    print(f"\n{sep}\nTABLE 2 — ENDOLYSIN EXTRACTION (representative per genome)\n{sep}")
    print(f"  Selected (per genome)   : {len(selected)}")
    print(f"  Unique sequences        : {len(unique_reps)}   <-- submit to InterPro")
    print(f"  Redundant (conserved)   : {len(selected) - len(unique_reps)}")
    print(f"  Genomes without endolysin: {len(no_pick)}")
    if no_pick:
        for a in no_pick:
            print(f"    • {a}  (curation/tBLASTn case or genuinely unannotated)")
    if interpro_summary:
        print(sep)
        print("  InterPro reconciliation (precision filter):")
        for k in ("confirmed-endolysin", "non-endolysin", "no-data",
                  "uncertain", "not-in-interpro"):
            if interpro_summary.get(k):
                print(f"    {k:20s}: {interpro_summary[k]}")
        print(f"  -> Final endolysin set = {interpro_summary.get('confirmed-endolysin', 0)} "
              f"InterPro-confirmed (see *_confirmed.faa)")
    print(sep)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="S4_endolysin_extractor.py",
        description="Endolysin classifier/extractor with protein-level dedup "
                    "for InterPro (phagecore engine).")
    p.add_argument("-i", "--input_dir", type=Path, required=True, metavar="DIR")
    p.add_argument("-o", "--output", type=Path,
                   default=Path("endolysin_candidates.faa"), metavar="FILE")
    p.add_argument("--output-unique", type=Path, default=None, metavar="FILE",
                   help="Deduplicated unique-sequence FASTA (submit to InterPro).")
    p.add_argument("--csv", type=Path, default=None, metavar="FILE")
    p.add_argument("--profile", required=True,
                   help="REQUIRED. Built-in name (staphylococcus_aureus | gram_negative_generic | mycobacterium) or path to a .yaml profile. No default (v3.0): prevents silent wrong-host runs.")
    p.add_argument("--reference", type=Path, default=None,
                   help="LysK reference FASTA for tBLASTn (else from profile).")
    p.add_argument("--run-tblastn", action="store_true")
    p.add_argument("--interpro", type=Path, default=None,
                   help="InterPro result (TSV or pasted summary) to reconcile "
                        "against; fills InterPro_Domain, writes *_confirmed.faa.")
    p.add_argument("--identity-threshold", type=float, default=90.0)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--recursive", action="store_true")
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--pharokka-triage", type=Path, default=None, metavar="FILE",
                   help="Write a Pharokka re-annotation triage report (flags genomes whose annotation is a GAP, not biology: 0-CDS large genomes, no lysis gene, "
                        "accessory-only, or unannotated).")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
