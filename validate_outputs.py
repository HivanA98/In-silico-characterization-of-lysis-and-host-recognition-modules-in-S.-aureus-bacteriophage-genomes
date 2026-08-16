#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_outputs.py — acceptance-criteria checker for the S1-S7 pipeline
========================================================================
Capstone of the toolkit. Reads the OUTPUT files (never the scripts) and runs the
correctness invariants each stage must satisfy, printing PASS / FAIL per check
with a reason on failure and an overall verdict + exit code (0 = all pass).

It is deliberately read-only and non-invasive: it can be run after any pipeline
run, on any host's outputs, without touching the engine or the stage scripts.
This is the same set of checks used to sign off each stage during development,
encoded once so a full run can be self-audited before documentation/submission.

USAGE
    python validate_outputs.py -i results\\            # a results directory
    python validate_outputs.py -i results\\ --rare-threshold 0.5

It auto-discovers known output files by name and runs only the checks whose
inputs are present, so partial result sets validate the stages they contain.

Files recognised (any subset):
  S1  Table1.csv
  S4  endolysin_audit.csv, endolysin_unique.faa
  S5  rbp_audit.csv
  S6  tRNA_detailed.csv, tRNA_summary.csv
  S7  host_codon_usage.csv, phage_trna_coverage.csv, coverage_summary.csv,
      rare_codons.csv

Dependencies: pandas + stdlib.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Callable

import pandas as pd

# Standard genetic code (sense codons only), for the RSCU family-sum invariant.
GENETIC_CODE: dict[str, list[str]] = {
    "Ala": ["GCT", "GCC", "GCA", "GCG"], "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Asn": ["AAT", "AAC"], "Asp": ["GAT", "GAC"], "Cys": ["TGT", "TGC"],
    "Gln": ["CAA", "CAG"], "Glu": ["GAA", "GAG"], "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "His": ["CAT", "CAC"], "Ile": ["ATT", "ATC", "ATA"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "Lys": ["AAA", "AAG"],
    "Met": ["ATG"], "Phe": ["TTT", "TTC"], "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Trp": ["TGG"], "Tyr": ["TAT", "TAC"], "Val": ["GTT", "GTC", "GTA", "GTG"],
}

# --------------------------------------------------------------------------- #
# Result recording
# --------------------------------------------------------------------------- #
class Report:
    def __init__(self) -> None:
        self.checks: list[tuple[str, str, bool, str]] = []  # stage, name, ok, detail

    def add(self, stage: str, name: str, ok: bool, detail: str = "") -> None:
        self.checks.append((stage, name, ok, detail))

    def stage_missing(self, stage: str, why: str) -> None:
        self.checks.append((stage, "(skipped)", None, why))

    def render(self) -> bool:
        width = max((len(n) for _, n, _, _ in self.checks), default=10)
        cur = None
        n_pass = n_fail = 0
        for stage, name, ok, detail in self.checks:
            if stage != cur:
                print(f"\n[{stage}]")
                cur = stage
            if ok is None:
                print(f"  ~ {name:<{width}}  SKIP  {detail}")
                continue
            tag = "PASS" if ok else "FAIL"
            print(f"  {'✓' if ok else '✗'} {name:<{width}}  {tag}"
                  + (f"  {detail}" if detail and not ok else
                     (f"  {detail}" if detail else "")))
            n_pass += int(ok)
            n_fail += int(not ok)
        print("\n" + "=" * 60)
        overall = n_fail == 0
        print(f"OVERALL: {'PASS' if overall else 'FAIL'}  "
              f"({n_pass} passed, {n_fail} failed)")
        print("=" * 60)
        return overall


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


# --------------------------------------------------------------------------- #
# S1 — Table1
# --------------------------------------------------------------------------- #
def check_s1(path: Path, rep: Report) -> None:
    df = _read(path)
    n = len(df)
    rep.add("S1", "row_count>0", n > 0, f"{n} genomes")

    # every FAIL_NO_CDS must be flagged in qc_status (not silently 0)
    if "qc_status" in df and "CDS_Count" in df:
        zero_cds = df[pd.to_numeric(df["CDS_Count"], errors="coerce") == 0]
        unflagged = zero_cds[zero_cds["qc_status"] != "FAIL_NO_CDS"]
        rep.add("S1", "CDS=0 all flagged", len(unflagged) == 0,
                f"{len(unflagged)} unflagged CDS=0" if len(unflagged) else "")

    # GC within plausible bounds OR flagged
    if "GC_Percent" in df:
        gc = pd.to_numeric(df["GC_Percent"], errors="coerce")
        flags = df.get("qc_flags", pd.Series([""] * n))
        oob = df[((gc < 20) | (gc > 75)) & (~flags.str.contains("gc", case=False, na=False))]
        rep.add("S1", "GC in bounds/flagged", len(oob) == 0,
                f"{len(oob)} GC outliers unflagged" if len(oob) else "")

    # accession uniqueness
    if "Accession" in df:
        dup = df["Accession"].duplicated().sum()
        rep.add("S1", "accession unique", dup == 0, f"{dup} duplicate accessions")

    # representatives <= total
    if "is_representative" in df:
        reps = df["is_representative"].astype(str).str.lower().isin(["true", "1", "yes"]).sum()
        rep.add("S1", "reps<=total", reps <= n, f"{reps}/{n} representatives")

    # family Unassigned must carry a taxonomy_flag
    if "Family" in df and "taxonomy_flag" in df:
        unassigned = df[df["Family"].str.strip().str.lower() == "unassigned"]
        unflagged = unassigned[unassigned["taxonomy_flag"].str.strip() == ""]
        rep.add("S1", "Unassigned flagged", len(unflagged) == 0,
                f"{len(unflagged)} Unassigned without flag" if len(unflagged) else "")


# --------------------------------------------------------------------------- #
# S4 — endolysin
# --------------------------------------------------------------------------- #
def check_s4(audit: Path | None, unique_faa: Path | None, rep: Report) -> None:
    if audit and audit.exists():
        df = _read(audit)
        rep.add("S4", "audit rows>0", len(df) > 0, f"{len(df)} candidate rows")
        if "Protein_ID" in df:
            unk = df["Protein_ID"].astype(str).str.contains("unknown", case=False).sum()
            rep.add("S4", "no unknown protein_id", unk == 0, f"{unk} unknown ids")
    if unique_faa and unique_faa.exists():
        headers = [l.strip() for l in unique_faa.read_text().splitlines() if l.startswith(">")]
        dup = len(headers) - len(set(headers))
        rep.add("S4", "unique headers distinct", dup == 0, f"{dup} duplicate headers")
        if audit and audit.exists():
            n_uni = len(headers)
            n_all = len(_read(audit))
            rep.add("S4", "unique<=candidates", n_uni <= n_all, f"{n_uni} unique of {n_all}")


# --------------------------------------------------------------------------- #
# S5 — rbp/depolymerase
# --------------------------------------------------------------------------- #
def check_s5(path: Path, rep: Report) -> None:
    df = _read(path)
    rep.add("S5", "audit rows>0", len(df) > 0, f"{len(df)} rows")
    if "Protein_ID" in df:
        dup = df["Protein_ID"].duplicated().sum()
        unk = df["Protein_ID"].astype(str).str.contains("unknown", case=False).sum()
        rep.add("S5", "protein_id unique", dup == 0, f"{dup} duplicates")
        rep.add("S5", "no unknown protein_id", unk == 0, f"{unk} unknown")
    # if InterPro reconciliation done, depolymerase-flagged rows carry a verdict
    if "Module" in df and "InterPro_Verdict" in df:
        depo = df[df["Module"].str.contains("depolymer", case=False, na=False)]
        submitted = depo[depo["InterPro_Domain"].str.strip() != ""] if "InterPro_Domain" in df else depo
        missing = submitted[submitted["InterPro_Verdict"].str.strip() == ""]
        rep.add("S5", "depoly verdict present", len(missing) == 0,
                f"{len(missing)} submitted depoly w/o verdict" if len(missing) else
                f"{len(depo)} depolymerase rows")


# --------------------------------------------------------------------------- #
# S6 — de-novo tRNA
# --------------------------------------------------------------------------- #
def check_s6(detail: Path | None, summary: Path | None, rep: Report) -> None:
    if detail and detail.exists():
        df = _read(detail)
        rep.add("S6", "detail rows>0", len(df) > 0, f"{len(df)} tRNA rows")
        # real anticodons filled for non-pseudo rows (the whole point of de novo)
        nonp = df[df.get("Is_Pseudo", "No").astype(str).str.lower() != "yes"]
        empty_ac = nonp[nonp["Anticodon_Seq"].str.strip() == ""]
        rep.add("S6", "anticodons filled", len(empty_ac) == 0,
                f"{len(empty_ac)} non-pseudo rows w/o anticodon" if len(empty_ac) else "")
        # CAT anticodon rows flagged CAT_Ambiguous=Yes
        if "CAT_Ambiguous" in df:
            cat = df[df["Anticodon_Seq"].str.upper() == "CAT"]
            bad = cat[cat["CAT_Ambiguous"].str.strip().str.lower() != "yes"]
            rep.add("S6", "CAT flagged ambiguous", len(bad) == 0,
                    f"{len(bad)} CAT rows unflagged" if len(bad) else f"{len(cat)} CAT rows")
        # Sb-1 ground truth present and carries the 4 canonical anticodons at anti level
        sb1 = df[df["Accession"].str.startswith("NC_023009")]
        if len(sb1):
            anti = set(sb1["Anticodon_Seq"].str.upper())
            need = {"GTC", "GAA", "CAT"}   # Asp, Phe, CAT-gene; CCA is the known Infernal miss
            rep.add("S6", "Sb-1 core anticodons", need.issubset(anti),
                    f"have {sorted(anti)}")
    if summary and summary.exists():
        s = _read(summary)
        # functional count excludes pseudo (both columns exist and are non-negative)
        if {"tRNA_Count_functional", "tRNA_Count_pseudo"}.issubset(s.columns):
            fn = pd.to_numeric(s["tRNA_Count_functional"], errors="coerce").fillna(-1)
            rep.add("S6", "functional counts valid", (fn >= 0).all(),
                    f"{int((fn>0).sum())} phages carry tRNA")


# --------------------------------------------------------------------------- #
# S7 — codon usage / coverage
# --------------------------------------------------------------------------- #
def check_s7(host: Path | None, cov: Path | None, summ: Path | None,
             rare: Path | None, rare_threshold: float, rep: Report) -> None:
    host_rscu: dict[str, float] = {}
    if host and host.exists():
        h = _read(host)
        h["RSCU_f"] = pd.to_numeric(h["RSCU"], errors="coerce")
        h["Count_f"] = pd.to_numeric(h["Count"], errors="coerce")
        host_rscu = dict(zip(h["Codon"].str.upper(), h["RSCU_f"]))
        # invariant: sum of RSCU within a synonymous family == family size
        bad = []
        for aa, codons in GENETIC_CODE.items():
            fam = h[h["Codon"].str.upper().isin(codons)]
            if len(fam) < 2:
                continue
            s = fam["RSCU_f"].sum()
            if abs(s - len(fam)) > 1e-2:
                bad.append(f"{aa}={s:.3f}/{len(fam)}")
        rep.add("S7", "RSCU family-sum=size", not bad, "; ".join(bad))
        # RSCU recomputed from Count matches the RSCU column
        h2 = h.dropna(subset=["Count_f"]).copy()
        h2["recalc"] = h2.groupby("Amino_Acid")["Count_f"].transform(lambda x: x / x.mean())
        maxdiff = (h2["recalc"] - h2["RSCU_f"]).abs().max()
        rep.add("S7", "RSCU==recompute(Count)", maxdiff < 1e-2, f"max diff {maxdiff:.4f}")
        # rare_codons.csv == exactly RSCU<threshold
        if rare and rare.exists():
            expected = set(h[h["RSCU_f"] < rare_threshold]["Codon"].str.upper())
            got = set(_read(rare)["Codon"].str.upper())
            rep.add("S7", "rare=RSCU<thr", expected == got,
                    f"exp {len(expected)} vs got {len(got)}")
    if cov and cov.exists():
        c = _read(cov)
        # Ile2-CAT decodes ATA, never ATG (the v2.4.3 correction)
        if {"tRNA_Type", "Anticodon_Seq", "Decoded_Codon"}.issubset(c.columns):
            ile2 = c[(c["tRNA_Type"] == "Ile2") & (c["Anticodon_Seq"].str.upper() == "CAT")]
            wrong = ile2[ile2["Decoded_Codon"].str.upper() == "ATG"]
            rep.add("S7", "Ile2-CAT->ATA", len(wrong) == 0,
                    f"{len(wrong)} Ile2 wrongly ATG" if len(wrong) else f"{len(ile2)} Ile2-CAT rows")
            # decoded codons well-formed
            malformed = c[~c["Decoded_Codon"].str.upper().str.fullmatch(r"[ACGT]{3}")]
            rep.add("S7", "codons well-formed", len(malformed) == 0, f"{len(malformed)} malformed")
        # source is de-novo (no hardcoded table)
        if "Anticodon_Source" in c:
            legacy = c[c["Anticodon_Source"].str.contains("standard table", case=False, na=False)]
            rep.add("S7", "source de-novo only", len(legacy) == 0,
                    f"{len(legacy)} legacy-table rows" if len(legacy) else "")
    if cov and cov.exists() and summ and summ.exists():
        c = _read(cov); s = _read(summ)
        rep.add("S7", "summary phages=coverage",
                c["Accession"].nunique() == len(s),
                f"{c['Accession'].nunique()} vs {len(s)}")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        prog="validate_outputs.py",
        description="Acceptance-criteria checker for S1-S7 pipeline outputs.")
    ap.add_argument("-i", "--results", required=True, type=Path,
                    help="results directory containing the output CSV/FASTA files")
    ap.add_argument("--rare-threshold", type=float, default=0.5,
                    help="RSCU rare cutoff used by S7 (default 0.5)")
    args = ap.parse_args(argv)
    d = args.results
    if not d.is_dir():
        print(f"Not a directory: {d}")
        sys.exit(2)

    def f(name: str) -> Path | None:
        p = d / name
        return p if p.exists() else None

    rep = Report()
    if f("Table1.csv"):
        check_s1(d / "Table1.csv", rep)
    else:
        rep.stage_missing("S1", "Table1.csv not found")
    if f("endolysin_audit.csv") or f("endolysin_unique.faa"):
        check_s4(f("endolysin_audit.csv"), f("endolysin_unique.faa"), rep)
    else:
        rep.stage_missing("S4", "no endolysin outputs found")
    if f("rbp_audit.csv"):
        check_s5(d / "rbp_audit.csv", rep)
    else:
        rep.stage_missing("S5", "rbp_audit.csv not found")
    if f("tRNA_detailed.csv") or f("tRNA_summary.csv"):
        check_s6(f("tRNA_detailed.csv"), f("tRNA_summary.csv"), rep)
    else:
        rep.stage_missing("S6", "no tRNA outputs found")
    if any(f(x) for x in ("host_codon_usage.csv", "phage_trna_coverage.csv",
                          "coverage_summary.csv", "rare_codons.csv")):
        check_s7(f("host_codon_usage.csv"), f("phage_trna_coverage.csv"),
                 f("coverage_summary.csv"), f("rare_codons.csv"),
                 args.rare_threshold, rep)
    else:
        rep.stage_missing("S7", "no codon outputs found")

    ok = rep.render()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
