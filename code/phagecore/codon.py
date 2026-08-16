#!/usr/bin/env python3
"""
phagecore.codon
===============
Universal codon-usage (RSCU) + phage-tRNA coverage engine (S7 engine).

Extracted VERBATIM from S7_codon_trna_coverage.py v2.4.5 — no logic change — so
that the detection engine lives ONCE and cannot drift per host. Everything here
is HOST-AGNOSTIC: the host genome is supplied by the caller (S7 --host), and the
only host input is that genome file. The host-specific anchor (which reference
genome represents the host, e.g. NC_007795.1 for S. aureus, NC_000964.3 for
B. subtilis) is CONFIG carried by the HostProfile (host_reference_acc); it is not
needed by this engine, which reads organism/length from the genome itself.

Universal invariants that MUST NOT change across hosts (all live here now):
  * Ile2-CAT -> ATA (lysidine k2C recoding), Met-CAT -> ATG   [decode_trna]
  * RSCU family-sum invariant (sum within a synonymous family == family size)
  * generic wobble decoding; undetermined-isotype / pseudogene skipping
"""

import logging
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Genetic code (DNA sense strand, 5->3) — universal
# ---------------------------------------------------------------------------
GENETIC_CODE: dict[str, list[str]] = {
    "Phe": ["TTC", "TTT"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Met": ["ATG"],
    "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Tyr": ["TAT", "TAC"],
    "His": ["CAT", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAT", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAT", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["TGT", "TGC"],
    "Trp": ["TGG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "Stop": ["TAA", "TAG", "TGA"],
}

CODON_TO_AA: dict[str, str] = {
    codon: aa
    for aa, codons in GENETIC_CODE.items()
    for codon in codons
}

# ---------------------------------------------------------------------------
# Isotype-aware anticodon -> codon decoding
# ---------------------------------------------------------------------------
# The anticodon comes from S6 (de-novo tRNAscan-SE), so no per-type anticodon
# table is needed. Decoding is generic wobble EXCEPT one biologically essential
# special case, UNIVERSAL across all bacterial hosts:
#
#   tRNA-Ile2 (anticodon CAT): C34 is post-transcriptionally modified to LYSIDINE
#   (k2C) in bacteria, which re-codes the tRNA to read ATA (Ile) and be charged
#   with Ile — NOT ATG (Met). Generic base-pairing would wrongly give CAT -> ATG.
#   tRNAscan-SE labels these "Ile2"; ARAGORN/GenBank labelled the same gene
#   "Met". Using the de-novo isotype call, Ile2-CAT is decoded as ATA. A genuine
#   Met-CAT (tRNAscan-SE type "Met") still decodes ATG.
#
# This rule is why S7 consumes the de-novo TYPE, not just the anticodon: the
# recurring phage CAT tRNA decodes ATA (a codon whose host-rarity is read from
# the host RSCU table), not ATG. This is host-independent biology.

LYSIDINE_ILE2_CODON: str = "ATA"      # Ile2-CAT reads ATA (Ile), not ATG


def decode_trna(trna_type: str, anticodon_seq: str) -> tuple[list[str], str]:
    """
    Return (decoded_codons, note) for a de-novo tRNA using its isotype + anticodon.
    Skips undetermined isotypes (Undet / NNN).
    """
    t = (trna_type or "").strip()
    ac = (anticodon_seq or "").strip().upper().replace("U", "T")
    if not ac or ac == "NNN" or t.lower().startswith("undet"):
        return [], "undetermined isotype/anticodon — not decoded"
    # Ile2 (or an Ile isotype carrying CAT) -> lysidine -> ATA
    if ac == "CAT" and (t == "Ile2" or t.startswith("Ile")):
        return [LYSIDINE_ILE2_CODON], ("Ile2-CAT: lysidine(k2C) modification re-codes "
                                       "C34 to read ATA (Ile), not ATG")
    # a genuine initiator/elongator Met with CAT reads ATG
    if ac == "CAT" and t == "Met":
        return ["ATG"], "Met-CAT decodes ATG (single-codon AA)"
    codons = anticodon_to_codons(ac.lower())
    return codons, "generic wobble decoding from de-novo anticodon"


# ===========================================================================
# Step 1 - Host codon usage and RSCU
# ===========================================================================

def count_codons_from_cds(record: SeqRecord) -> dict[str, int]:
    """Count sense codons from all CDS features. Excludes stop codons and ambiguous bases."""
    counts: dict[str, int] = defaultdict(int)
    stop_set = set(GENETIC_CODE["Stop"])

    for feature in record.features:
        if feature.type != "CDS":
            continue
        try:
            cds_seq = str(feature.extract(record.seq)).upper()
        except Exception:
            continue
        for i in range(0, len(cds_seq) - 2, 3):
            codon = cds_seq[i:i + 3]
            if len(codon) < 3 or "N" in codon or codon in stop_set:
                continue
            if codon in CODON_TO_AA:
                counts[codon] += 1

    return dict(counts)


def calculate_rscu(codon_counts: dict[str, int]) -> dict[str, float]:
    """
    RSCU(codon_i) = observed_count / mean_count_in_synonymous_family.
    Returns 0.0 if the family total is zero.
    """
    rscu: dict[str, float] = {}
    for aa, codons in GENETIC_CODE.items():
        if aa == "Stop":
            continue
        total = sum(codon_counts.get(c, 0) for c in codons)
        mean  = total / len(codons) if total > 0 else 0.0
        for c in codons:
            rscu[c] = round(codon_counts.get(c, 0) / mean, 4) if mean > 0 else 0.0
    return rscu


def _table_from_counts(counts: dict[str, int], rare_threshold: float) -> pd.DataFrame:
    """Build the 61-codon RSCU table from a codon-count dict (single genome or
    aggregated across strains). RSCU maths is identical either way."""
    total = sum(counts.values())
    rscu  = calculate_rscu(counts)
    rows = []
    for aa, codons in GENETIC_CODE.items():
        if aa == "Stop":
            continue
        for c in codons:
            cnt = counts.get(c, 0)
            rows.append({
                "Codon":                  c,
                "Amino_Acid":             aa,
                "Synonymous_Family_Size": len(codons),
                "Count":                  cnt,
                "Freq_per_1000":          round(cnt / total * 1000, 2) if total > 0 else 0.0,
                "RSCU":                   rscu.get(c, 0.0),
                "Rare":                   "Yes" if rscu.get(c, 0.0) < rare_threshold else "No",
            })
    return pd.DataFrame(rows).sort_values(["Amino_Acid", "Codon"]).reset_index(drop=True)


def build_host_codon_table(record: SeqRecord, rare_threshold: float) -> pd.DataFrame:
    """Build the 61-codon RSCU table for one host genome record."""
    return _table_from_counts(count_codons_from_cds(record), rare_threshold)


def _longest_record(path: Path) -> "SeqRecord":
    """Longest record in a GenBank file = the main chromosome (skips plasmids)."""
    recs = list(SeqIO.parse(str(path), "genbank"))
    if not recs:
        raise ValueError(f"empty GenBank file: {path}")
    return max(recs, key=lambda r: len(r.seq))


def aggregate_host_counts(host_dir: Path) -> "tuple[dict[str, int], list[str]]":
    """Sum codon counts across all host strains in a directory (one main
    chromosome per file). Codon usage is a SPECIES-level property; aggregating a
    few RefSeq-complete strains gives a slightly more robust table than any single
    strain, and the intra-species variance is small by design."""
    exts = (".gb", ".gbk", ".gbff", ".genbank", ".gbf")
    files = sorted(f for f in host_dir.iterdir()
                   if f.is_file() and f.suffix.lower() in exts)
    if not files:
        raise ValueError(f"no GenBank files in host dir: {host_dir}")
    total: dict[str, int] = defaultdict(int)
    used = []
    for f in files:
        rec = _longest_record(f)
        for cod, n in count_codons_from_cds(rec).items():
            total[cod] += n
        used.append(rec.id)
    return dict(total), used


# ===========================================================================
# Step 2 - Anticodon to codon mapping
# ===========================================================================

def anticodon_to_codons(anticodon_5to3: str) -> list[str]:
    """
    Map a 5->3 anticodon to the DNA codons it decodes using wobble rules.

    Anticodon pos34 (5'-most) pairs with codon pos3 (wobble):
      G -> C and T (NNY)
      C -> G only  (NNG)
      T/U -> A and G (NNA, NNG)
      A -> T only  (NNT)
    """
    if len(anticodon_5to3) != 3:
        return []
    ac = anticodon_5to3.upper()
    wc = {"A": "T", "T": "A", "C": "G", "G": "C"}
    c1 = wc.get(ac[2])
    c2 = wc.get(ac[1])
    if not c1 or not c2:
        return []
    wobble = ac[0]
    c3_opts = (["C", "T"] if wobble == "G" else
               ["G"]      if wobble == "C" else
               ["A", "G"] if wobble in ("T", "U") else
               ["T"]      if wobble == "A" else [])
    return [f"{c1}{c2}{c3}" for c3 in c3_opts]


def read_s6_trna(trna_csv: Path) -> "dict[tuple[str, str], list[dict]]":
    """
    Read S6's tRNA_detailed.csv (de-novo tRNAscan-SE output). Returns
    {(Phage, Accession): [ {tRNA_Type, Anticodon, Decodes, Note}, ... ]}.
    Pseudogenes and undetermined isotypes are skipped for decoding.
    """
    df = pd.read_csv(trna_csv, dtype=str).fillna("")
    grouped: dict[tuple[str, str], list[dict]] = {}
    for _, r in df.iterrows():
        if str(r.get("Is_Pseudo", "No")).strip().lower() == "yes":
            continue
        ttype = r.get("tRNA_Type", "")
        ac = r.get("Anticodon_Seq", "")
        decodes, note = decode_trna(ttype, ac)
        if not decodes:
            continue
        key = (r.get("Phage", ""), r.get("Accession", ""))
        grouped.setdefault(key, []).append({
            "tRNA_Type": ttype, "Anticodon": ac.upper(),
            "Decodes": decodes, "Note": note})
    return grouped


# ===========================================================================
# Step 3 - Coverage analysis
# ===========================================================================

def analyze_coverage(phage: str, accession: str, trna_entries: list[dict],
                     host_rscu: dict[str, float],
                     rare_threshold: float) -> list[dict]:
    """One row per (phage tRNA, decoded codon) with host RSCU and rare flag."""
    rows = []
    for t in trna_entries:
        for codon in t["Decodes"]:
            codon = codon.upper()
            aa = CODON_TO_AA.get(codon, "?")
            fam = len(GENETIC_CODE.get(aa, [])) if aa != "?" else "?"
            rscu_val = host_rscu.get(codon)
            rows.append({
                "Phage":                   phage,
                "Accession":               accession,
                "tRNA_Type":               t["tRNA_Type"],
                "Anticodon_Seq":           t["Anticodon"],
                "Anticodon_Source":        "de-novo tRNAscan-SE (S6)",
                "Decoded_Codon":           codon,
                "Amino_Acid":              aa,
                "Synonymous_Family_Size":  fam,
                "Host_RSCU":               rscu_val if rscu_val is not None else "N/A",
                "Is_Rare_in_Host":         "Yes" if rscu_val is not None and rscu_val < rare_threshold else "No",
                "Coverage_Note":           t["Note"],
            })
    return rows


def build_summary(df_cov: pd.DataFrame, rare_threshold: float) -> pd.DataFrame:
    """
    One-row-per-phage coverage summary.

    Interpretation distinguishes the two hypotheses:
      (A) translational supplement — a rare codon is covered WITHOUT its common
          synonym from the same amino-acid family being covered.
      (B) anti-defense / non-specific — a rare codon is covered TOGETHER with its
          common synonym (the whole synonymous family is covered).
      (none) no rare codon covered.

    Host-agnostic: which codons are rare comes from the host RSCU table, so the
    same logic applies to any host. The Specific_Rare_Codons column lists only
    rare codons covered WITHOUT their common synonym.
    """
    rows = []
    for (phage, acc), grp in df_cov.groupby(["Phage", "Accession"]):
        covered = set(grp["Decoded_Codon"].astype(str).str.upper())
        rare = grp[grp["Is_Rare_in_Host"] == "Yes"]["Decoded_Codon"].astype(str).str.upper().unique().tolist()
        comm = grp[grp["Is_Rare_in_Host"] == "No"]["Decoded_Codon"].astype(str).str.upper().unique().tolist()

        # A rare codon is "specifically supplemented" only if NO common synonym
        # from the same amino-acid family is also covered by the phage tRNAs.
        specific_rare = []
        for codon in rare:
            aa = CODON_TO_AA.get(codon)
            if aa is None:
                continue
            common_synonyms = set(GENETIC_CODE[aa]) - {codon}
            if not (common_synonyms & covered):
                specific_rare.append(codon)

        if specific_rare:
            interpretation = (
                "Covers rare codon(s) WITHOUT their common synonym ("
                + "; ".join(specific_rare)
                + ") — consistent with translational supplement"
            )
        elif rare:
            interpretation = (
                "Rare codon(s) covered but the full synonymous family is also "
                "covered — consistent with anti-defense or non-specific (NOT "
                "selective rare-codon supplementation)"
            )
        else:
            interpretation = "No rare codon covered"

        rows.append({
            "Phage":                    phage,
            "Accession":                acc,
            "tRNA_Types":               "; ".join(sorted(set(grp["tRNA_Type"]))),
            "Total_Codons_Decoded":     len(covered),
            "Rare_Host_Codons_Covered": len(rare),
            "Rare_Codons":              "; ".join(rare) if rare else "—",
            "Specific_Rare_Codons":     "; ".join(specific_rare) if specific_rare else "—",
            "Common_Codons_Covered":    len(comm),
            "Common_Codons":            "; ".join(comm) if comm else "—",
            "RSCU_Threshold":           rare_threshold,
            "Interpretation":           interpretation,
        })
    return pd.DataFrame(rows).sort_values("Accession").reset_index(drop=True)


# ===========================================================================
# Batch driver
# ===========================================================================

def run(host_path: Path, trna_csv: Path, output_dir: Path,
        rare_threshold: float) -> None:
# --- Host ---
    log.info(f"Reading host genome: '{host_path}'")
    try:
        if host_path.is_dir():
            # --host-dir: aggregate codon counts across all strains in the folder
            counts, used = aggregate_host_counts(host_path)
            host_name = host_path.name
            host_id   = f"{len(used)} strains: " + ", ".join(used)
            log.info(f"  [aggregate] summed codon counts across {len(used)} strain(s): {used}")
            df_host = _table_from_counts(counts, rare_threshold)
        else:
            # single genome; longest record = main chromosome (skips plasmids)
            host_rec  = _longest_record(host_path)
            host_name = host_rec.annotations.get("organism", host_rec.id)
            host_id   = host_rec.id
            recs_n    = len(list(SeqIO.parse(str(host_path), "genbank")))
            if recs_n > 1:
                log.info(f"  [Info] {recs_n} records; using longest: {host_rec.id}")
            df_host = build_host_codon_table(host_rec, rare_threshold)
    except Exception as exc:
        log.error(f"Cannot read host genome: {exc}")
        sys.exit(1)

    log.info(f"  {host_name} ({host_id})")

    host_rscu = dict(zip(df_host["Codon"], df_host["RSCU"]))
    n_rare    = (df_host["Rare"] == "Yes").sum()
    log.info(f"  Total sense codons: {df_host['Count'].sum():,} | "
             f"Rare (RSCU<{rare_threshold}): {n_rare}/61")

    # --- V&V: RSCU arithmetic invariant (sum of RSCU within a synonymous family
    #     must equal the family size). A violation means a computation bug. ---
    vv_ok = True
    for aa, codons in GENETIC_CODE.items():
        fam = [c for c in codons if c in host_rscu]
        if len(fam) < 2:
            continue
        s = sum(host_rscu[c] for c in fam)
        if abs(s - len(fam)) > 1e-3:
            vv_ok = False
            log.warning(f"  [V&V] RSCU sum for {aa} = {s:.4f} != family size {len(fam)}")
    log.info(f"  [V&V] RSCU family-sum invariant: {'PASS' if vv_ok else 'FAIL'}")

    # --- Phage tRNAs: consume S6 de-novo output (NOT GenBank features) ---
    trna_by_phage = read_s6_trna(trna_csv)
    if not trna_by_phage:
        log.error(f"No decodable tRNA rows in '{trna_csv}'.")
        sys.exit(1)
    log.info(f"\nS6 de-novo tRNA: {len(trna_by_phage)} phages with >=1 decodable tRNA")

    all_cov: list[dict] = []
    for (phage, acc), entries in sorted(trna_by_phage.items(), key=lambda kv: kv[0][1]):
        cov = analyze_coverage(phage, acc, entries, host_rscu, rare_threshold)
        all_cov.extend(cov)
        n_rare_cov = sum(1 for r in cov if r["Is_Rare_in_Host"] == "Yes")
        log.info(f"  {acc:<14}  {len(entries)} tRNA(s)  "
                 f"-> {len(set(r['Decoded_Codon'] for r in cov))} codons  "
                 f"| rare covered: {n_rare_cov}")

    # --- Write outputs ---
    output_dir.mkdir(parents=True, exist_ok=True)

    df_host.to_csv(output_dir / "host_codon_usage.csv", index=False, encoding="utf-8")
    df_rare = df_host[df_host["Rare"] == "Yes"].sort_values("RSCU").reset_index(drop=True)
    df_rare.to_csv(output_dir / "rare_codons.csv", index=False, encoding="utf-8")
    log.info(f"\nhost_codon_usage.csv  written ({len(df_host)} codons)")
    log.info(f"rare_codons.csv       written ({len(df_rare)} rare codons)")

    if all_cov:
        df_cov = pd.DataFrame(all_cov)
        df_sum = build_summary(df_cov, rare_threshold)
        df_cov.to_csv(output_dir / "phage_trna_coverage.csv", index=False, encoding="utf-8")
        df_sum.to_csv(output_dir / "coverage_summary.csv", index=False, encoding="utf-8")
        log.info(f"phage_trna_coverage.csv written ({len(df_cov)} rows)")
        log.info(f"coverage_summary.csv  written ({len(df_sum)} phages)")
    else:
        df_cov = df_sum = pd.DataFrame()
        log.info("No tRNA-positive phages found.")

    # --- Terminal print ---
    sep = "=" * 75
    print(f"\n{sep}")
    print("S7 CODON / tRNA COVERAGE — RESULTS")
    print(sep)
    print(f"  Host: {host_name} ({host_id})")
    print(f"  Rare codons (RSCU < {rare_threshold}): {n_rare} of 61")
    if not df_rare.empty:
        print(f"\n  Rare host codons:")
        for _, r in df_rare.iterrows():
            print(f"    {r['Codon']}  ({r['Amino_Acid']:3s})  RSCU={r['RSCU']:.3f}  "
                  f"Freq={r['Freq_per_1000']:.2f}/1000")
    if not df_sum.empty:
        print(f"\n  Coverage per phage:")
        for _, r in df_sum.iterrows():
            print(f"  {r['Accession']:<18} RareCovered={r['Rare_Host_Codons_Covered']}  "
                  f"{r['Interpretation']}")
    print(sep)
    print(f"  Output: {output_dir.resolve()}")
    print(sep)

    return {
        "host_name": host_name, "host_id": host_id,
        "n_rare": int(n_rare), "vv_ok": vv_ok,
    }
