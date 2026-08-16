#!/usr/bin/env python3
"""
S3_terl_extractor.py
====================
Terminase Large Subunit (TerL) Extractor for Phylogenetic Analysis

==============================================================================
MANUSCRIPT CONTRIBUTION
  Generates input for:  FIGURE 1 — Maximum Likelihood Phylogenetic Tree

  Complete Figure 1 pipeline:
    Step 1 — This script   → TerL_combined.faa  (multi-FASTA)
    Step 2 — MAFFT         → TerL_aligned.faa   (multiple sequence alignment)
    Step 3 — MEGA 12.1.2   → Figure 1           (Maximum Likelihood tree)
  See "Downstream Phylogenetic Workflow" below for exact commands and MEGA
  parameter settings.
==============================================================================

Associated manuscript:
    "Molecular Characterization of Lytic Bacteriophages Against Resistant
    Staphylococcus aureus Based on NCBI GenBank Sequences:
    A Bioinformatic Literature Review"

Description
-----------
Extracts terminase large subunit (TerL) amino acid sequences from a directory
of GenBank files and writes a single combined multi-FASTA file for direct
alignment and phylogenetic analysis.

TerL is the standard phylogenetic marker for tailed bacteriophage
classification: universally present, under strong purifying selection,
and provides clear inter-family resolution (Meier-Kolthoff & Goker, 2017).

Historical Background and Annotation Problem Resolved
------------------------------------------------------
During the original analysis, the extraction script was initially unable to
identify TerL in 9 of 22 genomes. Investigation revealed two distinct
situations:

  Case A — Kayvirus group (7 genomes: NC_047722, NC_047723, NC_047724,
  NC_047725, NC_047726, NC_047727, EU418428):
    The TerL protein in these genomes (605 aa) is annotated in the GenBank
    "product" qualifier with the 3-letter abbreviation "Ter" rather than the
    standard text "terminase large subunit". This non-standard annotation
    causes a failure in standard keyword searches.
    FIX: EXACT_NAMES = {"Ter", "ter"} — checks if the product qualifier
    EXACTLY EQUALS "Ter" or "ter" (case-sensitive). This catches all 7
    Kayvirus genomes and is the critical correction in this script.

  Case B — Staphylococcus phage Portland (MT926124) and
  Staphylococcus phage vB_SauP-436A1 (MN150710):
    These two genomes genuinely have no annotated TerL sequence. Detailed
    inspection reveals:
      Portland (MT926124):        "putative encapsidation protein" (415 aa)
      vB_SauP-436A1 (MN150710):  "DNA packaging protein" (415 aa)
    Both are micro-phages (~17–18 kb genome) that likely do not belong to
    Myovirus morphology and use a packaging mechanism not annotated under
    standard TerL nomenclature. Exclusion from phylogenetic analysis is
    scientifically justified.
    DECISION: These two genomes are CORRECTLY EXCLUDED from the output FASTA
    and from Figure 1.

Methods Statement for Manuscript
---------------------------------
The following text was incorporated into the Methods section:

  "Staphylococcus phage Portland (MT926124) and vB_SauP-436A1 (MN150710)
  were excluded from phylogenetic analysis due to the absence of annotated
  terminase large subunit sequences, consistent with their atypical small
  genome sizes (<20 kb) relative to the remaining dataset."

Detection Logic (Two Mechanisms)
----------------------------------
For each CDS feature in a GenBank record, two checks are applied:

  Mechanism 1 — Keyword match:
    Any string from TERL_KEYWORDS is found as a substring in the lowercased
    "product" qualifier. Catches standard annotations such as
    "terminase large subunit", "large terminase", "TerL protein", etc.

  Mechanism 2 — Exact product match:
    The product qualifier EXACTLY EQUALS "Ter" or "ter" (case-sensitive,
    checked against EXACT_NAMES). This is the fix for the Kayvirus group.

The FIRST CDS feature satisfying either mechanism AND containing a valid
"translation" qualifier is returned as the TerL for that genome.

Output FASTA Header Format
--------------------------
    >{accession}|{organism_no_spaces}|{product_annotation}
Example:
    >NC_047722.1|Staphylococcus_phage_Staph1N|Ter

This format is directly importable by MAFFT and MEGA 12.1.2.

Downstream Phylogenetic Workflow (Figure 1)
-------------------------------------------
Step 1 — Run this script:
    python S3_terl_extractor.py -i GenBank -o results\\TerL_combined.faa

Step 2 — Multiple Sequence Alignment with MAFFT (WEB SERVER, not local app):
    Open the MAFFT online server:
      https://mafft.cbrc.jp/alignment/server/

    a. Upload results\\TerL_combined.faa, OR paste the FASTA content.
    b. At ~150 sequences, use strategy:
         G-INS-i   (accurate, global homology; good for one conserved domain)
       L-INS-i is also acceptable but slower at this size. Either is fine as
       long as the SAME strategy is used for every rebuild (reproducibility).
    c. Click "Submit".
    d. Save the "Fasta format" result as results\\TerL_aligned.faa.

Step 3 — Maximum-Likelihood tree with IQ-TREE (the tree is BUILT here):
    At 150 taxa, IQ-TREE 2 is used instead of MEGA's built-in ML: it does
    automatic model selection (ModelFinder) and ultrafast bootstrap, which is
    faster and more rigorous than MEGA for a dataset this size.
    Download IQ-TREE 2: http://www.iqtree.org/
    Run (one line):
      iqtree2 -s results\\TerL_aligned.faa -m MFP -B 1000 -alrt 1000 -T AUTO
        -s        : the MAFFT alignment
        -m MFP    : ModelFinder Plus picks the best substitution model
                    (it will typically select LG+G+I or similar for TerL)
        -B 1000   : 1000 ultrafast bootstrap replicates
        -alrt 1000: SH-aLRT branch test (report both supports)
        -T AUTO   : auto-detect CPU threads
    Output: results\\TerL_aligned.faa.treefile  (Newick, with supports)
    Optional outgroup: add  -o EW_accession_label  (the EW tip label).

Step 4 — Display in MEGA 12 as a USER TREE (MEGA no longer computes the tree):
    Download MEGA 12: https://www.megasoftware.net/
    Menu: User Tree > Display Newick Tree
    Select the IQ-TREE **`.contree`** file (bootstrap consensus tree). IQ-TREE
    writes two trees: `.treefile` (the ML tree) and `.contree` (the extended
    majority-rule consensus of the ultrafast-bootstrap replicates, with support
    values); the `.contree` is used for display here. MEGA renders the imported
    topology + supports for figure layout — the tree is not recomputed, so the
    figure is exactly the IQ-TREE result.

Reproducibility Note (tree topology vs support values)
----------------------------------------------------------
The TREE TOPOLOGY is the reportable, stable result, reproducible given the same
input and the same MAFFT strategy. IQ-TREE ultrafast bootstrap is seeded, so
re-runs on the same alignment are far more reproducible than MEGA's resampling;
to fix the seed exactly, add  -seed 12345  to the IQ-TREE command. Report the
topology and approximate UFBoot/SH-aLRT supports, not exact integers, and always
reuse the same MAFFT strategy + IQ-TREE settings for any rebuild.

Reference
---------
Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I,
Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009). Biopython: freely
available Python tools for computational molecular biology and bioinformatics.
Bioinformatics, 25(11):1422-1423. doi:10.1093/bioinformatics/btp163

Meier-Kolthoff JP, Goker M (2017). VICTOR: genome-based phylogeny and
classification of prokaryotic viruses. Bioinformatics, 33(21):3396-3404.
doi:10.1093/bioinformatics/btx440

Tested Environment (Windows)
-----------------------------
    OS         : Windows 10 / 11
    Python     : 3.12.10
    biopython  : 1.87
    MAFFT      : WEB SERVER (https://mafft.cbrc.jp/alignment/server/),
                 strategy L-INS-i — no local MAFFT installation required
    MEGA       : 12.1.2 (https://www.megasoftware.net)

Installation (Command Prompt / PowerShell)
------------------------------------------
    pip install biopython==1.87

Usage (Windows Command Prompt)
-------------------------------
    python S3_terl_extractor.py -i GenBank -o results\\TerL_combined.faa
"""

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  [%(levelname)-8s]  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VALID_EXTENSIONS: frozenset[str] = frozenset({".gb", ".gbk", ".gbff"})

# ---- Mechanism 1: Keyword match on "product" qualifier ----
# Matched case-insensitively as substrings.
TERL_KEYWORDS: tuple[str, ...] = (
    "terminase large subunit",
    "large terminase",
    "large subunit terminase",
    "terminase, large subunit",
    "dna packaging terminase large subunit",
    "dna terminase large subunit",
    "large terminase subunit",
    "terl",                  # catches "TerL", "TerL protein", etc.
)

# ---- Mechanism 2: Exact match on "product" qualifier ----
# Checks if the product annotation is EXACTLY one of these strings
# (case-sensitive). Required to catch Kayvirus group annotation "Ter".
EXACT_NAMES: frozenset[str] = frozenset({
    "Ter",    # Kayvirus group (NC_047722–727 and EU418428): TerL is 605 aa,
    "ter",    # annotated as product="Ter" — non-standard 3-letter abbreviation.
})


# ===========================================================================
# Data container
# ===========================================================================

@dataclass
class TerLRecord:
    """Container for an extracted TerL sequence and detection metadata."""
    accession:         str
    organism:          str
    product:           str
    sequence:          str
    detection_method:  str    # "keyword_match" or "exact_product_match"

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def fasta_header(self) -> str:
        """
        Tip label, made Newick/IQ-TREE-safe: only [A-Za-z0-9_.] survive.

        IQ-TREE and Newick break on spaces, parentheses, colons, commas,
        semicolons, and quotes in tip labels, so the organism is sanitised and
        the (constant) product is dropped from the label — it adds nothing to a
        tree tip and carries spaces. Format: {accession}_{organism_sanitised},
        trimmed to a readable length. Accession alone keeps every tip unique.
        """
        import re as _re
        org = self.organism
        for prefix in ("Staphylococcus phage ", "Staphylococcus virus ",
                       "Staphylococcus "):
            if org.startswith(prefix):
                org = org[len(prefix):]
                break
        org_safe = _re.sub(r"[^A-Za-z0-9]+", "_", org).strip("_")[:32]
        label = f"{self.accession}_{org_safe}" if org_safe else self.accession
        return label.strip("_")


# ===========================================================================
# Detection logic
# ===========================================================================

def _clean_accession(record_id: str) -> str:
    """
    Normalise a record id to a clean accession.version. Pharokka/older GenBank
    records can carry a compound id like 'gi|2204821279|gb|CP062427.1|', which
    breaks Newick tip labels; this extracts 'CP062427.1'. Clean ids pass through.
    """
    rid = (record_id or "").strip()
    if "|" not in rid:
        return rid
    parts = rid.split("|")
    for tag in ("gb", "ref", "emb", "dbj", "tpg", "tpe", "tpd"):
        if tag in parts:
            i = parts.index(tag)
            if i + 1 < len(parts) and parts[i + 1]:
                return parts[i + 1]
    toks = [p for p in parts if p]
    return toks[-1] if toks else rid


def find_terl_in_record(record: SeqRecord,
                        min_length_aa: int = 0) -> Optional[TerLRecord]:
    """
    Search all CDS features in a SeqRecord for the TerL sequence.

    Applies Mechanism 1 (keyword match) then Mechanism 2 (exact product match)
    to each CDS feature, collects ALL valid matches (with a non-empty
    "translation"), and returns the LONGEST.

    v3.0 fix (correctness, not cosmetic): the TerL gene is frequently SPLIT into
    several ORFs annotated "terminase large subunit 1/2/3" (common in B. subtilis
    phages; also introns/HNH insertions). The previous code returned the FIRST
    match, which was often a 130-200 aa fragment while the ~400-650 aa full-length
    ORF sat later in the record — producing wrong branch positions in the TerL
    phylogeny. Selecting the longest match recovers the full-length TerL. Genomes
    with a single TerL (e.g. all Staph Kayvirus) are unaffected: longest == first,
    so this is byte-identical for single-TerL hosts by construction.

    Parameters
    ----------
    record : SeqRecord
        A Biopython SeqRecord from a GenBank flat file.
    min_length_aa : int, optional
        Explicit fragment gate. If > 0, a genome whose longest TerL match is
        shorter than this is treated as NOT FOUND (excluded from the phylogeny,
        like the Staph Portland/436A1 exclusions). Default 0 = OFF. This is a
        deliberate knob, NOT a silent default: a length threshold can be wrong for
        a host with a genuinely short TerL, so the operator turns it on knowingly.

    Returns
    -------
    TerLRecord or None
        Longest TerL match (subject to min_length_aa); None if no TerL CDS is
        found, all matches lack a translation, or the longest is below the gate.
    """
    from phagecore.genbank_io import resolve_organism
    organism = resolve_organism(record, getattr(record, "_source_path", None))

    candidates: list[TerLRecord] = []
    for feature in record.features:
        if feature.type != "CDS":
            continue

        product_text = " ".join(feature.qualifiers.get("product", [])).lower()
        product_raw  = feature.qualifiers.get("product", ["terminase large subunit"])[0]

        # Mechanism 1: keyword substring match
        if any(kw in product_text for kw in TERL_KEYWORDS):
            method = "keyword_match"
        # Mechanism 2: exact product name match (catches Kayvirus "Ter")
        elif product_raw in EXACT_NAMES:
            method = "exact_product_match"
        else:
            continue

        # Verify translation qualifier exists
        translation = feature.qualifiers.get("translation", [])
        if not translation or not translation[0]:
            log.warning(
                f"  {record.id}: TerL match ('{product_raw}', {method}) "
                f"has no 'translation' qualifier — skipped"
            )
            continue

        candidates.append(TerLRecord(
            accession=_clean_accession(record.id),
            organism=organism,
            product=product_raw,
            sequence=translation[0],
            detection_method=method,
        ))

    if not candidates:
        return None

    # Longest match wins (first-longest on ties: max() returns the first maximal
    # element in feature order, so selection is deterministic).
    best = max(candidates, key=lambda c: c.length)

    if len(candidates) > 1:
        log.info(
            f"  {record.id}: {len(candidates)} TerL matches "
            f"({', '.join(str(c.length) for c in candidates)} aa); "
            f"selected longest = {best.length} aa"
        )

    if min_length_aa and best.length < min_length_aa:
        log.info(
            f"  {record.id}: longest TerL {best.length} aa < gate "
            f"{min_length_aa} aa — excluded as fragment"
        )
        return None

    return best


def parse_file(gb_path: Path, min_length_aa: int = 0) -> Optional[TerLRecord]:
    """Parse one GenBank file and extract its (longest) TerL sequence."""
    try:
        record = SeqIO.read(str(gb_path), "genbank")
        return find_terl_in_record(record, min_length_aa=min_length_aa)
    except ValueError as exc:
        log.warning(f"  Skipped '{gb_path.name}': {exc}")
    except Exception as exc:
        log.warning(f"  Skipped '{gb_path.name}': unexpected error — {exc}")
    return None


# ===========================================================================
# Batch processing and I/O
# ===========================================================================

def run(input_dir: Path, output_path: Path,
        min_terl_aa: int = 0) -> list[TerLRecord]:
    """
    Process all GenBank files and write combined multi-FASTA for MAFFT/MEGA.

    Genomes without a detectable TerL (e.g. Staph Portland MT926124 and
    vB_SauP-436A1 MN150710) are excluded from the FASTA output and reported in
    the terminal log, consistent with their exclusion from Figure 1. With
    min_terl_aa > 0, genomes whose longest TerL is below the gate are excluded
    the same way (reported as fragments).

    Parameters
    ----------
    input_dir : Path
        Directory containing GenBank flat files.
    output_path : Path
        Output multi-FASTA file (.faa).
    min_terl_aa : int, optional
        Explicit fragment gate passed to find_terl_in_record. Default 0 = OFF.
    """
    gb_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in VALID_EXTENSIONS
    )

    if not gb_files:
        log.error(f"No GenBank files found in '{input_dir}'.")
        sys.exit(1)

    log.info(f"Found {len(gb_files)} GenBank file(s) in '{input_dir}'")
    if min_terl_aa:
        log.info(f"TerL fragment gate ON: excluding TerL shorter than {min_terl_aa} aa")

    found:     list[TerLRecord] = []
    not_found: list[str]        = []

    for gb in gb_files:
        log.info(f"  Processing {gb.name}")
        terl = parse_file(gb, min_length_aa=min_terl_aa)
        if terl:
            found.append(terl)
            log.info(
                f"    ✓  {terl.accession:<16}  {terl.length:>4} aa  "
                f"[{terl.detection_method}]  product: '{terl.product}'"
            )
        else:
            not_found.append(gb.stem)   # accession stem
            log.info(f"    ✗  {gb.stem}  — no TerL detected (correctly excluded)")

    if not found:
        log.error("No TerL sequences extracted. Check GenBank annotations.")
        sys.exit(1)

    # Write combined multi-FASTA sorted by accession
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        for t in sorted(found, key=lambda x: x.accession):
            fh.write(f">{t.fasta_header}\n{t.sequence}\n")

    log.info(f"\nOutput multi-FASTA written: '{output_path}'")
    log.info(f"  Sequences included : {len(found)}")
    log.info(f"  Genomes excluded   : {len(not_found)}")

    if not_found:
        log.info("  Excluded genomes (no annotated TerL):")
        for acc in not_found:
            log.info(f"    • {acc}")
        log.info(
            "  Methods statement: 'Staphylococcus phage Portland (MT926124) "
            "and vB_SauP-436A1 (MN150710) were excluded from phylogenetic "
            "analysis due to the absence of annotated terminase large subunit "
            "sequences, consistent with their atypical small genome sizes "
            "(<20 kb) relative to the remaining dataset.'"
        )

    return found


# ===========================================================================
# CLI
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="S3_terl_extractor.py",
        description=(
            "Extract TerL sequences from Staphylococcus phage GenBank files. "
            "Applies two detection mechanisms: (1) keyword match and "
            "(2) exact product match for 'Ter'/'ter' (Kayvirus group fix). "
            "Outputs combined multi-FASTA for MAFFT + MEGA 12.1.2 (Figure 1)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MANUSCRIPT OUTPUT
-----------------
  FIGURE 1 — Maximum Likelihood Phylogenetic Tree

  Pipeline for Figure 1:
    Step 1: python S3_terl_extractor.py -i GenBank -o results\\TerL_combined.faa
    Step 2: Align on MAFFT web server (https://mafft.cbrc.jp/alignment/server/)
            Upload TerL_combined.faa, set strategy = L-INS-i, Submit,
            save result as TerL_aligned.faa
    Step 3: Open results\\TerL_aligned.faa in MEGA 12.1.2
            Phylogeny > ML > LG+G+I > 1000 bootstrap > 80%% partial deletion
            Outgroup: Staphylococcus phage EW (NC_007056.1)

Detection mechanisms:
  Mechanism 1: keyword match in product qualifier (standard annotations)
  Mechanism 2: exact match for product="Ter"/"ter" (Kayvirus group fix)

Genomes excluded from FASTA (no annotated TerL; consistent with Figure 1):
  • Portland (MT926124)     — 17,471 bp, no TerL annotation
  • vB_SauP-436A1 (MN150710) — 18,028 bp, no TerL annotation

Examples (Windows Command Prompt)
----------------------------------
  python S3_terl_extractor.py -i GenBank -o results\\TerL_combined.faa
        """,
    )
    parser.add_argument(
        "--input_dir", "-i",
        type=Path, required=True, metavar="DIR",
        help="Directory containing GenBank files (.gb / .gbk / .gbff)",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path, default=Path("TerL_combined.faa"), metavar="FILE",
        help="Output multi-FASTA file. Default: TerL_combined.faa",
    )
    parser.add_argument(
        "--min-terl-aa",
        type=int, default=0, metavar="AA",
        help="OPTIONAL fragment gate: exclude any genome whose longest TerL is "
             "shorter than AA residues (like the Portland/436A1 exclusions). "
             "Default 0 = OFF. Turn on deliberately — a length threshold can be "
             "wrong for a host with a genuinely short TerL. Longest-match "
             "selection is always on and needs no flag.",
    )
    return parser


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    args  = build_parser().parse_args()
    found = run(args.input_dir, args.output, min_terl_aa=args.min_terl_aa)

    sep = "=" * 80
    kw_count    = sum(1 for t in found if t.detection_method == "keyword_match")
    exact_count = sum(1 for t in found if t.detection_method == "exact_product_match")

    print(f"\n{sep}")
    print("FIGURE 1 INPUT — TerL EXTRACTION SUMMARY")
    print(sep)
    print(f"{'Accession':<20} {'Organism':<36} {'Length':>6}  {'Method':<22}  {'Product'}")
    print("-" * 80)
    for t in sorted(found, key=lambda x: x.accession):
        print(f"{t.accession:<20} {t.organism[:35]:<36} {t.length:>6}  "
              f"{t.detection_method:<22}  {t.product}")
    print(sep)
    print(f"  Total sequences extracted : {len(found)}")
    print(f"  Via keyword_match         : {kw_count}  (standard annotations)")
    print(f"  Via exact_product_match   : {exact_count}  (Kayvirus 'Ter' fix)")
    print(f"\n  Output FASTA   : {args.output}")
    print(f"\nNEXT STEP:")
    print(f"  1. Align on MAFFT web: https://mafft.cbrc.jp/alignment/server/")
    print(f"     Upload {args.output}, set strategy = L-INS-i, Submit,")
    print(f"     save result as TerL_aligned.faa")
    print(f"  2. Open TerL_aligned.faa in MEGA 12.1.2 -> ML tree")
    print(f"     (LG+G+I, 1000 bootstrap, outgroup EW NC_007056.1)")
