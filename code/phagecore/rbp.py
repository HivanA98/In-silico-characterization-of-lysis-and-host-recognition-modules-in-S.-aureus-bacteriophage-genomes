#!/usr/bin/env python3
# phagecore.rbp (v3.0) — S5 RBP + depolymerase ENGINE.
# Universal detection lives here once; host CONTENT arrives via configure(profile).
# The standalone S5_rbp_extractor.py is now a thin CLI over this module.
"""
S5_rbp_extractor.py
===================
Receptor-Binding Protein (RBP) / Tail Fiber Extractor and Classifier

==============================================================================
MANUSCRIPT CONTRIBUTION
  Extends the host-recognition data of TABLE 2 (S2_holin_tailfiber_annotation.py)
  from a binary presence/absence flag to sequence-level classification:

  S2 output  →  Tail_Fiber_RBP_Present  (Yes / No)
  S5 output  →  Per genome:
                  • Classified multi-FASTA  → InterPro domain analysis
                  • Multi-FASTA             → MAFFT (G-INS-i) → MEGA RBP phylogeny
                  • Audit CSV with subtype, length, confidence, product source,
                    length flag, and InterPro_Domain column to fill after the run

  THREE DESIGN GOALS
  ──────────────────
  Goal 1 — Current 22-phage Journal 2 dataset:
    S2 detected RBP in 14 of 22 phages; 8 returned "No". S5 uses an expanded
    keyword set (Mechanism 2) and a length heuristic (Mechanism 3, opt-in) to
    distinguish:
      (a) ANNOTATION GAP (rescued)   — an RBP-like protein IS annotated under
          a non-standard product name the S2 keyword set did not capture.
      (b) ANNOTATION GAP (confirmed) — no RBP-like annotation found by any
          mechanism; requires tBLASTn / HHpred / Foldseek as next step.

  Goal 2 — Expanded dataset (target: 150 phages, ~40 with RBP):
    Annotation quality across GenBank records is variable. Three output
    quality controls address this for large-scale runs:
      • Product_Source column  — flags when the keyword was found in note/gene
        rather than the product qualifier; confidence is downgraded one tier
        when product = "hypothetical protein" regardless of mechanism.
      • Length_Flag column     — marks short_fragment (<100 aa) and very_long
        (>2000 aa) entries for downstream filtering before MAFFT alignment.
      • --deduplicate flag     — when the same organism appears under multiple
        accessions, keeps the RefSeq (NC_) record preferentially.

  Goal 3 — RBP phylogeny vs TerL phylogeny comparison:
    The confirmed multi-FASTA feeds MAFFT + MEGA 12.1.2 to produce an RBP
    phylogenetic tree. Topological incongruence between the RBP tree and the
    TerL tree (Figure 1, Journal 2) is evidence of RBP horizontal gene
    transfer — a key finding that supports the modularity argument underpinning
    the S2 proposal.

  WHAT S5 CANNOT DO
  ─────────────────
  S5 is annotation-dependent. If GenBank carries no RBP-like product name or
  qualifier for a given CDS, S5 will not detect it regardless of keyword
  breadth. For annotation-gap genomes the next steps are manual:
    • tBLASTn against characterised RBP references (Yehl et al., 2019)
    • HHpred / Foldseek for structural homology detection
  S5 flags these genomes explicitly so manual effort is targeted correctly.

==============================================================================

QUALITY CONTROL LOGIC (four layers)
--------------------------------------
Layer 1 — Structural exclusion filter (before all mechanisms):
  CDS whose PRODUCT qualifier names a known non-RBP structural component (tail
  sheath, tape measure protein, baseplate assembly wedge, etc.) are excluded
  before any keyword check. The filter is applied to the product field only,
  not to note/gene/function text, so a real RBP whose note merely mentions an
  excluded term is not falsely discarded. This prevents the most common
  false-positive class without creating false negatives.

Layer 2 — Product_Source tracking:
  get_qualifier_text() searches product + gene + note + function qualifiers.
  When a keyword match is found, the script records WHERE the keyword came from:
    "product_qualifier"  — keyword in the product field (more reliable)
    "note/gene_qualifier" — keyword NOT in product field; found in gene/note
  If product = "hypothetical protein" AND source = "note/gene_qualifier",
  confidence is downgraded one tier:
    confirmed → high
    high      → putative
    putative  → putative  (unchanged)
  This corrects the S5-v1 behaviour where "hypothetical protein" entries were
  classified as "confirmed" solely because the note field contained "tail fiber".

Layer 3 — Length flags (not filters; analyst decides):
  Every candidate receives a Length_Flag:
    normal         :  80 aa <= length <= 2000 aa
    short_fragment : length < 80 aa
    very_long      : length > 2000 aa
  Note on short sequences: legitimate short tail fiber subunits exist (e.g.,
  tail fiber proximal subunit ~110-135 aa in Siphoviridae-type phages, and
  verified Twort/EW tail fibers at 110-121 aa). The threshold is 80 aa (not
  100 aa) to avoid excluding these biologically real proteins. Sequences below
  80 aa are unlikely to fold independently as host-recognition domains.
  --min-length controls which Length_Flag values are written to FASTA output:
  default 80 aa; set to 200 aa for a phylogeny-ready filtered FASTA.

Layer 4 — Deduplication (opt-in, --deduplicate):
  When two or more GenBank records share the same organism name, only one is
  processed. Selection priority: RefSeq (NC_/NZ_) accession > lowest version
  number (e.g., .1 before .2). The discarded accession is logged at WARNING.
  Purpose: prevents phiSA12 (AB903967.1 and NC_023573.1) and similar cases
  from inflating sequence counts in the FASTA and CSV.

Detection Logic (Three Mechanisms)
-------------------------------------
Mechanism 1 — Primary RBP keywords (explicit annotation):
  Standard product names widely used in GenBank. A match here yields
  confidence "confirmed" (tail_fiber, tail_spike, rbp_annotated) or
  "high" (host_recognition subtype), subject to Layer 2 downgrade.

Mechanism 2 — Extended RBP keywords (annotation-gap rescue):
  Broader terms capturing non-standard but plausible RBP annotations.
  A match here yields confidence "putative"; InterPro validation required.

Mechanism 3 — Length heuristic (opt-in, --include-length-candidates):
  Proteins >= RBP_LENGTH_HEURISTIC_MIN_AA (600 aa) that carry any
  tail/baseplate context term but do not match Mechanism 1 or 2 are
  flagged as "length_candidate". Disabled by default; high false-positive
  rate requires manual InterPro review.

RBP Subtype Classification
---------------------------
  tail_fiber        : "tail fiber" / "tail fibre" (including long/short/lateral)
  tail_spike        : "tail spike" / "tailspike"
  rbp_annotated     : "receptor binding protein" / "rbp"
  host_recognition  : "host recognition", "host specificity", "adsorption",
                      "anti-receptor", "baseplate receptor-binding"
  putative_rbp      : Mechanism 2 match only, or Mechanism 1 confidence-downgraded
  length_candidate  : Mechanism 3 match only (opt-in)

Multiple RBP per genome
-------------------------
  Some phages encode multiple tail fiber or RBP genes. S5 collects ALL
  matching CDS per genome (unlike S3 which returns the first TerL). The
  audit CSV has one row per candidate; the FASTA includes all entries that
  pass the --min-length threshold.

Known annotation-gap accessions (22-phage Journal 2 dataset)
--------------------------------------------------------------
  The 8 accessions returning Tail_Fiber_RBP_Present = "No" in S2:
    KY779849  (qdsa002)     — Twortvirinae
    MN336261  (Sb1_8383)    — Sb-1 lineage
    MN336262  (Sb1M_6168)   — Sb-1 lineage
    MN336263  (Sb1M_9832)   — Sb-1 lineage
    NC_047724 (676Z)        — Kayvirus group
    NC_047725 (Fi200W)      — Kayvirus group
    NC_047726 (MSA6)        — Kayvirus group
    NC_047727 (P4W)         — Kayvirus group

Downstream Workflow (RBP Phylogeny)
--------------------------------------
  Step 1 — This script:
    python S5_rbp_extractor.py -i GenBank ^
        -o results\\rbp_candidates.faa --csv results\\rbp_audit.csv

  Step 2 — Filter FASTA before InterPro:
    Keep only rows where Confidence IN (confirmed, high) AND
    Length_Flag = normal. Remove putative and length_candidate unless
    InterPro confirms an RBP domain. This is the analyst's decision;
    the CSV provides all information needed to apply these filters.

  Step 3 — InterPro domain analysis:
    Submit filtered FASTA to https://www.ebi.ac.uk/interpro/search/sequence/
    Databases: Pfam, CDD, SUPERFAMILY, Gene3D.
    Fill the InterPro_Domain column of rbp_audit.csv from the result.
    Expected RBP domains: DUF3751, Phage_fiber, Tail_spike_N, baseplate_J,
    Receptor_bind, TSP_N, or structural folds (beta-propeller, beta-helix).

  Step 4 — Multiple sequence alignment (MAFFT web server):
    Open https://mafft.cbrc.jp/alignment/server/
    Strategy: G-INS-i (globally alignable sequences, moderate divergence).
    Use L-INS-i only for highly similar, short RBP sets.
    Save result as results\\rbp_aligned.faa.

  Step 5 — RBP Phylogenetic tree (MEGA 12.1.2):
    Open rbp_aligned.faa > Phylogeny > Maximum Likelihood
      Substitution model   : LG+G+I (or best BIC model)
      Bootstrap replicates : 1000
      Partial deletion     : 80% site-coverage cutoff
      Outgroup             : most divergent RBP lineage (EW or SA13)
    Compare topology against TerL tree (Figure 1, Journal 2):
      Congruent topology   → RBP co-evolves with TerL core genome
      Incongruent topology → RBP horizontal gene transfer; supports
                             modularity argument for S2 proposal

References
----------
Yehl K, Lemire S, Yang AC, Ando H, Mimee M et al. (2019). Engineering phage
host-range and suppressing bacterial resistance through phage tail fiber
mutagenesis. Cell 179(2):459-469. doi:10.1016/j.cell.2019.09.015

Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I,
Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009). Biopython: freely
available Python tools for computational molecular biology and bioinformatics.
Bioinformatics 25(11):1422-1423. doi:10.1093/bioinformatics/btp163

Tested Environment (Windows)
-----------------------------
    OS         : Windows 10 / 11
    Python     : 3.12.10
    biopython  : 1.87

Installation (Command Prompt / PowerShell)
------------------------------------------
    pip install biopython==1.87

Usage (Windows Command Prompt)
-------------------------------
    # Standard run — Mechanisms 1 and 2, min length 80 aa:
    python S5_rbp_extractor.py -i GenBank ^
        -o results\\rbp_candidates.faa --csv results\\rbp_audit.csv

    # Phylogeny-ready FASTA (min 200 aa, deduplication on):
    python S5_rbp_extractor.py -i GenBank ^
        -o results\\rbp_phylo.faa --csv results\\rbp_audit.csv ^
        --min-length 200 --deduplicate

    # Full sensitivity (all three mechanisms, no length filter):
    python S5_rbp_extractor.py -i GenBank ^
        -o results\\rbp_all.faa --csv results\\rbp_audit.csv ^
        --min-length 0 --include-length-candidates
"""

import argparse
import csv
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
INTERPRO_MAX_SEQUENCES: int = 100
COLABFOLD_MAX_AA: int = 1200   # fold targets above this are domain-split first (Colab OOM)
LARGE_CARRIER_SUBMIT_AA: int = 700   # always submit carriers this large to InterPro,
                                     # regardless of confidence tier (see the submission filter)

# ---- Mechanism 1: Primary RBP keywords (explicit product annotation) ----
# Matched case-insensitively as substrings in the full qualifier text.
# Yields confidence "confirmed" or "high" before Layer 2 downgrade.
PRIMARY_RBP_KEYWORDS: tuple[str, ...] = (
    "tail fiber",
    "tail fibre",
    "receptor binding protein",
    "receptor-binding protein",
    "tail spike",
    "tailspike",
    "tail-spike",
    "host recognition protein",
    "adsorption protein",
    "host specificity protein",
    "baseplate receptor-binding",
    "anti-receptor",
    "rbp",
)

# ---- Mechanism 2: Extended RBP keywords (annotation-gap rescue) ----
# Broader terms for non-standard RBP product names.
# Yields confidence "putative" regardless of Layer 2 — no downgrade needed.
EXTENDED_RBP_KEYWORDS: tuple[str, ...] = (
    "distal tail protein",
    "tail tip protein",
    "receptor recognizing",
    "host range determinant",
    "lateral tail fiber",
    "baseplate spike",
    "tail knob",
    "tail needle",
    "baseplate hub",
    "long tail fiber",
    "short tail fiber",
    "dit protein",
    "fiber protein",
    "tail associated receptor",
)

# ===========================================================================
# v2 ADDITION — DEPOLYMERASE AXIS (the MSc biofilm-direction variable)
# ---------------------------------------------------------------------------
# Rationale: S5-v1 detects the HOST-RECOGNITION role (tail fiber / spike / RBP =
# what the phage BINDS). The MSc hypothesis predicts on the DEPOLYMERASE
# genotype (what the phage ENZYMATICALLY DEGRADES = capsule / EPS / biofilm
# matrix). These are different functional classes that frequently co-occur in
# ONE protein (a tail spike often CARRIES a depolymerase domain), so the
# depolymerase signal is recorded on a SEPARATE axis (depolymerase_signal +
# module), never overwriting the structural subtype.
#
# IMPORTANT SCOPE NOTE for S. aureus: its biofilm matrix is PNAG/PIA
# (ica-dependent poly-N-acetylglucosamine) or proteinaceous. The PNAG-relevant
# enzymatic class is beta-N-acetylglucosaminidase (a glycoside hydrolase).
# Most catalogued phage depolymerases (pectin/pectate lyase, sialidase) come
# from CAPSULED Gram-negatives (Klebsiella); they are included for the Phase-2
# multi-host scaffold but are rarer in Staphylococcus. Detecting FEW
# depolymerases in the S. aureus set is therefore a real biological result, not
# a failure of S5 — record it, do not force matches.
# ===========================================================================

# Each entry maps a product/note substring -> the enzyme class label recorded
# in depolymerase_signal. Order matters only for which label is reported first.
DEPOLYMERASE_KEYWORDS: tuple[tuple[str, str], ...] = (
    ("depolymerase",                      "depolymerase_generic"),
    ("pectin lyase",                      "pectin_lyase"),
    ("pectate lyase",                     "pectate_lyase"),
    ("pectin/pectate lyase",              "pectin_lyase"),
    ("parallel beta-helix",               "beta_helix_lyase"),
    ("parallel beta helix",               "beta_helix_lyase"),
    ("sialidase",                         "sialidase"),
    ("neuraminidase",                     "sialidase"),
    ("levanase",                          "levanase"),
    ("dextranase",                        "dextranase"),
    ("hyaluronidase",                     "hyaluronidase"),
    ("rhamnosidase",                      "rhamnosidase"),
    ("alginate lyase",                    "alginate_lyase"),
    ("glucosaminidase",                   "glucosaminidase"),       # PNAG-relevant
    ("n-acetylglucosaminidase",           "glucosaminidase"),       # PNAG-relevant
    ("acetylgalactosaminidase",           "glycoside_hydrolase"),
    ("glycoside hydrolase",               "glycoside_hydrolase"),
    ("glycosyl hydrolase",                "glycoside_hydrolase"),
    ("glycosidase",                       "glycoside_hydrolase"),
    ("polysaccharide lyase",              "polysaccharide_lyase"),
    ("sgnh hydrolase",                    "sgnh_esterase"),
    ("sgnh",                              "sgnh_esterase"),
    ("carbohydrate esterase",             "esterase"),
    ("acetylesterase",                    "esterase"),
    ("tail spike",                        "tailspike_depolymerase"),
    ("tailspike",                         "tailspike_depolymerase"),
    ("tail-spike",                        "tailspike_depolymerase"),
)

# Length above which a tail/structural protein is almost certainly multi-domain
# and very likely carries a depolymerase domain even when annotated generically
# (e.g. Machias 3084 aa, Madawaska 2706 aa, vB_StaM_PB50 2781 aa). These get a
# multidomain_flag that MANDATES per-domain InterPro/structural analysis.
MULTIDOMAIN_CARRIER_MIN_AA: int = 1200


# ---- Exclusion filter: structural tail components that are NOT RBP ----
# Applied BEFORE all three mechanisms.
STRUCTURAL_EXCLUSIONS: tuple[str, ...] = (
    "tail sheath",
    "tail tube",
    "tape measure",
    "major tail protein",
    "tail terminator",
    "tail length",
    "baseplate assembly wedge",
    "baseplate wedge",
    "tail completion",
    "tail lysozyme",
    "tail associated lysin",
    "phage_lysozyme2",
    "contractile tail",
    "tail baseplate protein",
    "tail sheath stabilizer",
)

# "adhesin" removed from EXTENDED_RBP_KEYWORDS: the term is too broad and
# generates false positives (e.g., "adhesin/Ig-like protein" at 170 aa in
# Maine is not a canonical RBP). Phage-specific adhesins that are genuine
# RBPs will typically carry an additional "tail" or "fiber" qualifier and
# will be captured by PRIMARY_RBP_KEYWORDS.

# ---- Mechanism 3: Length heuristic parameters (opt-in only) ----
RBP_LENGTH_HEURISTIC_MIN_AA: int = 600
RBP_LENGTH_HEURISTIC_CONTEXT_TERMS: tuple[str, ...] = (
    "tail", "baseplate", "structural protein", "fiber",
)

# ---- Length flag thresholds ----
# short_fragment: likely a tail fiber accessory subunit or annotation fragment.
# Note: 80 aa (not 100 aa) is the threshold because legitimate short tail
# fiber proteins exist at 100-135 aa (e.g., EW 110 aa, Twort 121 aa,
# tail fiber proximal subunit ~125-127 aa in multiple Siphoviridae-type phages).
# Sequences 80-100 aa are borderline; InterPro validation is recommended.
LENGTH_SHORT_FRAGMENT_MAX_AA: int = 80
LENGTH_VERY_LONG_MIN_AA: int = 2000

# ---- RefSeq prefix strings — used by deduplication to prefer NC_/NZ_ ----
REFSEQ_PREFIXES: tuple[str, ...] = ("NC_", "NZ_", "NG_", "NW_", "NR_")

# ---- Host CONTENT injected from the profile via configure() (v3.0) ----
# Was the module constant _S2_GAP; now profile.s2_annotation_gap_accessions.
_S2_GAP: frozenset = frozenset()


# ===========================================================================
# Data container
# ===========================================================================

@dataclass
class RBPCandidate:
    """Container for one RBP / tail fiber candidate CDS and its metadata."""
    protein_id:      str
    accession:       str
    organism:        str
    product:         str
    sequence:        str
    subtype:         str    # tail_fiber | tail_spike | rbp_annotated |
                            # host_recognition | putative_rbp | length_candidate
    confidence:      str    # confirmed | high | putative | length_candidate
    mechanism:       str    # primary | extended | length_heuristic
    product_source:  str    # product_qualifier | note/gene_qualifier
    length_flag:     str    # normal | short_fragment | very_long
    was_s2_no:       bool
    evidence:        str
    interpro_domain: str = ""   # filled after InterPro run
    # ---- v2 ADDITIONS (depolymerase axis + InterPro reconciliation) ----
    depolymerase_signal: str = "none"   # none | <enzyme class> (annotation-level)
    module:              str = "host_recognition"  # host_recognition | depolymerase | both
    multidomain_flag:    str = ""       # multidomain_carrier (>MULTIDOMAIN_CARRIER_MIN_AA)
    deposcope_score:     float = -1.0   # filled by --deposcope: ESM-2 probability 0-1 (-1 = not run)
    deposcope_call:      str = ""       # depolymerase | non_depolymerase | not_evaluated
    method_agreement:    str = ""       # agree_positive | agree_negative | S5_only | DepoScope_only
    interpro_verdict:    str = ""       # filled by --interpro: rbp | depolymerase |
                                        # structural_tail | no_domain | uncertain

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def fasta_header(self) -> str:
        """
        FASTA header format:
          {protein_id}|{accession}|{organism_no_spaces}|
          subtype={subtype}|conf={confidence}|src={product_source}|{product}

        The src tag allows downstream filtering of note/gene_qualifier
        detections before MAFFT alignment.
        """
        org  = self.organism.replace(" ", "_")
        prod = self.product.replace("|", "/")
        return (
            f"{self.protein_id}|{self.accession}|{org}|"
            f"subtype={self.subtype}|conf={self.confidence}|"
            f"src={self.product_source}|{prod}"
        )


# ===========================================================================
# Quality control helpers
# ===========================================================================

def clean_accession(record_id: str) -> str:
    """
    v2 fix — normalise a record id to a clean accession.version.

    Pharokka / older GenBank records can carry a compound NCBI identifier such
    as 'gi|2204821279|gb|CP062427.1|'. Using that verbatim in a FASTA header is
    ugly and, combined with a missing protein_id, produces colliding headers.
    This extracts 'CP062427.1' from such strings; clean ids pass through.
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


def robust_protein_id(feature, accession: str) -> str:
    """
    v2 fix — never return a non-unique 'unknown'.

    Order: protein_id → locus_tag → '{accession_stem}_cds@{start}'.
    The location-based fallback is unique within a genome because each CDS has a
    distinct start coordinate, so headers from Pharokka records (which lack
    protein_id) no longer collide.
    """
    for key in ("protein_id", "locus_tag"):
        vals = feature.qualifiers.get(key, [])
        if vals and vals[0]:
            return vals[0]
    try:
        start = int(feature.location.start)
    except Exception:
        start = 0
    return f"{accession.split('.')[0]}_cds@{start}"


def get_qualifier_text(feature) -> tuple[str, str, str]:
    """
    Concatenate CDS qualifier values for keyword matching.

    Returns
    -------
    tuple[str, str, str]
        (full_text_lower, product_text_lower, product_raw)
        full_text_lower : product + gene + note + function, lowercased
        product_text_lower : product qualifier only, lowercased
        product_raw : product qualifier as-is (original case)
    """
    product_vals  = feature.qualifiers.get("product", [])
    product_raw   = product_vals[0] if product_vals else "hypothetical protein"
    product_lower = product_raw.lower()

    extra_parts = []
    for key in ("gene", "note", "function"):
        extra_parts.extend(feature.qualifiers.get(key, []))
    extra_lower = " ".join(extra_parts).lower()

    full_lower = f"{product_lower} {extra_lower}".strip()
    return full_lower, product_lower, product_raw


def resolve_product_source(product_lower: str, matched_keyword: str) -> str:
    """
    Determine whether the matched keyword came from the product qualifier.

    Returns "product_qualifier" if the keyword is in the product field,
    "note/gene_qualifier" otherwise. This drives the Layer 2 confidence
    downgrade when product = "hypothetical protein".
    """
    if matched_keyword in product_lower:
        return "product_qualifier"
    return "note/gene_qualifier"


def apply_confidence_downgrade(
    confidence: str,
    subtype: str,
    product_source: str,
    product_lower: str,
) -> tuple[str, str]:
    """
    Apply Layer 2 confidence downgrade.

    Rule: if product = "hypothetical protein" AND keyword was found in
    note/gene (not in the product field itself), downgrade one tier:
      confirmed → high
      high      → putative
      putative  → putative  (floor; no further downgrade)
    The subtype remains unchanged; it reflects what the note/gene says.

    Returns
    -------
    tuple[str, str]
        (adjusted_confidence, adjusted_subtype)
    """
    is_hypo = product_lower.strip() == "hypothetical protein"
    is_note_source = product_source == "note/gene_qualifier"

    if is_hypo and is_note_source:
        if confidence == "confirmed":
            return "high", subtype
        if confidence == "high":
            return "putative", subtype
        # putative → stays putative
    return confidence, subtype


def assign_length_flag(length: int) -> str:
    """
    Assign a length flag to a candidate sequence.

    Returns
    -------
    str
        "short_fragment" (<= LENGTH_SHORT_FRAGMENT_MAX_AA),
        "very_long"      (>= LENGTH_VERY_LONG_MIN_AA),
        "normal"         (all others)
    """
    if length < LENGTH_SHORT_FRAGMENT_MAX_AA:
        return "short_fragment"
    if length >= LENGTH_VERY_LONG_MIN_AA:
        return "very_long"
    return "normal"


def is_structural_exclusion(text: str) -> bool:
    """
    Return True if the annotation matches a known non-RBP structural component.
    Applied before Mechanism 1–3 to eliminate the most common false-positive class.
    """
    return any(excl in text for excl in STRUCTURAL_EXCLUSIONS)


def detect_depolymerase(full_text: str, length: int) -> tuple[str, str]:
    """
    v2 — Scan for a depolymerase enzymatic signal on the SEPARATE depolymerase
    axis. This never overwrites the structural subtype; a tail spike that also
    carries a depolymerase keyword is recorded as subtype=tail_spike AND
    depolymerase_signal=<class>.

    Returns
    -------
    tuple[str, str]
        (depolymerase_signal, multidomain_flag)
        depolymerase_signal : "none" or the enzyme-class label.
        multidomain_flag    : "multidomain_carrier" if length is large enough to
                              almost certainly be multi-domain, else "".

    Caveat encoded here: annotation-level detection only. A genomic depolymerase
    keyword/domain is a GENOTYPE signal, not proof of active matrix degradation;
    the MSc biofilm phenotype must be measured in the wet lab. The multidomain
    flag deliberately fires on LENGTH alone (independent of keywords) because the
    giant tail proteins that carry depolymerase domains are frequently annotated
    only as "putative tail fiber".
    """
    signal = "none"
    for kw, label in DEPOLYMERASE_KEYWORDS:
        if kw in full_text:
            signal = label
            break
    flag = "multidomain_carrier" if length >= MULTIDOMAIN_CARRIER_MIN_AA else ""
    return signal, flag


def resolve_module(subtype: str, depolymerase_signal: str) -> str:
    """Combine the structural role and the enzymatic axis into one module label."""
    has_depo = depolymerase_signal not in ("none", "")
    is_recognition = subtype in (
        "tail_fiber", "tail_spike", "rbp_annotated",
        "host_recognition", "putative_rbp", "length_candidate",
    )
    if has_depo and is_recognition:
        return "both"
    if has_depo:
        return "depolymerase"
    return "host_recognition"


# ===========================================================================
# Subtype and confidence assignment
# ===========================================================================

def classify_subtype_primary(text: str) -> tuple[str, str]:
    """
    Assign RBP subtype and initial confidence for a Mechanism 1 match.

    Parameters
    ----------
    text : str
        Lowercased full qualifier text.

    Returns
    -------
    tuple[str, str]
        (subtype, confidence)  — before Layer 2 downgrade.
    """
    if any(k in text for k in (
        "tail fiber", "tail fibre", "long tail fiber",
        "short tail fiber", "lateral tail fiber",
    )):
        return "tail_fiber", "confirmed"
    if any(k in text for k in ("tail spike", "tailspike", "tail-spike")):
        return "tail_spike", "confirmed"
    if any(k in text for k in ("receptor binding protein", "receptor-binding protein")):
        return "rbp_annotated", "confirmed"
    # "rbp" alone — only classify if it is clearly a label, not embedded in a word
    if " rbp " in f" {text} ":
        return "rbp_annotated", "confirmed"
    if any(k in text for k in (
        "host recognition", "host specificity", "adsorption protein",
        "anti-receptor", "baseplate receptor-binding",
    )):
        return "host_recognition", "high"
    # Fallback for any remaining primary keyword match
    return "putative_rbp", "high"


# ===========================================================================
# Detection logic
# ===========================================================================

def scan_cds_for_rbp(
    feature,
    accession: str,
    organism: str,
    include_length_candidates: bool,
) -> Optional[RBPCandidate]:
    """
    Apply the three-mechanism detection pipeline to a single CDS feature.

    Returns an RBPCandidate if any mechanism matches, else None.

    Detection order: exclusion filter → Mechanism 1 → Mechanism 2 →
    Mechanism 3 (opt-in). Each candidate is then processed through
    Layer 2 (confidence downgrade) and Layer 3 (length flag).

    Parameters
    ----------
    feature : SeqFeature
        A Biopython CDS SeqFeature.
    accession : str
        GenBank accession string (with version, e.g., NC_047724.1).
    organism : str
        Organism name from record.annotations["organism"].
    include_length_candidates : bool
        Whether to activate Mechanism 3.

    Returns
    -------
    RBPCandidate or None
    """
    # Require a translated sequence
    translation = feature.qualifiers.get("translation", [])
    if not translation or not translation[0]:
        return None

    full_text, product_lower, product_raw = get_qualifier_text(feature)

    # Layer 1 — Structural exclusion filter.
    # Applied to the PRODUCT qualifier ONLY (not the full note/gene/function
    # text). Rationale: a genuine RBP whose note merely mentions an excluded
    # term (e.g., "located downstream of the tail tube") must not be discarded.
    # The product field is the authoritative functional annotation; if it names
    # a structural non-RBP component (tail sheath, tape measure, etc.) the CDS
    # is excluded. The previous version tested full_text and could drop real
    # RBP candidates whose note text happened to contain an excluded term.
    if is_structural_exclusion(product_lower):
        return None

    pid      = robust_protein_id(feature, accession)
    seq      = translation[0]
    acc_stem = accession.split(".")[0]

    # ---- Mechanism 1: Primary keywords ----
    for kw in PRIMARY_RBP_KEYWORDS:
        if kw in full_text:
            subtype, confidence = classify_subtype_primary(full_text)
            product_source      = resolve_product_source(product_lower, kw)
            confidence, subtype = apply_confidence_downgrade(
                confidence, subtype, product_source, product_lower
            )
            length_flag = assign_length_flag(len(seq))
            depo_signal, mdflag = detect_depolymerase(full_text, len(seq))
            return RBPCandidate(
                protein_id=pid, accession=accession, organism=organism,
                product=product_raw, sequence=seq,
                subtype=subtype, confidence=confidence, mechanism="primary",
                product_source=product_source, length_flag=length_flag,
                was_s2_no=(acc_stem in _S2_GAP),
                evidence=f"primary keyword: '{kw}'",
                depolymerase_signal=depo_signal,
                module=resolve_module(subtype, depo_signal),
                multidomain_flag=mdflag,
            )

    # ---- Mechanism 2: Extended keywords (universal + host vocabulary) ----
    for kw in tuple(EXTENDED_RBP_KEYWORDS) + _HOST_RBP_TERMS:
        if kw in full_text:
            product_source = resolve_product_source(product_lower, kw)
            length_flag    = assign_length_flag(len(seq))
            depo_signal, mdflag = detect_depolymerase(full_text, len(seq))
            return RBPCandidate(
                protein_id=pid, accession=accession, organism=organism,
                product=product_raw, sequence=seq,
                subtype="putative_rbp", confidence="putative",
                mechanism="extended",
                product_source=product_source, length_flag=length_flag,
                was_s2_no=(acc_stem in _S2_GAP),
                evidence=f"extended keyword: '{kw}'",
                depolymerase_signal=depo_signal,
                module=resolve_module("putative_rbp", depo_signal),
                multidomain_flag=mdflag,
            )

    # ---- Mechanism 2b (v2): STANDALONE DEPOLYMERASE route ----
    # A depolymerase can be a free enzyme with no RBP keyword (no tail/fiber/spike
    # in the annotation). Mechanism 1/2 would miss it. Capture it here on the
    # depolymerase axis so the thesis variable is not silently dropped.
    depo_signal, mdflag = detect_depolymerase(full_text, len(seq))
    if depo_signal != "none":
        product_source = resolve_product_source(product_lower, "")
        # find which keyword fired, for the evidence string
        fired = next((lab for kw, lab in DEPOLYMERASE_KEYWORDS if kw in full_text), depo_signal)
        return RBPCandidate(
            protein_id=pid, accession=accession, organism=organism,
            product=product_raw, sequence=seq,
            subtype="depolymerase", confidence="putative",
            mechanism="depolymerase",
            product_source=product_source,
            length_flag=assign_length_flag(len(seq)),
            was_s2_no=(acc_stem in _S2_GAP),
            evidence=f"depolymerase keyword: class='{depo_signal}'",
            depolymerase_signal=depo_signal,
            module="depolymerase",
            multidomain_flag=mdflag,
        )

    # ---- Mechanism 3: Length heuristic (opt-in) ----
    if include_length_candidates and len(seq) >= RBP_LENGTH_HEURISTIC_MIN_AA:
        if any(term in full_text for term in RBP_LENGTH_HEURISTIC_CONTEXT_TERMS):
            product_source = resolve_product_source(product_lower, "")
            depo_signal, mdflag = detect_depolymerase(full_text, len(seq))
            return RBPCandidate(
                protein_id=pid, accession=accession, organism=organism,
                product=product_raw, sequence=seq,
                subtype="length_candidate", confidence="length_candidate",
                mechanism="length_heuristic",
                product_source=product_source,
                length_flag=assign_length_flag(len(seq)),
                was_s2_no=(acc_stem in _S2_GAP),
                evidence=(
                    f"length heuristic: {len(seq)} aa >= "
                    f"{RBP_LENGTH_HEURISTIC_MIN_AA}, "
                    f"tail/baseplate context term present"
                ),
                depolymerase_signal=depo_signal,
                module=resolve_module("length_candidate", depo_signal),
                multidomain_flag=mdflag,
            )

    return None


def collect_candidates(
    record: SeqRecord,
    include_length_candidates: bool,
) -> list[RBPCandidate]:
    """
    Collect ALL RBP / tail fiber candidates from a SeqRecord.

    All matching CDS features are returned (not just the first), because
    some phages encode multiple distinct tail fiber or RBP genes.

    Parameters
    ----------
    record : SeqRecord
        A Biopython SeqRecord from a GenBank flat file.
    include_length_candidates : bool
        Whether to activate Mechanism 3.

    Returns
    -------
    list[RBPCandidate]
        All detected candidates in feature order.
    """
    from phagecore.genbank_io import resolve_organism
    organism = resolve_organism(record, getattr(record, "_source_path", None))
    accession = clean_accession(record.id)
    out: list[RBPCandidate] = []

    for feature in record.features:
        if feature.type != "CDS":
            continue
        cand = scan_cds_for_rbp(
            feature, accession, organism, include_length_candidates
        )
        if cand is not None:
            out.append(cand)

    return out


# ===========================================================================
# v2 ADDITION — InterPro reconciliation (self-contained; no phagecore import)
# ---------------------------------------------------------------------------
# Mirrors the S4 endolysin --interpro step but with host-recognition and
# depolymerase domain vocabularies. Three jobs:
#   1. Fill InterPro_Domain + InterPro_Verdict for every candidate.
#   2. Reclassify: a "tail_fiber" whose InterPro shows a pectin-lyase/beta-helix
#      or glycoside-hydrolase domain is promoted to module=both/depolymerase.
#   3. Flag false positives: a candidate whose ONLY InterPro domains are
#      structural-tail (sheath/tube/tape-measure/baseplate/portal) is marked
#      structural_tail and should be dropped from the confirmed set.
# Recovery of depolymerases annotated as bare "hypothetical protein"
# (the SAP6A class) is possible ONLY if those CDS were submitted to InterPro;
# this function reports any InterPro IDs it could not map back to a candidate so
# the analyst can add genuine recoveries manually.
# ===========================================================================

# InterPro/Pfam description substrings that CONFIRM a host-recognition RBP.
RBP_DOMAIN_TERMS: tuple[str, ...] = (
    "tail fiber", "tail fibre", "phage fiber", "phage tail fiber",
    "receptor-binding", "receptor binding", "tail spike protein",
    "t4-like", "gp37", "gp38", "short tail fiber", "fibritin",
    "phage tail protein", "tail collar", "side tail fiber",
)
# Descriptions that confirm a BIOFILM-MATRIX depolymerase (EPS / capsule / PNAG)
# — THIS is the MSc thesis variable (matrix degradation → biofilm effect).
# These are actual GLYCAN-DEGRADING ENZYME signatures.
_MATRIX_TERMS: tuple = ()   # profile.matrix_depolymerase_terms (set by configure())
# The tail-spike beta-helix FOLD is a structural scaffold that CAN carry a
# depolymerase OR an endopeptidase OR neither. Its presence alone is NOT proof of
# a matrix depolymerase — it flags a protein for structural/substrate resolution.
# v3.1.0 UNIVERSAL precision guard. These phrases denote a CHAPERONE of an enzyme
# or a pure assembly/trimerization module - never the enzyme itself. They are
# stripped from the InterPro domain text before depolymerase matching, so
# "Chaperone of endosialidase" no longer scores as a sialidase.
DEPOLYMERASE_CHAPERONE_DENYLIST: tuple[str, ...] = (
    "chaperone of endosialidase",
    "intramolecular chaperone auto-processing",
    "intramolecular chaperone",
    "endosialidase chaperone",
    "tailspike chaperone",
    "tail spike chaperone",
    "chaperone domain",
    "trimerization",
    "consensus disorder prediction",
)

TAILSPIKE_FOLD_TERMS: tuple[str, ...] = (
    "tail spike domain", "phage tail spike", "tailspike",
    # v3.0 calibration: UNIVERSAL InterPro structural signatures of polysaccharide-
    # degrading tailspikes/depolymerases (host-agnostic folds). A bare fold hit with
    # no cell-wall term -> ambiguous_depolymerase (a ColabFold/Dali target), not
    # 'uncertain'. Staph is unaffected: it carries these in _MATRIX_TERMS and
    # has_matrix is evaluated first (-> matrix_depolymerase, unchanged).
    "pectin lyase fold",        # IPR012334
    "pectin lyase-like",
    "parallel beta-helix",      # IPR006626 / IPR011050
    "right-handed beta-helix",
    "beta-helix repeat",
    "parallel beta-helix repeat",
)
# Descriptions that are CELL-WALL peptidoglycan hydrolases / tail endopeptidases /
# virion-associated lysins (VAPH). These degrade PEPTIDOGLYCAN for infection or
# lysis — they are NOT biofilm-matrix depolymerases and must NOT feed the thesis
# axis. They belong with the lysis module (S4-adjacent).
CELL_WALL_HYDROLASE_TERMS: tuple[str, ...] = (
    "endopeptidase", "peptidoglycan hydrolase", "cell-wall hydrolase",
    "cell wall hydrolase", "muramidase", "n-acetylmuramoyl", "amidase_2",
    "amidase_3", "flgj", "phosphodiester glycosidase",
    "phosphodiester alpha-n-acetylglucosaminidase", "lysozyme",
    "cysteine protease", "nlpc",
)
# Glycoside-hydrolase / glucosaminidase is AMBIGUOUS for S. aureus: PNAG is a
# beta-1,6-GlcNAc polymer, so a glucosaminidase COULD act on the biofilm matrix,
# but the same family (GH73) is the canonical cell-wall hydrolase. Flag, resolve
# case-by-case (structure + substrate), never auto-assign to the thesis axis.
_AMBIGUOUS_TERMS: tuple = ()   # profile.ambiguous_matrix_terms (set by configure())
_KNOWN_CASES: dict = {}         # profile.depolymerase_known_cases (set by configure())
_HOST_RBP_TERMS: tuple = ()     # profile.rbp_keywords (set by configure()) — v3.2 round 2:
                                # the engine previously IGNORED profile.rbp_keywords entirely,
                                # so a host whose receptor-binding module uses a different
                                # vocabulary (mycobacteriophage 'minor tail protein') produced
                                # almost no RBP candidates. Host terms are ADDITIVE to the
                                # universal EXTENDED list, so validated hosts do not regress.

# Descriptions that mark a STRUCTURAL tail/capsid/baseplate protein — NOT an RBP.
# Expanded after the first calibration run surfaced BppU / minor-tail / tail-tip /
# Dit / major-tail proteins hiding in the keyword "tail fiber" output.
STRUCTURAL_DOMAIN_DENYLIST: tuple[str, ...] = (
    "tail sheath", "tail tube", "tape measure", "baseplate", "portal",
    "major capsid", "major tail", "head-tail", "terminase", "tail terminator",
    "neck", "prohead", "scaffold", "tail completion", "tail assembly chaperone",
    "tail knob", "bppu", "baseplate upper", "minor tail protein",
    "phage tail tip", "tail tip", "distal tail", "dit-like", "gp31",
    "tail fibre assembly", "tail fiber assembly",
    # NOTE: DUF2977 is the conserved Twortvirinae TAIL-FIBRE domain (it recurs in
    # the ~302 aa tail fibers across the set) — it is NOT a structural component
    # and must NOT be denylisted, or genuine tail fibers get excluded.
)


def _interpro_verdict_from_domains(domain_text: str) -> str:
    """
    Classify a concatenated InterPro domain string for an RBP candidate.

    Verdicts (priority order):
      matrix_depolymerase : EPS/capsule/PNAG degrader  -> THESIS axis
      cell_wall_hydrolase : peptidoglycan/VAPH lysin    -> NOT thesis (lysis)
      ambiguous_gh        : glycoside hydrolase/GH73    -> resolve case-by-case
      rbp                 : receptor-binding / tail-fibre tip domain
      structural_tail     : baseplate/tail-tube/etc.    -> false positive
      no_domain / uncertain
    Structural is checked BEFORE ambiguous_gh/rbp so a baseplate protein that also
    has a coil is not mis-promoted.
    """
    d = domain_text.lower()
    if not d or "no data" in d or "no hit" in d:
        return "no_domain"
    has_struct = any(t in d for t in STRUCTURAL_DOMAIN_DENYLIST)
    # v3.1.0 UNIVERSAL precision guard: a CHAPERONE of an enzyme, or a pure
    # trimerization/assembly module, is never itself a depolymerase. Naive substring
    # matching otherwise turns "Chaperone of endosialidase" into a sialidase hit
    # (observed on E. coli T4-like long tail fibres). The chaperone PHRASES are
    # removed before enzyme matching, so a genuine depolymerase domain occurring
    # elsewhere in the same text still counts. Host-agnostic: chaperones are not
    # enzymes in any host.
    d_enz = d
    for _cp in DEPOLYMERASE_CHAPERONE_DENYLIST:
        d_enz = d_enz.replace(_cp, " ")
    has_matrix = any(t in d_enz for t in _MATRIX_TERMS)
    has_fold   = any(t in d_enz for t in TAILSPIKE_FOLD_TERMS)
    has_cw     = any(t in d for t in CELL_WALL_HYDROLASE_TERMS)
    has_gh     = any(t in d for t in _AMBIGUOUS_TERMS)
    has_rbp    = any(t in d for t in RBP_DOMAIN_TERMS)
    # A clear matrix depolymerase: a glycan-degrading enzyme signature with no
    # competing cell-wall-hydrolase signal.
    if has_matrix and not has_cw:
        return "matrix_depolymerase"
    # Ambiguous, needs structural/substrate resolution:
    #  (a) a glycan enzyme AND a cell-wall hydrolase both present, OR
    #  (b) a TAIL-SPIKE FOLD carrying a peptidase — a tail-spike endopeptidase may
    #      degrade a PROTEINACEOUS biofilm matrix (thesis-relevant) or peptidoglycan
    #      (cell-wall); the fold context distinguishes it from a free VAPH, OR
    #  (c) a bare tail-spike fold with no decisive enzyme term.
    if (has_matrix and has_cw) or (has_fold and has_cw) or \
       (has_fold and not has_cw and not has_matrix):
        return "ambiguous_depolymerase"
    # A free cell-wall hydrolase / endopeptidase with NO tail-spike fold = VAPH/lysis
    if has_cw:
        return "cell_wall_hydrolase"
    if has_gh:
        return "ambiguous_gh"
    if has_rbp:
        return "rbp"
    if has_struct:
        return "structural_tail"
    return "uncertain"


def parse_interpro_for_rbp(path: Path) -> dict[str, str]:
    """
    Parse an InterProScan TSV (15-col) OR the pasted "header = domains" summary.
    Returns {sequence_id_or_accession: concatenated_domain_text}.
    Keyed by BOTH the full first-column id and the bare accession so the
    reconciliation can match candidates by protein_id or accession.

    `path` may be a single file, a DIRECTORY (all *.tsv/*.txt inside are merged),
    or a comma-separated list of files — so a submission split into several InterPro
    batches (>100 seqs) is reconciled in one pass.
    """
    # collect the concrete files to read
    files: list[Path] = []
    if isinstance(path, Path) and path.is_dir():
        files = sorted(f for f in path.iterdir()
                       if f.suffix.lower() in (".tsv", ".txt"))
    elif "," in str(path):
        files = [Path(p.strip()) for p in str(path).split(",") if p.strip()]
    else:
        files = [path]

    out: dict[str, list[str]] = {}
    lines = []
    for f in files:
        lines += [l.rstrip("\n") for l in open(f, encoding="utf-8") if l.strip()]
    is_tsv = any(("\t" in l and " = " not in l) for l in lines[:5])

    def _push(key: str, desc: str):
        if not key or not desc or desc == "-":
            return
        out.setdefault(key, [])
        if desc not in out[key]:
            out[key].append(desc)

    if is_tsv:
        for l in lines:
            cols = l.split("\t")
            if len(cols) < 6:
                continue
            seq_id = cols[0].strip()
            # signature description (col 5) + InterPro entry name (col 12)
            for idx in (5, 12):
                if idx < len(cols):
                    _push(seq_id, cols[idx].strip())
            acc = seq_id.split("|")[0].split()[0]
            for idx in (5, 12):
                if idx < len(cols):
                    _push(acc, cols[idx].strip())
    else:
        for l in lines:
            if " = " not in l:
                continue
            head, _, dom = l.partition(" = ")
            head = head.strip(); dom = dom.strip()
            for piece in dom.split(","):
                _push(head, piece.strip())
            acc = head.split("|")[0].split()[0]
            for piece in dom.split(","):
                _push(acc, piece.strip())

    return {k: ", ".join(v) for k, v in out.items()}


def reconcile_interpro(
    candidates: list[RBPCandidate], interpro_path: Path
) -> tuple[int, int, int]:
    """
    Fill interpro_domain + interpro_verdict on each candidate, promote
    depolymerase hits to the depolymerase axis, and count outcomes.

    Returns (n_confirmed_rbp_or_depo, n_structural_false_positive, n_no_domain).
    """
    imap = parse_interpro_for_rbp(interpro_path)
    n_ok = n_fp = n_nodata = 0
    for c in candidates:
        # v3.0 D4: curated depolymerase_known_cases override (from the profile).
        # A protein_id confirmed by prior structural work (e.g. Dali) wins over the
        # keyword/InterPro verdict — this is how the B6 giant-carrier depolymerase
        # DOMAINS (invisible to keyword scan) are recorded as matrix_depolymerase.
        kc = _KNOWN_CASES.get(c.protein_id) or _KNOWN_CASES.get(c.protein_id.split(".")[0])
        if kc:
            c.interpro_verdict = kc["verdict"]
            c.depolymerase_signal = kc.get("signal", "known_case_depolymerase")
            c.module = "depolymerase" if c.subtype == "depolymerase" else "both"
            c.evidence = (c.evidence + " | " if c.evidence else "") + \
                         f"KNOWN_CASE[{kc['source']}]: {kc.get('evidence','')}"
            n_ok += 1
            continue
        # match by protein_id, accession, or accession-with-version
        domain = ""
        for key in (c.protein_id, c.accession, c.accession.split(".")[0]):
            if key in imap:
                domain = imap[key]; break
        c.interpro_domain = domain
        verdict = _interpro_verdict_from_domains(domain)
        c.interpro_verdict = verdict
        if verdict == "matrix_depolymerase":
            # the thesis-relevant class: EPS/capsule/PNAG degrader
            if c.depolymerase_signal == "none":
                c.depolymerase_signal = "interpro_matrix_depolymerase"
            c.module = "both" if c.subtype != "depolymerase" else "depolymerase"
            n_ok += 1
        elif verdict in ("cell_wall_hydrolase", "ambiguous_gh",
                         "ambiguous_depolymerase"):
            # NOT confirmed thesis axis: either lysis/VAPH, or needs case-by-case
            # structural/substrate resolution. Record, do not auto-promote.
            if c.depolymerase_signal == "none":
                c.depolymerase_signal = f"interpro_{verdict}"
            n_ok += 1
        elif verdict == "rbp":
            n_ok += 1
        elif verdict == "structural_tail":
            n_fp += 1
        elif verdict == "no_domain":
            n_nodata += 1
    return n_ok, n_fp, n_nodata


# ===========================================================================
# Deduplication
# ===========================================================================

def is_refseq(accession: str) -> bool:
    """Return True if the accession stem has a RefSeq prefix (NC_, NZ_, etc.)."""
    return any(accession.startswith(p) for p in REFSEQ_PREFIXES)


def deduplicate_files(gb_files: list[Path], length_tol: float = 0.02) -> list[Path]:
    """
    When multiple GenBank files represent the same organism, keep one PER
    DISTINCT GENOME.

    Two-step logic:
      1. Group files by organism name (GenBank ORGANISM line).
      2. Within each name group, sub-cluster by genome LENGTH. Records whose
         lengths are within `length_tol` of each other are treated as true
         duplicates (e.g., an original GenBank record and its RefSeq copy).
         Records that share the organism name but differ materially in length
         are a NAME COLLISION — genuinely different genomes deposited under the
         same name (e.g., Staphylococcus phage StAP1: KX532239 ~135.5 kb vs
         OQ025229 ~144.7 kb) — and are NOT deduplicated against each other.

    Within each true-duplicate cluster, selection priority:
      1. RefSeq accession (NC_/NZ_ prefix) — most curated
      2. Lowest accession stem alphabetically (stable, reproducible)

    Discarded files and detected name collisions are logged at WARNING level.

    Parameters
    ----------
    gb_files : list[Path]
        Sorted list of GenBank file paths.
    length_tol : float
        Fractional genome-length tolerance for treating two same-named records
        as the same genome. Default 0.02 (2%).

    Returns
    -------
    list[Path]
        Deduplicated file list; one representative per distinct genome.
    """
    from collections import defaultdict
    organism_to_entries: dict[str, list[tuple[Path, int]]] = defaultdict(list)

    for gb in gb_files:
        try:
            record = next(SeqIO.parse(str(gb), "genbank"), None)  # multi-record safe (RefSeq gbff)
            if record is None: continue
            organism = resolve_organism(record, gb).strip().lower()
            organism_to_entries[organism].append((gb, len(record.seq)))
        except Exception as exc:
            log.warning(f"  Dedup pre-scan skipped '{gb.name}': {exc}")
            organism_to_entries[gb.stem].append((gb, 0))

    kept: list[Path] = []
    for organism, entries in organism_to_entries.items():
        if len(entries) == 1:
            kept.append(entries[0][0])
            continue

        # Sub-cluster same-named records by genome length.
        clusters: list[list[tuple[Path, int]]] = []
        for path, length in sorted(entries, key=lambda e: e[1]):
            placed = False
            for cl in clusters:
                ref_len = cl[0][1]
                if ref_len > 0 and abs(length - ref_len) <= length_tol * ref_len:
                    cl.append((path, length))
                    placed = True
                    break
            if not placed:
                clusters.append([(path, length)])

        if len(clusters) > 1:
            lengths = ", ".join(f"{cl[0][1]:,}" for cl in clusters)
            log.warning(
                f"  [dedup] organism '{organism}' maps to {len(clusters)} "
                f"distinct genome lengths ({lengths} bp) — NAME COLLISION; "
                f"keeping one record PER distinct genome (not merged)"
            )

        # Keep one representative per distinct-genome cluster.
        for cl in clusters:
            files  = [p for p, _ in cl]
            refseq = [f for f in files if is_refseq(f.stem)]
            chosen = sorted(refseq, key=lambda f: f.stem)[0] if refseq \
                else sorted(files, key=lambda f: f.stem)[0]
            kept.append(chosen)
            for discarded in files:
                if discarded != chosen:
                    log.warning(
                        f"  [dedup] '{discarded.name}' duplicates "
                        f"'{chosen.name}' (organism '{organism}', "
                        f"~{cl[0][1]:,} bp) — discarded"
                    )

    return sorted(kept)


# ===========================================================================
# Batch processing and I/O
# ===========================================================================

# ---------------------------------------------------------------------------
# DepoScope reconciliation (v3.4)
# ---------------------------------------------------------------------------
# DepoScope (Boeckaerts et al., PLOS Comput Biol 2024) is an ESM-2 sequence model
# that scores every PHANOTATE-called gene for depolymerase activity and delineates
# the domain. It is ORTHOGONAL to S5: S5 detects depolymerase signal from
# InterPro/Pfam domain signatures (homology), DepoScope from learned sequence
# representations. Two methods with independent failure modes agreeing is a
# stronger claim than either alone; disagreement is a FLAG for inspection, never an
# automatic correction. Neither method overrides the other.
#
# Matching: DepoScope calls its own ORFs with PHANOTATE, so start codons can differ
# from the GenBank CDS. Primary key is the full protein-sequence md5; the fallback
# is the C-terminal 60 aa, which is stable across start-codon disagreements.

DEPOSCOPE_DEFAULT_THRESHOLD = 0.5   # the script's own internal cut-off; fix it BEFORE
                                    # looking at candidates (post-hoc tuning = selection bias)
DEPOSCOPE_CTERM_KEY_AA = 60


def _cterm_key(seq: str) -> str:
    return seq[-DEPOSCOPE_CTERM_KEY_AA:] if len(seq) >= DEPOSCOPE_CTERM_KEY_AA else seq


def parse_deposcope(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Read DepoScope output CSV(s). `path` may be one CSV or a DIRECTORY (one CSV
    per genome, which is how DepoScope must be run: single sequence per call).
    Returns (by_md5, by_cterm) score maps."""
    import csv as _csv, hashlib as _hl
    files = sorted(f for f in path.iterdir() if f.suffix.lower() == ".csv") \
        if path.is_dir() else [path]
    by_md5: dict[str, float] = {}
    by_cterm: dict[str, float] = {}
    for f in files:
        with open(f, newline="", encoding="utf-8", errors="replace") as fh:
            for row in _csv.DictReader(fh):
                prot = (row.get("protein_sequence") or "").strip().rstrip("*").upper()
                if not prot:
                    continue
                try:
                    score = float(row.get("scores_DepoScope", "nan"))
                except ValueError:
                    continue
                if score != score:      # NaN
                    continue
                m = _hl.md5(prot.encode()).hexdigest()
                by_md5[m] = max(score, by_md5.get(m, -1.0))
                ck = _cterm_key(prot)
                by_cterm[ck] = max(score, by_cterm.get(ck, -1.0))
    log.info(f"DepoScope: parsed {len(files)} file(s), {len(by_md5)} scored protein(s)")
    return by_md5, by_cterm


def reconcile_deposcope(candidates: list, deposcope_path: Path,
                        threshold: float = DEPOSCOPE_DEFAULT_THRESHOLD) -> dict:
    """Attach the DepoScope score/call to each candidate and record whether the two
    independent methods agree. Does NOT change interpro_verdict."""
    import hashlib as _hl
    by_md5, by_cterm = parse_deposcope(deposcope_path)
    depol_verdicts = ("matrix_depolymerase", "ambiguous_depolymerase")
    stats = {"scored": 0, "unmatched": 0, "agree_positive": 0, "agree_negative": 0,
             "S5_only": 0, "DepoScope_only": 0}
    for c in candidates:
        seq = (c.sequence or "").rstrip("*").upper()
        score = by_md5.get(_hl.md5(seq.encode()).hexdigest())
        if score is None:
            score = by_cterm.get(_cterm_key(seq))
        if score is None:
            c.deposcope_call = "not_evaluated"
            stats["unmatched"] += 1
            continue
        c.deposcope_score = round(score, 4)
        c.deposcope_call = "depolymerase" if score >= threshold else "non_depolymerase"
        stats["scored"] += 1
        s5_pos = c.interpro_verdict in depol_verdicts
        dp_pos = score >= threshold
        if s5_pos and dp_pos:
            c.method_agreement = "agree_positive"
        elif not s5_pos and not dp_pos:
            c.method_agreement = "agree_negative"
        elif s5_pos:
            c.method_agreement = "S5_only"
        else:
            c.method_agreement = "DepoScope_only"
        stats[c.method_agreement] += 1
    log.info("DepoScope reconciliation (threshold %.2f): %d scored, %d unmatched | "
             "agree+ %d, agree- %d, S5-only %d, DepoScope-only %d",
             threshold, stats["scored"], stats["unmatched"], stats["agree_positive"],
             stats["agree_negative"], stats["S5_only"], stats["DepoScope_only"])
    return stats


def run(
    input_dir: Path,
    fasta_out: Path,
    csv_out: Optional[Path],
    min_length: int,
    include_length_candidates: bool,
    deduplicate: bool,
    interpro_path: Optional[Path] = None,
    interpro_ready_out: Optional[Path] = None,
    colabfold_targets_out: Optional[Path] = None,
    large_rbp_min_aa: int = 0,
    rbp_triage_out: Optional[Path] = None,
    rbp_triage_carrier_aa: int = 700,
    deposcope_path: Optional[Path] = None,
    deposcope_threshold: float = DEPOSCOPE_DEFAULT_THRESHOLD,
) -> tuple[dict[str, list[RBPCandidate]], list[str], list[str]]:
    """
    Process all GenBank files and write the combined FASTA and audit CSV.

    Parameters
    ----------
    input_dir : Path
        Directory containing GenBank flat files (.gb / .gbk / .gbff).
    fasta_out : Path
        Output multi-FASTA for InterPro + MAFFT.
    csv_out : Path or None
        Audit CSV path; skipped if None.
    min_length : int
        Minimum sequence length (aa) to include in FASTA output. Does NOT
        filter the CSV — all candidates appear in the CSV regardless of length.
    include_length_candidates : bool
        Whether to activate Mechanism 3 (length heuristic).
    deduplicate : bool
        Whether to discard duplicate organism accessions before processing.

    Returns
    -------
    tuple of:
      per_genome    : dict  accession → list[RBPCandidate]
      rescued       : list  accession stems rescued from S2 "No" list
      confirmed_gap : list  accession stems still "No" after S5
    """
    gb_files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in VALID_EXTENSIONS
    )
    if not gb_files:
        log.error(f"No GenBank files found in '{input_dir}'.")
        sys.exit(1)
    log.info(f"Found {len(gb_files)} GenBank file(s) in '{input_dir}'")

    if deduplicate:
        gb_files = deduplicate_files(gb_files)
        log.info(f"After deduplication: {len(gb_files)} file(s) retained")

    per_genome:    dict[str, list[RBPCandidate]] = {}
    all_candidates: list[RBPCandidate]           = []
    rescued:        list[str]                    = []
    confirmed_gap:  list[str]                    = []

    for gb in gb_files:
        try:
            record = next(SeqIO.parse(str(gb), "genbank"), None)  # multi-record safe (RefSeq gbff)
            if record is None: continue
        except ValueError as exc:
            log.warning(f"  Skipped '{gb.name}': {exc}")
            continue
        except Exception as exc:
            log.warning(f"  Skipped '{gb.name}': unexpected error — {exc}")
            continue

        cands    = collect_candidates(record, include_length_candidates)
        acc_stem = record.id.split(".")[0]

        if cands:
            per_genome[record.id] = cands
            all_candidates.extend(cands)

            mech_summary = ", ".join(sorted({c.mechanism for c in cands}))
            rescue_tag   = ""
            if acc_stem in _S2_GAP:
                rescued.append(acc_stem)
                rescue_tag = "  [RESCUED from S2 'No']"

            log.info(
                f"  {gb.name:<22}  {len(cands):>2} candidate(s)  "
                f"[{mech_summary}]{rescue_tag}"
            )
            for c in cands:
                src_tag = "*" if c.product_source == "note/gene_qualifier" else " "
                log.info(
                    f"    {src_tag} {c.protein_id:<16}  {c.length:>5} aa  "
                    f"[{c.length_flag:<14}]  conf={c.confidence:<18}  "
                    f"'{c.product[:40]}'"
                )
        else:
            if acc_stem in _S2_GAP:
                confirmed_gap.append(acc_stem)
                log.info(
                    f"  {gb.name:<22}  no RBP detected  "
                    f"[annotation_gap_confirmed → next: tBLASTn / HHpred]"
                )
            else:
                log.info(f"  {gb.name:<22}  no RBP detected")

    # ---- v2: InterPro reconciliation (before writing outputs) ----
    if interpro_path is not None:
        if not interpro_path.exists():
            log.error(f"--interpro file not found: '{interpro_path}'")
            sys.exit(1)
        n_ok, n_fp, n_nodata = reconcile_interpro(all_candidates, interpro_path)

    # Orthogonal ESM-2 check (independent failure modes; neither overrides the other)
    if deposcope_path is not None:
        reconcile_deposcope(all_candidates, deposcope_path, deposcope_threshold)
        log.info(
            f"InterPro reconciliation: {n_ok} confirmed (rbp/depolymerase), "
            f"{n_fp} structural false-positive, {n_nodata} no-domain"
        )

    # ---- Write combined FASTA (respects --min-length) ----
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    ordered    = sorted(all_candidates, key=lambda c: (c.accession, c.subtype))
    fasta_kept = [c for c in ordered if c.length >= min_length]

    def _write_unique(path: Path, cands: list[RBPCandidate]) -> int:
        """Write FASTA guaranteeing unique headers (InterPro rejects duplicates)."""
        seen: dict[str, int] = {}
        with open(path, "w", encoding="utf-8") as fh:
            for c in cands:
                h = c.fasta_header
                if h in seen:
                    seen[h] += 1
                    h = f"{h}__dup{seen[h]}"
                else:
                    seen[h] = 1
                fh.write(f">{h}\n{c.sequence}\n")
        return len(cands)

    _write_unique(fasta_out, fasta_kept)

    excluded_by_length = len(ordered) - len(fasta_kept)
    log.info(f"\nFASTA written: '{fasta_out}'")
    log.info(f"  Sequences written     : {len(fasta_kept)}")
    log.info(f"  Excluded (< {min_length} aa)  : {excluded_by_length}")

    # ---- v2: InterPro-ready submission set (Step B2 of the calibration doc) ----
    if interpro_ready_out is not None:
        submission = [
            c for c in fasta_kept
            if (
                c.multidomain_flag == "multidomain_carrier"   # always submit giants
                or c.module in ("depolymerase", "both")        # always submit depo
                # v3.2 round 3: ALWAYS submit a LARGE carrier, whatever its confidence
                # tier. Candidates found through a host-specific vocabulary are
                # "putative" BY CONSTRUCTION (they arrive via the extended mechanism),
                # so the confirmed/high gate below excluded every one of them. On
                # M. smegmatis that meant 12 of 181 candidates reached InterPro and
                # none of the 38 carriers >=700 aa did — the exact proteins that can
                # hide a depolymerase domain. 700 aa is the same literature-grounded
                # threshold used by --rbp-triage-carrier-aa (above the 576-630 aa
                # characterised standalone depolymerase range).
                or c.length >= LARGE_CARRIER_SUBMIT_AA
                or (c.length_flag == "normal"
                    and c.confidence in ("confirmed", "high"))
            )
        ]
        # v3.0 calibration: collapse IDENTICAL sequences (md5) so the same protein
        # is never folded/submitted twice (e.g. YP_008770036 == YP_002300377). This
        # enforces the 'minimal sequences' rule for the 4 h/day ColabFold budget.
        import hashlib as _hashlib
        _seen_md5: dict[str, str] = {}
        _dedup_submission = []
        _collapsed = 0
        for c in submission:
            m = _hashlib.md5(c.sequence.encode()).hexdigest()
            if m in _seen_md5:
                _collapsed += 1
                continue
            _seen_md5[m] = c.protein_id
            _dedup_submission.append(c)
        if _collapsed:
            log.info(f"  [dedup] collapsed {_collapsed} identical-sequence duplicate(s) "
                     f"-> {len(_dedup_submission)} unique for folding/InterPro")
        submission = _dedup_submission
        if len(submission) <= INTERPRO_MAX_SEQUENCES:
            _write_unique(interpro_ready_out, submission)
            log.info(
                f"InterPro-ready FASTA: '{interpro_ready_out}'  "
                f"({len(submission)} seqs)  ✓ <= {INTERPRO_MAX_SEQUENCES}, "
                f"paste-ready in one submission."
            )
        else:
            # auto-chunk into <=100-sequence files: name_1.faa, name_2.faa, ...
            n_chunks = (len(submission) + INTERPRO_MAX_SEQUENCES - 1) // INTERPRO_MAX_SEQUENCES
            stem = interpro_ready_out.with_suffix("")
            suf  = interpro_ready_out.suffix or ".faa"
            for i in range(n_chunks):
                chunk = submission[i * INTERPRO_MAX_SEQUENCES:(i + 1) * INTERPRO_MAX_SEQUENCES]
                part  = Path(f"{stem}_{i + 1}{suf}")
                _write_unique(part, chunk)
                log.info(f"  InterPro chunk {i + 1}/{n_chunks}: '{part.name}'  ({len(chunk)} seqs)")
            log.info(
                f"InterPro-ready: {len(submission)} seqs split into {n_chunks} "
                f"files of <= {INTERPRO_MAX_SEQUENCES}. Submit each, then "
                f"concatenate the TSV results before --interpro."
            )

    if len(fasta_kept) > INTERPRO_MAX_SEQUENCES:
        log.warning(
            f"  ⚠ {len(fasta_kept)} sequences > InterPro limit "
            f"({INTERPRO_MAX_SEQUENCES}). Filter to conf=confirmed+high "
            f"and length_flag=normal, or split the FASTA into chunks."
        )
    else:
        log.info(
            f"  ✓ {len(fasta_kept)} <= {INTERPRO_MAX_SEQUENCES}: "
            f"paste-ready for InterPro in one submission."
        )

    # ---- Write audit CSV (ALL candidates, regardless of --min-length) ----
    if csv_out is not None:
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        cols = [
            "Phage", "Accession", "Protein_ID", "Product", "Length_aa",
            "Subtype", "Module", "Depolymerase_Signal", "Multidomain_Flag",
            "Confidence", "Mechanism", "Product_Source",
            "Length_Flag", "InterPro_Domain", "InterPro_Verdict",
            "Was_S2_No", "Evidence", "Sequence_MD5",
            "DepoScope_Score", "DepoScope_Call", "Method_Agreement",
        ]
        import hashlib as _hl
        with open(csv_out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for c in ordered:
                w.writerow([
                    c.organism, c.accession, c.protein_id, c.product,
                    c.length, c.subtype, c.module, c.depolymerase_signal,
                    c.multidomain_flag,
                    c.confidence, c.mechanism,
                    c.product_source, c.length_flag,
                    c.interpro_domain, c.interpro_verdict,
                    "Yes" if c.was_s2_no else "No",
                    c.evidence,
                    _hl.md5(c.sequence.encode()).hexdigest()[:12],
                    c.deposcope_score if c.deposcope_score >= 0 else "",
                    c.deposcope_call, c.method_agreement,
                ])
        log.info(f"CSV written:   '{csv_out}'  ({len(ordered)} row(s), all candidates)")
        log.info(
            "  Note: CSV includes ALL lengths; FASTA respects --min-length. "
            "Fill InterPro_Domain after the InterPro run."
        )

    # ------------------------------------------------------------------ #
    # ColabFold-ready fold targets (directly usable in AlphaFold2/ColabFold)
    # ------------------------------------------------------------------ #
    if colabfold_targets_out is not None:
        import hashlib as _hl
        import re as _re
        # only the depolymerase-axis verdicts are worth GPU time
        want = {"matrix_depolymerase", "ambiguous_depolymerase"}
        # NEVER re-fold a protein that is already structurally confirmed and
        # recorded in the profile's depolymerase_known_cases: that work is done
        # and permanent. This is what stops V.2 -> V.3 -> V.4 re-folding the same
        # proteins and protects the 4 h/day GPU budget.
        already_confirmed = set(_KNOWN_CASES) | {p.split(".")[0] for p in _KNOWN_CASES}
        targets = [c for c in all_candidates
                   if c.interpro_verdict in want
                   and c.protein_id not in already_confirmed
                   and c.protein_id.split(".")[0] not in already_confirmed]
        # RESCUE (opt-in): a GIANT receptor-binding carrier can hide a depolymerase
        # DOMAIN that neither the keyword scan nor a whole-protein InterPro verdict
        # can see — this is exactly how the S. aureus B6 depolymerases (2706-3084 aa
        # carriers) were found. For a host whose annotation yields no depolymerase
        # verdict at all, folding these carriers is the evidence-driven next step
        # (the STRUCTURE decides; no keyword is widened to manufacture a positive).
        # DepoScope_only positives: the ESM-2 model sees a depolymerase that the
        # homology/HMM route missed. These are the highest-value fold targets — a
        # candidate NOVEL depolymerase with no recognisable InterPro fold.
        dp_only = [c for c in all_candidates
                   if c.method_agreement == "DepoScope_only"
                   and c not in targets
                   and c.protein_id not in already_confirmed
                   and c.protein_id.split(".")[0] not in already_confirmed]
        if dp_only:
            log.info(f"  [DepoScope] +{len(dp_only)} target(s) the homology route missed "
                     "(candidate novel depolymerase, no recognisable InterPro fold)")
            targets += dp_only
        if large_rbp_min_aa:
            # NOTE (v3.2 round 3): the product-term re-check that used to sit here was
            # both REDUNDANT and HARMFUL. Redundant because every member of
            # all_candidates already passed the RBP scan — it IS a receptor-binding
            # candidate by construction. Harmful because the list was hardcoded
            # ("tail fiber", "tailspike", "spike", "receptor") and therefore repeated,
            # for the third time in this module, the assumption that every host uses
            # that vocabulary. On M. smegmatis all 38 carriers >=700 aa are named
            # "minor tail protein", so the rescue matched none of them and returned 0
            # targets at every threshold the operator tried (700/900/1200/1600).
            # Length is the only additional criterion the rescue needs.
            rescue = [c for c in all_candidates
                      if c.length >= large_rbp_min_aa
                      and c not in targets
                      and c.protein_id not in already_confirmed
                      and c.protein_id.split(".")[0] not in already_confirmed]
            if rescue:
                log.info(f"  [large-RBP rescue] +{len(rescue)} carrier(s) >= "
                         f"{large_rbp_min_aa} aa added as fold targets "
                         "(depolymerase domains hide inside giant carriers)")
                targets += rescue
        n_skipped = sum(1 for c in all_candidates
                        if c.protein_id in already_confirmed
                        or c.protein_id.split(".")[0] in already_confirmed)
        if n_skipped:
            log.info(f"  [known_cases] {n_skipped} already-confirmed protein(s) excluded "
                     "from folding (structural evidence already recorded in the profile)")
        # collapse identical sequences (never fold the same protein twice)
        seen: dict[str, str] = {}
        uniq = []
        for c in sorted(targets, key=lambda x: x.length):
            m = _hl.md5(c.sequence.encode()).hexdigest()
            if m in seen:
                continue
            seen[m] = c.protein_id
            uniq.append(c)
        fold = [c for c in uniq if c.length <= COLABFOLD_MAX_AA]
        giant = [c for c in uniq if c.length > COLABFOLD_MAX_AA]

        def _clean_id(pid: str) -> str:               # ColabFold-safe job name
            return _re.sub(r"[^A-Za-z0-9_]", "_", pid)

        colabfold_targets_out.parent.mkdir(parents=True, exist_ok=True)
        with open(colabfold_targets_out, "w", encoding="utf-8") as fh:
            for c in fold:
                fh.write(f">{_clean_id(c.protein_id)}\n{c.sequence}\n")
        log.info(f"ColabFold targets: '{colabfold_targets_out}' "
                 f"({len(fold)} unique depolymerase seqs <= {COLABFOLD_MAX_AA} aa, "
                 f"clean headers — feed directly to ColabFold, one per run)")
        if giant:
            gpath = colabfold_targets_out.with_suffix(".oversized.txt")
            with open(gpath, "w", encoding="utf-8") as fh:
                fh.write(f"# > {COLABFOLD_MAX_AA} aa — split by domain before folding "
                         "(fold the depolymerase-domain region only)\n")
                for c in giant:
                    fh.write(f"{_clean_id(c.protein_id)}\t{c.length} aa\t{c.product}\t"
                             f"{c.interpro_verdict}\n")
            log.info(f"  {len(giant)} oversized carrier(s) -> '{gpath.name}' "
                     "(domain-split before folding)")

    # ------------------------------------------------------------------ #
    # RBP-axis Pharokka triage (v3.3): which genomes should be re-annotated
    # BEFORE any GPU is spent on structure prediction.
    # ------------------------------------------------------------------ #
    if rbp_triage_out is not None:
        from phagecore.triage import run_rbp_triage
        run_rbp_triage(per_genome, rbp_triage_out, rbp_triage_carrier_aa)

    return per_genome, rescued, confirmed_gap


# ===========================================================================
# CLI
# ===========================================================================


# ===========================================================================
# v3.0 — host content injection (thin CLI calls this before run())
# ===========================================================================
def configure(profile) -> None:
    """Load host CONTENT from the profile into the engine (Tier-3 boundary).

    Universal DETECTION vocabulary (RBP keywords, structural denylist, cell-wall
    terms, tailspike fold, depolymerase enzyme classes) stays in this module. Only
    the host-specific THESIS-axis content is injected: which enzyme signatures are
    biofilm-matrix depolymerases for this host, which are ambiguous, and the
    per-host S2 annotation-gap list. Passing a profile whose values equal the old
    Staph constants reproduces v2.x behaviour byte-for-byte (D4 gate).
    """
    global _S2_GAP, _MATRIX_TERMS, _AMBIGUOUS_TERMS, _KNOWN_CASES, _HOST_RBP_TERMS
    _S2_GAP = frozenset(getattr(profile, "s2_annotation_gap_accessions", ()) or ())
    _MATRIX_TERMS = tuple(getattr(profile, "matrix_depolymerase_terms", ()) or ())
    _AMBIGUOUS_TERMS = tuple(getattr(profile, "ambiguous_matrix_terms", ()) or ())
    _HOST_RBP_TERMS = tuple(
        t.lower() for t in (getattr(profile, 'rbp_keywords', ()) or ())
        if t.lower() not in {k.lower() for k in PRIMARY_RBP_KEYWORDS}
        and t.lower() not in {k.lower() for k in EXTENDED_RBP_KEYWORDS})
    _KNOWN_CASES = {}
    for case in (getattr(profile, "depolymerase_known_cases", ()) or ()):
        pid = str(case.get("protein_id", "")).strip()
        if pid:
            _KNOWN_CASES[pid] = {
                "verdict": case.get("verdict", "matrix_depolymerase"),
                "evidence": case.get("evidence", ""),
                "source": case.get("source", "curated"),
                "signal": case.get("signal", "known_case_depolymerase"),
            }

