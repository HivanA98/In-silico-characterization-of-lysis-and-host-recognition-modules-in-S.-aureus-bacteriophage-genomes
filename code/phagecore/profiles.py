"""
phagecore.profiles
==================
Host-profile system — the seam along which Phase 2 (other bacteria) plugs in.

A HostProfile carries every piece of HOST-SPECIFIC knowledge that the engine
needs but must not hard-code:

  * taxonomy_overrides     subfamily -> family map (curated, citable)
  * size/GC/CDS bounds      plausibility envelope for QC
  * keyword sets            holin / RBP / lysis collection vocabularies
  * length priors           endolysin length windows used by ranking
  * VAPH length threshold   above which a lytic CDS is a virion-associated enzyme
  * reference accession     tBLASTn LysK reference for keyword-missed endolysins
  * active_lysis_modules    which detection modules run (Phase 2: +spanin/+lysinB)
  * curation registry       KNOWN_CASES — manually validated per-host exceptions

v3.0 ships ONE CALIBRATED built-in profile, `staphylococcus_aureus` (from the
22-phage validated set + 150/105-phage stress tests), PLUS TWO UNVALIDATED
SCAFFOLDS (`gram_negative_generic`, `mycobacterium`) whose lysis biology still
needs the spanin / lysinB modules (registered as stubs that fail loudly).
Additional host profiles are supplied as YAML files at
runtime (--profile path/to/ecoli.yaml). The YAML loader is OPTIONAL: if PyYAML
is not installed, the built-in Staph profile still works, so Phase 1 has zero
new dependencies. PyYAML is only needed to load external profiles in Phase 2.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

log = logging.getLogger("phagecore.profiles")


# ===========================================================================
# Profile container
# ===========================================================================

@dataclass
class HostProfile:
    name: str
    display_name: str
    gram: str                                   # "positive" | "negative" | "acid-fast"

    # --- taxonomy ---
    default_class: str
    taxonomy_overrides: dict                     # subfamily -> family

    # --- QC plausibility envelope ---
    min_cds_genome_size_bp: int                  # >= this with 0 CDS => FAIL
    size_bounds_bp: tuple                        # (lo, hi) warn outside
    gc_bounds_pct: tuple                         # (lo, hi) warn outside

    # --- detection vocabularies ---
    holin_keywords: tuple
    rbp_keywords: tuple
    lysis_keywords: tuple                        # broad collection net
    nonlytic_markers: tuple
    vaph_markers: tuple
    lytic_domain_terms: tuple

    # --- endolysin length priors / ranking ---
    vaph_length_threshold_aa: int
    endolysin_canonical_window_aa: tuple         # strong bonus
    endolysin_plausible_window_aa: tuple         # weak bonus

    # --- tBLASTn fallback ---
    reference_accession: str                     # auto-extracted LysK reference

    # --- special-case runtime flag sets ---
    intron_split_accessions: frozenset
    hnh_fragmented_accessions: frozenset

    # --- curation registry (validated exceptions) ---
    known_cases: dict

    # --- Phase 2 module switchboard ---
    # Modules registered in phagecore.lysis; only those listed here run.
    active_lysis_modules: tuple = ("endolysin", "holin")

    # --- v3.0: host naming of neutral evidence keys ---
    # Maps an engine evidence key (phagecore.lysis.EV_*) -> the display label
    # this host wants in the audit CSV. Engine falls back to a neutral default
    # for any key omitted. This is where "LysK-type" (a Staph name) now lives,
    # OUT of the engine. NOTE (v3.0 semantics change): `lytic_domain_terms` is
    # now the set of EVIDENCE KEYS that count as a free endolysin for this host
    # (consumed by classify_candidate), NOT free-text domain substrings as in
    # v2.x where the field was defined but never read.
    domain_labels: dict = field(default_factory=dict)

    # --- v3.0 S5 (RBP/depolymerase) host CONTENT (Tier-3 boundary) ---
    biofilm_matrix: str = ""                       # e.g. "PNAG/PIA" | "EPS/TasA/gamma-PGA"
    matrix_depolymerase_terms: tuple = ()          # enzyme signatures on the THESIS axis
    ambiguous_matrix_terms: tuple = ()             # host-specific ambiguity (Staph: GH73)
    s2_annotation_gap_accessions: tuple = ()       # per-host S2="No" list (was S2_NO_ACCESSIONS)
    depolymerase_known_cases: tuple = ()           # curated dicts {accession,verdict,evidence,source}

    # --- v3.0 S5-S7 host anchors (CONFIG; content stays per-host data) ---
    host_reference_acc: str = ""          # S7: reference host genome accession
                                          #     (e.g. NC_007795.1 aureus, NC_000964.3 subtilis)
    trna_ground_truth_acc: str = ""       # S6: primary anchor (complete isotype/anticodon annotation)
    trna_secondary_anchor_acc: str = ""   # S6: secondary anchor for concordance cross-check
    trna_canonical_isotypes: tuple = ()   # S6: DERIVED from the anchor(s), never assumed

    # --- precision denylist ---
    # Products/domains that match a broad collection keyword (e.g. "hydrolase")
    # but are NOT peptidoglycan hydrolases. Caught at classification time so they
    # are never selected as the endolysin representative (dUTPase, chaperonin,
    # ribonuclease, NTPase, tail spike, ...). Empty by default for safety.
    non_endolysin_markers: tuple = ()

    # --- dedup behaviour ---
    prefer_refseq_representative: bool = True


# ===========================================================================
# Built-in Staphylococcus aureus profile (Phase 1)
# Calibrated against the 22-phage validated set + 105-phage stress test.
# ===========================================================================

STAPH_AUREUS = HostProfile(
    name="staphylococcus_aureus",
    display_name="Staphylococcus aureus",
    gram="positive",
    default_class="Caudoviricetes",
    # Resolved from lineage first; this map is the curated fallback only.
    # Confirmed entries from the validated work are kept; entries that newly
    # appear at 105-scale are left to the lineage or flagged for ICTV check
    # rather than guessed.
    taxonomy_overrides={
        "Twortvirinae": "Herelleviridae",
        "Rakietenvirinae": "Rountreeviridae",
        "Azeredovirinae": "Unassigned",      # no family in NCBI lineage (EW, SA13)
        # Wallmarkvirinae / Bronfenbrennervirinae / Vequintavirinae intentionally
        # NOT guessed here: resolved from lineage when present, else flagged.
    },
    min_cds_genome_size_bp=10_000,           # any phage >=10 kb must have CDS
    size_bounds_bp=(15_000, 300_000),        # 16.9 kb micro .. 274 kb jumbo (observed)
    gc_bounds_pct=(28.0, 38.0),              # Staph-phage typical; FV3 43.5% -> flag
    holin_keywords=(
        "holin", "putative holin", "class i holin", "class ii holin",
        "phage holin",
    ),
    rbp_keywords=(
        "tail fiber", "tail fibre", "receptor binding protein",
        "receptor-binding protein", "tail spike", "tailspike", "tail-spike",
        "host recognition protein", "adsorption protein",
        "baseplate receptor-binding", "host specificity protein",
        # bare "rbp" handled with word-boundary matching in lysis.py
    ),
    lysis_keywords=(
        "lysin", "endolysin", "lysk", "amidase", "chap", "peptidase",
        "hydrolase", "muramidase", "glucosaminidase", "glycosidase",
        "lysozyme", "peptidoglycan", "nlpc", "p60",
    ),
    nonlytic_markers=("ybia", "nadar"),
    vaph_markers=(
        "tail", "baseplate", "structural", "virion", "lid_weld",
        "phage_lysozyme2", "tail-anchored", "tail associated", "tail-associated",
    ),
    # v3.0: EVIDENCE KEYS that count as a free endolysin (consumed by
    # classify_candidate). This set is the former hardcoded classify tuple,
    # moved out of the engine; it reproduces the v2.x free/not-free decision
    # byte-for-byte on the 150-set (hydrolase_generic / peptidase_generic are
    # deliberately EXCLUDED -> they stay 'uncertain', exactly as before).
    lytic_domain_terms=(
        "chap", "amidase", "chap_amidase", "peptidase_m15", "glucosaminidase",
        "nlpc_p60", "muramidase", "glycosidase", "endolysin_generic",
        "peptidoglycan_hydrolase",
    ),
    vaph_length_threshold_aa=700,            # Twort 1269-aa tail lysozyme is VAPH
    endolysin_canonical_window_aa=(440, 520),# full LysK
    endolysin_plausible_window_aa=(200, 600),
    reference_accession="MN336261",          # Sb1_8383 LysK
    intron_split_accessions=frozenset({"MN047438", "MF398190"}),
    hnh_fragmented_accessions=frozenset({"MN336262", "MN336263"}),
    known_cases={
        "MN047438": dict(status="intron-split", identity=None,
            note="lysK.1 N-terminal amidase moiety (~209 aa); full LysK ~495 aa "
                 "across a self-splicing intron (Kornienko 2020, Viruses, TU16)."),
        "MF398190": dict(status="intron-split", identity=None,
            note="lysK.1 N-terminal amidase moiety (~209 aa); same architecture "
                 "as MN047438."),
        "MN045228": dict(status="tblastn-recovered", identity=99.0,
            note="free LysK missed by keyword scan; recovered by tBLASTn vs "
                 "Sb1_8383 at 99% (full ~495 aa, CHAP+Amidase+SH3b)."),
        "NC_019726": dict(status="tblastn-recovered", identity=99.0,
            note="free 495-aa LysK recovered by tBLASTn vs Sb1_8383 at 99%; "
                 "keyword scan returned a 295-aa virion-associated NlpC/P60."),
        "MN336262": dict(status="hnh-disrupted", identity=None,
            note="no intact endolysin ORF (HNH endonuclease insertion); 141-aa "
                 "CHAP fragment recovered by tBLASTn (Kornienko 2023)."),
        "MN336263": dict(status="hnh-disrupted", identity=None,
            note="no intact endolysin ORF (HNH insertion); 141-aa CHAP fragment "
                 "recovered by tBLASTn."),
        "NC_007021": dict(status="vaph", identity=45.0,
            note="1269-aa phage tail lysozyme (Phage_lysozyme2), tail-anchored — "
                 "virion-associated peptidoglycan hydrolase, not a free endolysin."),
        "NC_007056": dict(status="divergent-endolysin", identity=None,
            note="576-aa NAGPA phosphodiester glycosidase, no SH3b CBD — a "
                 "divergent free endolysin (outgroup lineage)."),
    },
    active_lysis_modules=("endolysin", "holin"),  # Phase 2 adds "spanin","lysinB"
    # v3.0: reproduces the exact v2.x infer_domain() display strings so the
    # audit CSV Inferred_Domain column is byte-identical. "LysK-type" and
    # "NADAR/YbiA" — the only Staph-specific names — now live HERE, not in the
    # engine. Keys omitted here fall back to the engine's neutral defaults
    # (which already equal the v2.x strings for CHAP/Amidase/etc.).
    host_reference_acc="NC_007795.1",       # NCTC 8325 (S7 host)
    trna_ground_truth_acc="NC_023009",      # Sb-1 (S6 primary anchor, complete anticodons)
    trna_secondary_anchor_acc="",           # Staph used a single anchor
    trna_canonical_isotypes=("Ile", "Phe", "Asp"),  # Sb-1-derived, CAT-corrected: the CAT tRNA is Ile2 (ATA via lysidine), NOT Met; Trp dropped (not recovered)
    domain_labels={
        "nonlytic": "NADAR/YbiA (non-lytic)",
        "tail_lysozyme": "Phage_lysozyme2 (tail lysozyme)",
        "chap_amidase": "CHAP + Amidase (LysK-type)",
        "chap": "CHAP",
        "amidase": "Amidase",
        "peptidase_m15": "Peptidase_M15",
        "glucosaminidase": "Glucosaminidase",
        "nlpc_p60": "NlpC/P60",
        "muramidase": "Muramidase/Phage_lysozyme",
        "glycosidase": "Glycosidase (phosphodiester/NAGPA)",
        "endolysin_generic": "LysK-type (unspecified)",
        "peptidoglycan_hydrolase": "Peptidoglycan hydrolase",
        "peptidase_generic": "Peptidase (unspecified — verify)",
        "hydrolase_generic": "Hydrolase (unspecified — verify; may be non-PG)",
        "unknown": "unknown",
    },
    # Denylist calibrated from InterPro confirmation of the 59-unique set: these
    # matched the broad "hydrolase"/"peptidase" net but are not PG hydrolases.
    # Checked against the product/gene identity (not free-text notes) to avoid
    # excluding a real endolysin whose note merely mentions a neighbouring gene.
    non_endolysin_markers=(
        "dutpase", "deoxyuridine", "nucleotidohydrolase",
        "pyrophosphohydrolase", "pyrophosphatase", "dctp",      # dUTPase family
        "chaperonin", "groel", "cpn60", "chaperone",            # GroEL + chaperones
        # ("chap" keyword substring-matches "chaperone"/"chaperonin"; the denylist
        #  is checked before classification so these never become endolysins)
        "ribonuclease", "rnase", "hydroxyacylglutathione",      # RNase Z
        "p-loop", "nucleoside triphosphate hydrolase",          # P-loop NTPase
        "tail spike", "tailspike", "tail-spike",                # tail spike / depolymerase
    ),
    prefer_refseq_representative=True,
    biofilm_matrix="PNAG/PIA",
    matrix_depolymerase_terms=(
        "pectin lyase", "pectate lyase", "pectin/pectate", "parallel beta-helix",
        "right-handed beta-helix", "pectin lyase fold", "polysaccharide lyase",
        "sialidase", "neuraminidase", "bnr/asp-box", "sgnh", "gdsl",
        "carbohydrate esterase", "dextranase", "levanase", "hyaluronidase",
        "rhamnosidase", "alginate lyase",
    ),
    ambiguous_matrix_terms=(
        "glycoside hydrolase", "glycosyl hydrolase", "family 73", "glucosaminidase",
        "n-acetylglucosaminidase", "glycoside hydrolase/deacetylase",
    ),
    s2_annotation_gap_accessions=(
        "KY779849", "MN336261", "MN336262", "MN336263",
        "NC_047724", "NC_047725", "NC_047726", "NC_047727",
    ),
    # depolymerase_known_cases (D4): B6 giant-carrier depolymerase DOMAINS confirmed
    # by Dali. Keyword scan cannot see a depolymerase domain inside a 2700-3084 aa
    # tail-fibre carrier, so these are recorded as curated known_cases (source-tagged).
    depolymerase_known_cases=(
        {"accession": "MW349129.1", "protein_id": "QQO92708.1",
         "verdict": "matrix_depolymerase",
         "evidence": "Dali Z=27.4 pectate lyase fold (1ru4); domain in 2706-aa Madawaska carrier",
         "source": "Dali"},
        {"accession": "OR770614.1", "protein_id": "WPH64123.1",
         "verdict": "matrix_depolymerase",
         "evidence": "Dali Z=13.1 hyaluronate lyase fold (4d0q); domain in 2781-aa PB50 carrier",
         "source": "Dali"},
        {"accession": "NC_054982.1", "protein_id": "YP_010080028.1",
         "verdict": "matrix_depolymerase",
         "evidence": "Dali Z=30.3 SGNH esterase fold (8gkd); Sebago",
         "source": "Dali"},
        {"accession": "MW349128.1", "protein_id": "QQO92446.1",
         "verdict": "matrix_depolymerase",
         "evidence": "locked B6 lineage (v2.3.2 Dali); hyaluronate-lyase/xylanase in 3084-aa Machias carrier",
         "source": "Dali"},
        {"accession": "NC_030652.1", "protein_id": "YP_009268692.1",
         "verdict": "matrix_depolymerase",
         "evidence": "SGNH esterase (keyword + B6 lineage); phage 80; substrate pending wet-lab",
         "source": "NCBI"},
    ),

)


# ===========================================================================
# Phase-2 SCAFFOLD profiles (runnable now). active_lysis_modules deliberately
# excludes spanin/lysinB so S1/S2/S4 run on non-aureus genomes for the
# host-AGNOSTIC layers (genome stats, holin/RBP, endolysin). Full Gram-negative
# / Mycobacterium lysis characterisation awaits the spanin / lysinB modules.
# These are intentionally GENERIC starting points; tune per genus via YAML.
# ===========================================================================

_DENYLIST_COMMON = (
    "dutpase", "deoxyuridine", "nucleotidohydrolase", "pyrophosphohydrolase",
    "pyrophosphatase", "dctp", "chaperonin", "groel", "cpn60", "chaperone",
    "ribonuclease", "rnase", "hydroxyacylglutathione", "p-loop",
    "nucleoside triphosphate hydrolase", "tail spike", "tailspike", "tail-spike",
)

GRAM_NEGATIVE_GENERIC = HostProfile(
    name="gram_negative_generic",
    display_name="Gram-negative host (generic)",
    gram="negative",
    default_class="Caudoviricetes",
    taxonomy_overrides={
        # A few high-confidence current placements; everything else is resolved
        # from the lineage or flagged for ICTV verification (never guessed).
        "Tevenvirinae": "Straboviridae",
        "Studiervirinae": "Autographiviridae",
        "Vequintavirinae": "Straboviridae",
    },
    min_cds_genome_size_bp=8_000,
    size_bounds_bp=(30_000, 360_000),     # T7 ~40 kb .. jumbo phiKZ ~280 kb
    gc_bounds_pct=(33.0, 66.0),
    holin_keywords=("holin", "putative holin", "phage holin", "pinholin",
                    "antiholin"),
    rbp_keywords=(
        "tail fiber", "tail fibre", "long tail fiber", "short tail fiber",
        "receptor binding protein", "receptor-binding protein", "tail spike",
        "tailspike", "host specificity protein",
    ),
    lysis_keywords=(
        "endolysin", "lysin", "lysozyme", "amidase", "transglycosylase",
        "muramidase", "peptidase", "hydrolase", "peptidoglycan", "glycoside",
        "spanin",
    ),
    nonlytic_markers=(),
    vaph_markers=("tail", "baseplate", "structural", "virion",
                  "tail-associated"),
    # v3.0: EVIDENCE KEYS (scaffold — tune per genus; not validated).
    lytic_domain_terms=(
        "amidase", "muramidase", "glucosaminidase", "glycosidase",
        "endolysin_generic", "peptidoglycan_hydrolase", "chap", "nlpc_p60",
    ),
    vaph_length_threshold_aa=400,
    endolysin_canonical_window_aa=(140, 220),   # lambda R / T4 e / T7
    endolysin_plausible_window_aa=(120, 400),
    reference_accession="",                       # supply per host via --reference
    intron_split_accessions=frozenset(),
    hnh_fragmented_accessions=frozenset(),
    known_cases={},
    active_lysis_modules=("endolysin", "holin"),  # spanin NOT yet implemented
    non_endolysin_markers=_DENYLIST_COMMON,
    prefer_refseq_representative=True,
)

MYCOBACTERIUM = HostProfile(
    name="mycobacterium",
    display_name="Mycobacterium host (acid-fast)",
    gram="acid-fast",
    default_class="Caudoviricetes",
    taxonomy_overrides={},
    min_cds_genome_size_bp=8_000,
    size_bounds_bp=(40_000, 165_000),     # D29 ~49 kb, jumbo myco up to ~160 kb
    gc_bounds_pct=(56.0, 72.0),           # mycobacteriophages are high-GC
    holin_keywords=("holin", "putative holin", "phage holin"),
    rbp_keywords=("tail fiber", "tail fibre", "receptor binding protein",
                  "minor tail protein", "tail tip"),
    lysis_keywords=(
        "endolysin", "lysin", "lysin a", "lysa", "lysin b", "lysb", "amidase",
        "peptidase", "muramidase", "lysozyme", "hydrolase", "peptidoglycan",
        "esterase", "cutinase",
    ),
    nonlytic_markers=(),
    vaph_markers=("tail", "structural", "virion"),
    # v3.0: EVIDENCE KEYS (scaffold — tune per genus; not validated).
    lytic_domain_terms=(
        "amidase", "muramidase", "peptidoglycan_hydrolase", "endolysin_generic",
        "glucosaminidase", "chap",
    ),
    vaph_length_threshold_aa=600,
    endolysin_canonical_window_aa=(250, 480),   # LysA is larger/modular
    endolysin_plausible_window_aa=(150, 600),
    reference_accession="",
    intron_split_accessions=frozenset(),
    hnh_fragmented_accessions=frozenset(),
    known_cases={},
    active_lysis_modules=("endolysin", "holin"),  # lysinB NOT yet implemented
    non_endolysin_markers=_DENYLIST_COMMON,
    prefer_refseq_representative=True,
)


_BUILTINS = {p.name: p for p in (STAPH_AUREUS, GRAM_NEGATIVE_GENERIC,
                                 MYCOBACTERIUM)}


# ===========================================================================
# Loader
# ===========================================================================

def load_profile(spec: Optional[str]) -> HostProfile:
    """
    Resolve a profile from a name (built-in) or a path to a YAML file.

    spec is None                             -> ValueError (no silent default).
    spec == a built-in name                  -> that built-in profile.
    spec ending in .yaml/.yml                -> external profile.
    """
    if spec is None:
        raise ValueError(
            "No host profile specified. v3.0 removes the silent Staphylococcus "
            "aureus default to prevent wrong-host runs. Pass --profile <name|path>. "
            f"Built-ins: {', '.join(sorted(_BUILTINS))} (only staphylococcus_aureus "
            "is calibrated). Or give a path to a .yaml profile.")
    if spec in _BUILTINS:
        return _BUILTINS[spec]

    path = Path(spec)
    if path.suffix.lower() in (".yaml", ".yml"):
        return _load_yaml_profile(path)

    raise ValueError(
        f"Unknown profile '{spec}'. Use a built-in name "
        f"({', '.join(sorted(_BUILTINS))}) or a path to a .yaml profile."
    )


def _load_yaml_profile(path: Path) -> HostProfile:
    """Load an external profile. Requires PyYAML (Phase 2 only)."""
    try:
        import yaml
    except ImportError as exc:                          # pragma: no cover
        raise ImportError(
            "Loading an external profile needs PyYAML: pip install pyyaml. "
            "The built-in 'staphylococcus_aureus' profile needs no extra deps."
        ) from exc
    data = yaml.safe_load(path.read_text(encoding="utf-8"))

    def tup(key, default=()):
        return tuple(data.get(key, default) or default)

    return HostProfile(
        name=data["name"],
        display_name=data.get("display_name", data["name"]),
        gram=data.get("gram", "unknown"),
        default_class=data.get("default_class", "Caudoviricetes"),
        taxonomy_overrides=dict(data.get("taxonomy_overrides", {})),
        min_cds_genome_size_bp=int(data.get("min_cds_genome_size_bp", 10_000)),
        size_bounds_bp=tuple(data.get("size_bounds_bp", (10_000, 500_000))),
        gc_bounds_pct=tuple(data.get("gc_bounds_pct", (20.0, 70.0))),
        holin_keywords=tup("holin_keywords"),
        rbp_keywords=tup("rbp_keywords"),
        lysis_keywords=tup("lysis_keywords"),
        nonlytic_markers=tup("nonlytic_markers"),
        vaph_markers=tup("vaph_markers"),
        lytic_domain_terms=tup("lytic_domain_terms"),
        vaph_length_threshold_aa=int(data.get("vaph_length_threshold_aa", 700)),
        endolysin_canonical_window_aa=tuple(
            data.get("endolysin_canonical_window_aa", (440, 520))),
        endolysin_plausible_window_aa=tuple(
            data.get("endolysin_plausible_window_aa", (150, 600))),
        reference_accession=data.get("reference_accession", ""),
        intron_split_accessions=frozenset(data.get("intron_split_accessions", [])),
        hnh_fragmented_accessions=frozenset(data.get("hnh_fragmented_accessions", [])),
        known_cases=dict(data.get("known_cases", {})),
        active_lysis_modules=tuple(data.get("active_lysis_modules",
                                            ("endolysin", "holin"))),
        domain_labels=dict(data.get("domain_labels", {})),
        host_reference_acc=data.get("host_reference_acc", ""),
        trna_ground_truth_acc=data.get("trna_ground_truth_acc", ""),
        trna_secondary_anchor_acc=data.get("trna_secondary_anchor_acc", ""),
        trna_canonical_isotypes=tuple(data.get("trna_canonical_isotypes", ())),
        biofilm_matrix=data.get("biofilm_matrix", ""),
        matrix_depolymerase_terms=tuple(data.get("matrix_depolymerase_terms", ())),
        ambiguous_matrix_terms=tuple(data.get("ambiguous_matrix_terms", ())),
        s2_annotation_gap_accessions=tuple(data.get("s2_annotation_gap_accessions", ())),
        depolymerase_known_cases=tuple(data.get("depolymerase_known_cases", ())),
        non_endolysin_markers=tup("non_endolysin_markers"),
        prefer_refseq_representative=bool(
            data.get("prefer_refseq_representative", True)),
    )


def available_builtin_profiles() -> list[str]:
    return sorted(_BUILTINS)
