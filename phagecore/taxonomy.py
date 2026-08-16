"""
phagecore.taxonomy
=================
ICTV-aware taxonomy resolution from the GenBank LINEAGE field.

Improvement over the original
-----------------------------
The original resolved family ONLY from a 3-entry hard-coded dict whenever the
'-viridae' token was absent, and wrote "Unassigned" silently. At 105 genomes
that produced 50/105 "Unassigned" families, many of them simply because new
subfamilies (Wallmarkvirinae, Bronfenbrennervirinae, ...) were not in the dict.

This version resolves in priority order:
  1. the lineage tokens themselves (suffix-based; most authoritative)
  2. the profile's taxonomy_overrides map (small, documented, citable)
  3. fallback to "Unassigned" — but ALWAYS with a taxonomy_flag telling the user
     to verify against the current ICTV VMR. Silent settledness is the bug; an
     explicit "verify this" is the fix.

The override map lives in the HostProfile, NOT here: E. coli/Salmonella phages
have entirely different families, so Phase 2 supplies its own map without
touching this engine code.
"""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord


def resolve_taxonomy(record: SeqRecord, profile) -> tuple[str, str, str, str]:
    """
    Return (klass, family, subfamily, taxonomy_flag).

    taxonomy_flag is "" when family came from the lineage or a curated override,
    or "family_unresolved_verify_ICTV" when it fell back to Unassigned.
    """
    lineage = record.annotations.get("taxonomy", [])

    klass = family = subfamily = None
    for taxon in lineage:
        t = taxon.strip()
        if t.endswith("viricetes"):
            klass = t
        elif t.endswith("viridae"):
            family = t
        elif t.endswith("virinae"):
            subfamily = t

    if klass is None:
        klass = profile.default_class

    if subfamily is None:
        subfamily = "Unassigned"

    flag = ""
    if family is None:
        mapped = profile.taxonomy_overrides.get(subfamily)
        if mapped and mapped != "Unassigned":
            family = mapped
        else:
            family = "Unassigned"
            flag = "family_unresolved_verify_ICTV"

    return klass, family, subfamily, flag


def infer_ncbi_status(record: SeqRecord) -> str:
    """Completeness from KEYWORDS / DEFINITION (unchanged convention)."""
    keywords = record.annotations.get("keywords", [])
    description = record.description.lower()
    if any("complete" in kw.lower() for kw in keywords) or "complete" in description:
        return "Complete Genome"
    return "Draft/Partial"
