"""
phagecore — shared engine for NCBI bacteriophage genome characterization
========================================================================

This package is the INVARIANT CORE used by the three data-collection scripts
S1 (genome statistics), S2 (holin/RBP annotation) and S4 (endolysin extraction).
Everything in here is host-agnostic: GenBank parsing, QC, deduplication,
taxonomy resolution, manifest/provenance, and the detection harness.

Host-specific knowledge (length priors, keyword sets, taxonomy overrides,
curation registry, which lysis modules are active) lives OUTSIDE the engine,
in a HostProfile (see phagecore.profiles). Phase 1 ships one profile,
`staphylococcus_aureus`. Phase 2 adds profiles for other hosts WITHOUT editing
the engine — that is the whole point of this split.

Design contract
----------------
1. One bad input file never aborts a run (per-file error isolation).
2. Biologically impossible output (e.g. a 40 kb genome with 0 CDS) is FLAGGED,
   never emitted silently.
3. Identical sequences are FLAGGED with provenance, never silently deleted.
4. Every run writes a manifest (versions, checksums, per-file status) so any
   table is reproducible and auditable years later.

Versioning
----------
__version__ tracks the engine. Bump the minor for new capabilities that keep
backward-compatible output columns; bump the major for column/format changes.
"""

from __future__ import annotations

__version__ = "3.0.0"
__phase__ = "phase-2-multihost"

# Accepted GenBank flat-file extensions (case-insensitive).
VALID_EXTENSIONS = frozenset({".gb", ".gbk", ".gbff", ".genbank"})

# Canonical classification vocabulary (LOWERCASE — single source of truth).
# The 22-phage validated set uses exactly these labels; the 105-phage run
# accidentally introduced an upper-case "VAPH" duplicate. The engine normalises
# every label through this set so downstream group-by never fragments.
CLASS_FREE = "free-endolysin"
CLASS_VAPH = "vaph"
CLASS_NONLYTIC = "non-lytic"
CLASS_INTRON = "intron-split"
CLASS_HNH = "hnh-disrupted"
CLASS_TBLASTN = "tblastn-recovered"
CLASS_DIVERGENT = "divergent-endolysin"
CLASS_UNCERTAIN = "uncertain"
CLASS_NONENDO = "non-endolysin"      # matched a broad keyword but is not a PG hydrolase

CANONICAL_CLASSES = frozenset({
    CLASS_FREE, CLASS_VAPH, CLASS_NONLYTIC, CLASS_INTRON, CLASS_HNH,
    CLASS_TBLASTN, CLASS_DIVERGENT, CLASS_UNCERTAIN, CLASS_NONENDO,
})


def normalise_class(label: str) -> str:
    """Map any-case classification label to the canonical lowercase token."""
    return (label or "").strip().lower()
