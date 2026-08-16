"""
phagecore.genbank_io
====================
File discovery, streaming GenBank parsing, content checksums, per-file error
isolation, and the run manifest. Built to ingest MASSIVE batches of GenBank
files without loading them all into memory and without letting a single
malformed file abort the run.

Why streaming + per-file isolation
----------------------------------
A name-only batch loop that calls SeqIO.read() inside one try/except is fine for
22 files but fragile at 10^3-10^4 files: one corrupt record kills the whole run,
memory grows if records are retained, and there is no record of what was
processed. This module yields one ParsedGenome at a time, isolates failures into
the manifest, and never holds more than the current record in memory.
"""

from __future__ import annotations

import hashlib
import logging
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import VALID_EXTENSIONS, __version__, __phase__

log = logging.getLogger("phagecore.io")


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass

class ParsedGenome:
    """One successfully parsed GenBank record plus engine-side metadata."""
    record: SeqRecord
    source_file: Path
    accession: str            # full, e.g. NC_023573.1
    accession_base: str       # version-stripped, e.g. NC_023573
    organism: str
    seq_md5: str              # md5 of upper-case nucleotide sequence (dedup key)
    feature_census: dict      # {feature_type: count} — diagnostic for QC
    is_refseq: bool           # NC_/NZ_/NG_... — provenance for dedup


@dataclass
class FileOutcome:
    """Per-file processing record for the manifest."""
    source_file: str
    status: str               # "ok" | "parse_error" | "empty"
    accession: str = ""
    organism: str = ""
    seq_len: int = 0
    n_cds: int = 0
    seq_md5: str = ""
    message: str = ""


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def discover_genbank_files(input_dir: Path, recursive: bool = False) -> list[Path]:
    """
    Return a deterministically sorted list of GenBank files.

    recursive=True walks sub-directories (useful when a download tool shards
    thousands of files into nested folders). Sorting guarantees reproducible
    ordering regardless of filesystem enumeration order.
    """
    it = input_dir.rglob("*") if recursive else input_dir.iterdir()
    files = [p for p in it if p.is_file() and p.suffix.lower() in VALID_EXTENSIONS]
    return sorted(files)


# ---------------------------------------------------------------------------
# Checksums & provenance
# ---------------------------------------------------------------------------

def sequence_md5(seq: str) -> str:
    """md5 of the upper-cased nucleotide sequence — the genome-level dedup key."""
    return hashlib.md5(seq.upper().encode("ascii", "ignore")).hexdigest()


_REFSEQ_PREFIXES = ("NC_", "NZ_", "NG_", "AC_", "NW_", "NT_")


def is_refseq_accession(accession: str) -> bool:
    """True if the accession is an NCBI RefSeq copy (vs a primary INSDC deposit)."""
    return accession.upper().startswith(_REFSEQ_PREFIXES)


# ---------------------------------------------------------------------------
# Feature census (diagnostic that distinguishes the two CDS=0 root causes)
# ---------------------------------------------------------------------------

def census_features(record: SeqRecord) -> dict:
    """
    Count features by type. This separates the two failure modes that both
    present as CDS_Count=0:
      - annotation truly absent      → {'source': 1} only
      - CDS stored under a different
        structure / not parsed        → e.g. {'gene': 60, ...} with no 'CDS'
    The manifest carries this so the cause is diagnosable without reopening files.
    """
    census: dict[str, int] = {}
    for f in record.features:
        census[f.type] = census.get(f.type, 0) + 1
    return census


# ---------------------------------------------------------------------------
# Streaming parser with per-file isolation
# ---------------------------------------------------------------------------

def parse_genomes(files: list[Path]) -> Iterator[tuple[Optional[ParsedGenome], FileOutcome]]:
    """
    Yield (ParsedGenome | None, FileOutcome) for each file.

    Uses SeqIO.parse (not read) so multi-record files do not raise; the first
    complete record per file is taken (one genome per file is the convention).
    Any exception is captured into FileOutcome with status != "ok" and the run
    continues. Memory stays flat: only the current record is held.
    """
    for path in files:
        try:
            record = next(SeqIO.parse(str(path), "genbank"), None)
            if record is None:
                yield None, FileOutcome(str(path), "empty",
                                        message="no GenBank record found")
                continue
            seq = str(record.seq)
            census = census_features(record)
            pg = ParsedGenome(
                record=record,
                source_file=path,
                accession=record.id,
                accession_base=record.id.split(".")[0],
                organism=resolve_organism(record, path),
                seq_md5=sequence_md5(seq),
                feature_census=census,
                is_refseq=is_refseq_accession(record.id),
            )
            outcome = FileOutcome(
                str(path), "ok", record.id, pg.organism,
                len(seq), census.get("CDS", 0), pg.seq_md5,
            )
            yield pg, outcome
        except Exception as exc:                       # noqa: BLE001 (isolation is intentional)
            log.warning("  Skipped '%s': %s", path.name, exc)
            yield None, FileOutcome(str(path), "parse_error", message=str(exc))


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def tool_versions() -> dict:
    """Capture the environment so any output table is reproducible later."""
    try:
        import Bio
        biover = Bio.__version__
    except Exception:                                  # noqa: BLE001
        biover = "unknown"
    try:
        import pandas
        pdver = pandas.__version__
    except Exception:                                  # noqa: BLE001
        pdver = "unknown"
    return {
        "phagecore": f"{__version__} ({__phase__})",
        "python": sys.version.split()[0],
        "biopython": biover,
        "pandas": pdver,
        "platform": platform.platform(),
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def write_manifest(path: Path, outcomes: list[FileOutcome], profile_name: str,
                   script: str, extra: Optional[dict] = None) -> None:
    """Write a per-run manifest CSV (provenance + reproducibility record)."""
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    versions = tool_versions()
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["# phage characterization run manifest"])
        w.writerow(["# script", script])
        w.writerow(["# profile", profile_name])
        for k, v in versions.items():
            w.writerow([f"# {k}", v])
        if extra:
            for k, v in extra.items():
                w.writerow([f"# {k}", v])
        w.writerow([])
        w.writerow(["source_file", "status", "accession", "organism",
                    "seq_len", "n_cds", "seq_md5", "message"])
        for o in outcomes:
            w.writerow([o.source_file, o.status, o.accession, o.organism,
                        o.seq_len, o.n_cds, o.seq_md5, o.message])
    log.info("Manifest written: '%s' (%d files)", path, len(outcomes))

def resolve_organism(record, source_path=None) -> str:
    """Organism name with a documented fallback chain.

    Pharokka-re-annotated GenBank files carry `ORGANISM  .` — Pharokka does not
    preserve the source organism. Reading record.annotations["organism"] then
    faithfully yields "." and the Phage/Organism column looks empty. That is a DATA
    gap, not a parsing bug, so the chain is explicit and ordered:
        1. annotations["organism"]      (the real value when present)
        2. record.description           (Pharokka sometimes writes a name here)
        3. the source FILENAME stem     (the operator's own naming, e.g.
                                         CP062445_Aureus_phage_ECel-2020o)
        4. record.name / record.id      (last resort)
    Values that are placeholders (".", "unknown", "") are treated as absent.
    """
    _PLACEHOLDER = {"", ".", "unknown", "unclassified", "n/a", "na"}
    org = str(record.annotations.get("organism", "") or "").strip()
    if org.lower() not in _PLACEHOLDER:
        return org
    desc = str(getattr(record, "description", "") or "").strip()
    if desc.lower() not in _PLACEHOLDER:
        return desc
    if source_path is not None:
        from pathlib import Path as _P
        stem = _P(str(source_path)).stem.strip()
        if stem.lower() not in _PLACEHOLDER:
            return stem
    return str(getattr(record, "name", "") or getattr(record, "id", "") or "unknown")
