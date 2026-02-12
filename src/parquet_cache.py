"""
Automatic Parquet caching for mzML spectra.

First run loads from mzML and auto-saves a Parquet cache.
Second run detects the cache and loads instantly (~100x faster).

Pure Python module — no Streamlit imports.
"""
import hashlib
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


def compute_cache_fingerprint(
    idxml_paths: List[str],
    mzml_paths: List[str],
    target_peptides: List[str],
    reference_substring: str = "",
) -> str:
    """
    Compute a deterministic fingerprint for the current inputs.

    Uses SHA-256 of: sorted file paths + mtime + size, sorted uppercase peptides,
    lowercase reference_substring.

    Returns 12-char hex string.
    """
    h = hashlib.sha256()

    # File identity: sorted paths with mtime + size (avoids reading multi-GB files)
    file_entries = []
    for path_str in sorted(idxml_paths + mzml_paths):
        p = Path(path_str)
        if p.exists():
            stat = p.stat()
            file_entries.append(f"{p.resolve()}|{stat.st_mtime_ns}|{stat.st_size}")
        else:
            file_entries.append(f"{p.resolve()}|missing")

    for entry in file_entries:
        h.update(entry.encode("utf-8"))

    # Peptides: sorted, uppercased, stripped
    for pep in sorted(p.strip().upper() for p in target_peptides if p.strip()):
        h.update(pep.encode("utf-8"))

    # Reference substring
    h.update(reference_substring.strip().lower().encode("utf-8"))

    return h.hexdigest()[:12]


def find_existing_cache(cache_dir: str, fingerprint: str) -> Optional[str]:
    """
    Check for an existing cache file matching the fingerprint.

    Args:
        cache_dir: Directory to search for cache files
        fingerprint: 12-char hex fingerprint

    Returns:
        Path to cache file if found, None otherwise
    """
    cache_path = Path(cache_dir) / f"{fingerprint}.parquet"
    if cache_path.exists() and cache_path.stat().st_size > 0:
        logger.info(f"Cache hit: {cache_path}")
        return str(cache_path)
    return None


def get_parquet_row_counts(parquet_path: str) -> Tuple[int, int]:
    """
    Get reference and query counts from a Parquet file without loading full data.

    Reads only the is_reference column for speed.

    Returns:
        Tuple of (reference_count, query_count)
    """
    import pyarrow.parquet as pq

    try:
        table = pq.read_table(str(parquet_path), columns=["is_reference"])
        is_ref = table.column("is_reference").to_pylist()
        refs = sum(1 for v in is_ref if v)
        queries = len(is_ref) - refs
        return refs, queries
    except Exception as e:
        logger.warning(f"Could not read row counts from {parquet_path}: {e}")
        return 0, 0


def save_spectra_to_parquet(
    reference_spectra: Dict[str, list],
    query_spectra: List[Tuple[str, "LoadedSpectrum"]],
    output_path: str,
    cache_metadata: Optional[dict] = None,
) -> Optional[str]:
    """
    Save already-loaded spectra to a Parquet cache file.

    Uses the same PyArrow schema as tools/convert_mzml_to_parquet.py.

    Args:
        reference_spectra: Dict of variant_key -> list of LoadedSpectrum
        query_spectra: List of (variant_key, LoadedSpectrum) tuples
        output_path: Path for the output .parquet file
        cache_metadata: Optional dict of metadata to write to sidecar JSON

    Returns:
        Path to written file, or None on failure
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    rows = []

    # Add reference spectra
    for variant_key, specs in reference_spectra.items():
        for spec in specs:
            rows.append(_spectrum_to_row(spec, variant_key, is_reference=True))

    # Add query spectra
    for variant_key, spec in query_spectra:
        rows.append(_spectrum_to_row(spec, variant_key, is_reference=False))

    if not rows:
        logger.warning("No spectra to save to cache")
        return None

    # Build PyArrow table with same schema as convert_mzml_to_parquet.py
    schema = pa.schema([
        ("source_mzml", pa.string()),
        ("spectrum_id", pa.string()),
        ("peptidoform", pa.string()),
        ("base_sequence", pa.string()),
        ("variant_key", pa.string()),
        ("charge", pa.int32()),
        ("precursor_mz", pa.float64()),
        ("rt", pa.float64()),
        ("is_reference", pa.bool_()),
        ("mz", pa.list_(pa.float64())),
        ("intensity", pa.list_(pa.float64())),
        ("score", pa.float64()),
        ("idxml_file", pa.string()),
    ])

    arrays = {
        "source_mzml": pa.array([r["source_mzml"] for r in rows], type=pa.string()),
        "spectrum_id": pa.array([r["spectrum_id"] for r in rows], type=pa.string()),
        "peptidoform": pa.array([r["peptidoform"] for r in rows], type=pa.string()),
        "base_sequence": pa.array([r["base_sequence"] for r in rows], type=pa.string()),
        "variant_key": pa.array([r["variant_key"] for r in rows], type=pa.string()),
        "charge": pa.array([r["charge"] for r in rows], type=pa.int32()),
        "precursor_mz": pa.array([r["precursor_mz"] for r in rows], type=pa.float64()),
        "rt": pa.array([r["rt"] for r in rows], type=pa.float64()),
        "is_reference": pa.array([r["is_reference"] for r in rows], type=pa.bool_()),
        "mz": pa.array([r["mz"] for r in rows], type=pa.list_(pa.float64())),
        "intensity": pa.array([r["intensity"] for r in rows], type=pa.list_(pa.float64())),
        "score": pa.array([r["score"] for r in rows], type=pa.float64()),
        "idxml_file": pa.array([r["idxml_file"] for r in rows], type=pa.string()),
    }

    table = pa.table(arrays, schema=schema)

    # Add provenance metadata to Parquet
    provenance = {
        "created_by": "MSpecViz auto-cache",
        "created_at": datetime.now().isoformat(),
        "total_spectra": len(rows),
        "reference_spectra": sum(1 for r in rows if r["is_reference"]),
        "query_spectra": sum(1 for r in rows if not r["is_reference"]),
    }
    existing_meta = table.schema.metadata or {}
    existing_meta[b"mspecviz_provenance"] = json.dumps(provenance).encode()
    table = table.replace_schema_metadata(existing_meta)

    # Write Parquet
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, str(output))

    file_size_mb = output.stat().st_size / (1024 * 1024)
    logger.info(f"Saved cache: {len(rows)} spectra to {output} ({file_size_mb:.1f} MB)")

    # Write sidecar JSON for inspectability
    sidecar_path = output.with_suffix(".json")
    sidecar = {
        **provenance,
        "parquet_file": str(output),
        "file_size_mb": round(file_size_mb, 2),
    }
    if cache_metadata:
        sidecar.update(cache_metadata)

    try:
        with open(sidecar_path, "w") as f:
            json.dump(sidecar, f, indent=2)
        logger.info(f"Wrote sidecar metadata: {sidecar_path}")
    except Exception as e:
        logger.warning(f"Could not write sidecar JSON: {e}")

    return str(output)


def _spectrum_to_row(spec, variant_key: str, is_reference: bool) -> dict:
    """Convert a LoadedSpectrum to a row dict matching the Parquet schema."""
    score = 0.0
    idxml_file = ""
    if hasattr(spec, "match") and spec.match is not None:
        score = getattr(spec.match, "score", 0.0) or 0.0
        idxml_file = getattr(spec.match, "idxml_file", "") or ""

    return {
        "source_mzml": spec.source_file,
        "spectrum_id": spec.spectrum_id,
        "peptidoform": spec.peptidoform,
        "base_sequence": spec.base_sequence,
        "variant_key": variant_key,
        "charge": spec.precursor_charge,
        "precursor_mz": spec.precursor_mz,
        "rt": spec.rt,
        "is_reference": is_reference,
        "mz": spec.mz.tolist() if isinstance(spec.mz, np.ndarray) else list(spec.mz),
        "intensity": spec.intensity.tolist() if isinstance(spec.intensity, np.ndarray) else list(spec.intensity),
        "score": score,
        "idxml_file": idxml_file,
    }
