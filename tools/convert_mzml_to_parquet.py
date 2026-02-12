#!/usr/bin/env python3
"""
Standalone conversion script: mzML + idXML → Parquet.

Extracts only matching spectra from mzML files based on idXML search results,
writing a compact Parquet file with all data needed for MSpecViz analysis.

Usage (CLI):
    python tools/convert_mzml_to_parquet.py \
        --idxml /path/to/file1.idXML /path/to/file2.idXML \
        --mzml-dir /path/to/mzml_files/ \
        --peptides "PEPTIDEK,ANOTHERSEQ" \
        --output /path/to/output.parquet \
        --reference-substring "Synpep"

Usage (Python API):
    from tools.convert_mzml_to_parquet import convert_mzml_to_parquet
    output_path = convert_mzml_to_parquet(
        idxml_paths=[...],
        mzml_dir="/path/to/mzmls",
        target_peptides=["PEPTIDEK"],
        output_path="output.parquet",
    )

Dependencies: pyopenms, pyarrow, numpy (+ src/peptide_search.py, src/spectrum_loader.py)
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np

# Add project root to path so we can import src modules
_project_root = str(Path(__file__).resolve().parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import pyopenms as oms
import pyarrow as pa
import pyarrow.parquet as pq

from src.peptide_search import search_peptides_in_idxmls
from src.spectrum_loader import (
    build_mzml_to_matches_map,
    open_indexed_mzml_fast,
    find_spectrum_index,
    clean_spectrum_id,
)

logger = logging.getLogger(__name__)


def convert_mzml_to_parquet(
    idxml_paths: List[str],
    mzml_dir: str,
    target_peptides: List[str],
    output_path: str,
    reference_substring: str = "",
    progress_callback: Optional[Callable[[str], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None,
) -> Optional[str]:
    """
    Convert idXML+mzML data into a compact Parquet file.

    Searches idXML files for target peptides, resolves matching spectra
    in mzML files, extracts only those spectra, and writes them to Parquet.

    Args:
        idxml_paths: List of idXML file paths to search
        mzml_dir: Directory containing mzML files referenced by idXML
        target_peptides: List of peptide base sequences to find
        output_path: Path for output Parquet file
        reference_substring: Substring to identify reference mzML files
        progress_callback: Optional callback(message) for progress updates
        error_callback: Optional callback(message) for errors

    Returns:
        Path to written Parquet file, or None on failure
    """
    def _progress(msg: str):
        logger.info(msg)
        if progress_callback:
            progress_callback(msg)

    def _error(msg: str):
        logger.error(msg)
        if error_callback:
            error_callback(msg)

    # =========================================================================
    # Step 1: Search idXML files for target peptides
    # =========================================================================
    _progress(f"Searching {len(idxml_paths)} idXML files for {len(target_peptides)} peptides...")

    variant_groups = search_peptides_in_idxmls(
        idxml_paths=idxml_paths,
        target_peptides=target_peptides,
        error_callback=_error,
    )

    if not variant_groups:
        _error("No matching peptides found in idXML files")
        return None

    total_psms = sum(len(vg.matches) for vg in variant_groups.values())
    _progress(f"Found {total_psms} PSMs across {len(variant_groups)} variants")

    # =========================================================================
    # Step 2: Discover available mzML files
    # =========================================================================
    mzml_dir_path = Path(mzml_dir)
    mzml_paths = sorted(mzml_dir_path.glob("*.mzML")) + sorted(mzml_dir_path.glob("*.mzml"))

    if not mzml_paths:
        _error(f"No mzML files found in {mzml_dir}")
        return None

    mzml_path_map = {p.name.lower(): p for p in mzml_paths}
    available_names = set(mzml_path_map.keys())

    _progress(f"Found {len(mzml_paths)} mzML files in {mzml_dir}")

    # =========================================================================
    # Step 3: Build mzML → matches mapping
    # =========================================================================
    mzml_to_matches = build_mzml_to_matches_map(variant_groups, available_names)

    if not mzml_to_matches:
        _error("No mzML files could be mapped to matches. Check id_merge_index resolution.")
        return None

    mapped_spectra = sum(len(m) for m in mzml_to_matches.values())
    _progress(f"Mapped {mapped_spectra} spectra across {len(mzml_to_matches)} mzML files")

    # =========================================================================
    # Step 4: Extract spectra from mzML files
    # =========================================================================
    rows = []
    validation_warnings = []
    files_processed = 0

    for mzml_lower, matches_for_file in mzml_to_matches.items():
        mzml_path = mzml_path_map[mzml_lower]
        is_reference = bool(
            reference_substring
            and reference_substring.lower() in mzml_path.name.lower()
        )

        files_processed += 1
        _progress(f"[{files_processed}/{len(mzml_to_matches)}] Extracting from {mzml_path.name} "
                  f"({len(matches_for_file)} spectra, {'reference' if is_reference else 'query'})...")

        try:
            needed_ids = [match.spectrum_id for match, _ in matches_for_file]
            indexed_exp, id_to_index = open_indexed_mzml_fast(str(mzml_path), needed_ids=needed_ids)
            n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

            for match, variant_key in matches_for_file:
                try:
                    idx = find_spectrum_index(match.spectrum_id, id_to_index, n_spectra)

                    if idx is None:
                        warn = f"Spectrum {match.spectrum_id} not found in {mzml_path.name}"
                        validation_warnings.append(warn)
                        _error(warn)
                        continue

                    spectrum = indexed_exp.getSpectrum(idx)

                    # Validate: MS level should be 2
                    ms_level = spectrum.getMSLevel()
                    if ms_level != 2:
                        warn = f"Spectrum {match.spectrum_id} in {mzml_path.name} is MS{ms_level}, expected MS2"
                        validation_warnings.append(warn)
                        logger.warning(warn)

                    # Extract peak data
                    mz_array = np.array([p.getMZ() for p in spectrum], dtype=np.float64)
                    intensity_array = np.array([p.getIntensity() for p in spectrum], dtype=np.float64)

                    # Get precursor info from spectrum
                    precursor_mz = 0.0
                    precursor_charge = match.charge
                    if spectrum.getPrecursors():
                        prec = spectrum.getPrecursors()[0]
                        precursor_mz = prec.getMZ()
                        if prec.getCharge() > 0:
                            precursor_charge = prec.getCharge()

                    # Provenance: cross-check precursor m/z if available
                    if precursor_mz > 0 and match.charge > 0:
                        # Very rough sanity check (> 1 Da difference = likely wrong spectrum)
                        expected_mz = precursor_mz  # from mzML
                        # We don't have match.precursor_mz but we can log the value for traceability

                    rt = spectrum.getRT() if spectrum.getRT() else 0.0

                    rows.append({
                        "source_mzml": mzml_path.name,
                        "spectrum_id": match.spectrum_id,
                        "peptidoform": match.peptidoform,
                        "base_sequence": match.base_sequence,
                        "variant_key": variant_key,
                        "charge": precursor_charge,
                        "precursor_mz": precursor_mz,
                        "rt": rt,
                        "is_reference": is_reference,
                        "mz": mz_array.tolist(),
                        "intensity": intensity_array.tolist(),
                        "score": match.score,
                        "idxml_file": match.idxml_file,
                    })

                except Exception as e:
                    warn = f"Error extracting {match.spectrum_id} from {mzml_path.name}: {e}"
                    validation_warnings.append(warn)
                    _error(warn)

        except Exception as e:
            _error(f"Error opening {mzml_path.name}: {e}")
            continue

    if not rows:
        _error("No spectra could be extracted")
        return None

    _progress(f"Extracted {len(rows)} spectra from {files_processed} files")

    # =========================================================================
    # Step 5: Write Parquet file
    # =========================================================================
    _progress(f"Writing Parquet to {output_path}...")

    # Build PyArrow arrays for type safety
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

    # Add provenance metadata
    provenance = {
        "created_by": "MSpecViz convert_mzml_to_parquet",
        "created_at": datetime.now().isoformat(),
        "idxml_files": [str(p) for p in idxml_paths],
        "mzml_dir": str(mzml_dir),
        "target_peptides": target_peptides,
        "reference_substring": reference_substring,
        "total_spectra": len(rows),
        "total_variants": len(set(r["variant_key"] for r in rows)),
        "reference_spectra": sum(1 for r in rows if r["is_reference"]),
        "query_spectra": sum(1 for r in rows if not r["is_reference"]),
        "validation_warnings": len(validation_warnings),
    }

    existing_meta = table.schema.metadata or {}
    existing_meta[b"mspecviz_provenance"] = json.dumps(provenance).encode()
    table = table.replace_schema_metadata(existing_meta)

    # Write
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, str(output_path_obj))

    file_size_mb = output_path_obj.stat().st_size / (1024 * 1024)
    _progress(f"Done! Wrote {len(rows)} spectra to {output_path} ({file_size_mb:.1f} MB)")

    if validation_warnings:
        _progress(f"Warnings: {len(validation_warnings)}")
        for w in validation_warnings[:5]:
            logger.warning(f"  {w}")
        if len(validation_warnings) > 5:
            logger.warning(f"  ... and {len(validation_warnings) - 5} more")

    return str(output_path_obj)


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Convert idXML+mzML to Parquet for fast spectrum loading",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert with specific peptides
  python tools/convert_mzml_to_parquet.py \\
    --idxml data/*.idXML \\
    --mzml-dir data/mzml/ \\
    --peptides "PEPTIDEK,ANOTHERSEQ" \\
    --output output.parquet

  # Convert with reference identification
  python tools/convert_mzml_to_parquet.py \\
    --idxml data/*.idXML \\
    --mzml-dir data/mzml/ \\
    --peptides "PEPTIDEK" \\
    --output output.parquet \\
    --reference-substring "Synpep"
        """,
    )

    parser.add_argument(
        "--idxml", nargs="+", required=True,
        help="Path(s) to idXML file(s)",
    )
    parser.add_argument(
        "--mzml-dir", required=True,
        help="Directory containing mzML files",
    )
    parser.add_argument(
        "--peptides", required=True,
        help="Comma-separated list of target peptide sequences",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output Parquet file path",
    )
    parser.add_argument(
        "--reference-substring", default="",
        help="Substring to identify reference mzML files (e.g., 'Synpep')",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    # Parse peptides
    peptides = [p.strip() for p in args.peptides.split(",") if p.strip()]
    if not peptides:
        parser.error("No valid peptides provided")

    # Validate inputs
    for p in args.idxml:
        if not Path(p).exists():
            parser.error(f"idXML file not found: {p}")

    if not Path(args.mzml_dir).is_dir():
        parser.error(f"mzML directory not found: {args.mzml_dir}")

    # Run conversion
    result = convert_mzml_to_parquet(
        idxml_paths=args.idxml,
        mzml_dir=args.mzml_dir,
        target_peptides=peptides,
        output_path=args.output,
        reference_substring=args.reference_substring,
        progress_callback=lambda msg: print(f"  {msg}"),
        error_callback=lambda msg: print(f"  ERROR: {msg}", file=sys.stderr),
    )

    if result:
        print(f"\nSuccess: {result}")
        sys.exit(0)
    else:
        print("\nConversion failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
