"""High-level peptide-search orchestration.

Coordinates parallel/sequential execution of the per-file workers in
:mod:`mspecviz.search.workers` and aggregates results into ``VariantGroup``
dictionaries that the rest of the pipeline consumes.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import polars as pl

from mspecviz.common.parallel import run_tasks
from mspecviz.common.path_utils import disambiguate_paths
from mspecviz.search.models import VariantGroup, format_inline_sequence
from mspecviz.search.workers import (
    _scan_single_idxml_peptides,
    _search_single_idxml,
)

logger = logging.getLogger(__name__)


def search_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[callable] = None,
    error_callback: Optional[callable] = None,
    workers: int = 1,
) -> Dict[str, VariantGroup]:
    """
    Search for target peptides across multiple idXML files using PyOpenMS.

    Extracts all metadata needed for mzML resolution including id_merge_index.

    Args:
        idxml_paths: List of idXML file paths to search
        target_peptides: List of peptide sequences to find (base sequences)
        mzml_paths: Optional list of selected mzML file paths (for disambiguation)
        progress_callback: Optional callback(current, total, message) for progress
        error_callback: Optional callback(error_message) for immediate error reporting
        workers: Number of parallel workers. 1 = sequential (default).

    Returns:
        Dict mapping variant_key -> VariantGroup
    """
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}

    if not target_set:
        if error_callback:
            error_callback("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    idxml_disambig = disambiguate_paths(idxml_paths)
    mzml_path_strs = [str(p) for p in mzml_paths] if mzml_paths else []

    total_files = len(idxml_paths)
    variant_groups: Dict[str, VariantGroup] = {}
    total_psms_found = 0

    tasks = [
        (p, idxml_disambig.get(p, Path(p).name), target_set, mzml_path_strs)
        for p in idxml_paths
    ]

    if workers > 1:
        logger.info(f"Parallel idXML search: {total_files} files, {workers} workers")

    def _on_result(task, result):
        nonlocal total_psms_found
        _, worker_matches, worker_errors = result
        for variant_key, match in worker_matches:
            mod_key = variant_key.split("|", 1)[1] if "|" in variant_key else ""
            if variant_key not in variant_groups:
                variant_groups[variant_key] = VariantGroup(
                    base_sequence=match.base_sequence,
                    modification_key=mod_key,
                    display_label=format_inline_sequence(match.base_sequence, mod_key),
                )
            variant_groups[variant_key].matches.append(match)
            total_psms_found += 1
        for error_message in worker_errors:
            if error_callback:
                error_callback(error_message)
            logger.error(error_message)

    def _on_error(task, exc):
        error_message = f"Worker crashed processing {Path(task[0]).name}: {exc}"
        if error_callback:
            error_callback(error_message)
        logger.error(error_message)

    def _on_progress(completed, total):
        if progress_callback:
            progress_callback(completed, total, "Searching idXML files")

    run_tasks(
        tasks=tasks,
        worker=_search_single_idxml,
        workers=workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
    )

    if progress_callback:
        progress_callback(total_files, total_files, "Search complete")

    logger.info(f"Found {total_psms_found} PSMs across {len(variant_groups)} variants")

    return variant_groups


def get_unique_peptides_from_idxmls(
    idxml_paths: List[str],
    limit: int = 100,
    progress_callback: Optional[callable] = None,
    workers: int = 1,
) -> List[str]:
    """
    Extract unique peptide sequences from idXML files.

    Args:
        idxml_paths: List of idXML file paths
        limit: Maximum number of unique peptides to return
        progress_callback: Optional callback for progress
        workers: Number of parallel workers. 1 = sequential (default).

    Returns:
        List of unique base sequences (sorted alphabetically)
    """
    unique_peptides: Set[str] = set()
    total_files = len(idxml_paths)

    if workers > 1:
        logger.info(f"Parallel peptide scan: {total_files} files, {workers} workers")

    def _on_result(task, result):
        _, file_peptides = result
        unique_peptides.update(file_peptides)

    def _on_error(task, exc):
        logger.warning(f"Worker error scanning peptides: {exc}")

    def _on_progress(completed, total):
        if progress_callback:
            progress_callback(completed, total, "Scanning idXML files")

    run_tasks(
        tasks=idxml_paths,  # single-arg → auto-wrapped by run_tasks
        worker=_scan_single_idxml_peptides,
        workers=workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
        stop_when=lambda: len(unique_peptides) >= limit,
    )

    if progress_callback:
        progress_callback(total_files, total_files, "Scan complete")

    return sorted(list(unique_peptides))[:limit]


def summarize_search_results(variant_groups: Dict[str, VariantGroup]) -> pl.LazyFrame:
    """
    Create a variant-level summary LazyFrame for the openms_insight Table component.

    One row per unique peptide variant (base_sequence + modification combination).

    Args:
        variant_groups: Dict of variant_key -> VariantGroup

    Returns:
        Polars LazyFrame with columns:
        - id (Int32): Row identifier required by openms_insight Table
        - base_sequence (String): Unmodified amino acid sequence
        - modification (String): Modification key or "(unmodified)"
        - display_label (String): Inline-annotated sequence (e.g. "PEPTM[Ox]DEK")
        - psm_count (Int32): Number of PSMs found for this variant
        - idxml_files (Int32): Number of distinct idXML files with hits
        - mzml_resolved (Int32): PSMs with a resolved mzML path
        - unresolved (Int32): PSMs where mzML could not be resolved
    """
    rows = []
    for row_id, (variant_key, variant_group) in enumerate(variant_groups.items()):
        idxml_files = set(m.idxml_file for m in variant_group.matches)
        unresolved = sum(1 for m in variant_group.matches if not m.mzml_filename)
        rows.append({
            "id": row_id,
            "base_sequence": variant_group.base_sequence,
            "modification": variant_group.modification_key or "(unmodified)",
            "display_label": variant_group.display_label,
            "psm_count": len(variant_group.matches),
            "idxml_files": len(idxml_files),
            "mzml_resolved": len(variant_group.matches) - unresolved,
            "unresolved": unresolved,
        })

    if not rows:
        schema = {
            "id": pl.Int32, "base_sequence": pl.String, "modification": pl.String,
            "display_label": pl.String, "psm_count": pl.Int32,
            "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
        }
        return pl.DataFrame(schema=schema).lazy()

    return pl.DataFrame(rows).with_columns([
        pl.col("id").cast(pl.Int32),
        pl.col("psm_count").cast(pl.Int32),
        pl.col("idxml_files").cast(pl.Int32),
        pl.col("mzml_resolved").cast(pl.Int32),
        pl.col("unresolved").cast(pl.Int32),
    ]).lazy()
