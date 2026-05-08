"""High-level peptide-search orchestration.

Coordinates parallel/sequential execution of the per-file workers in
:mod:`mspecviz.search.idxml_workers` and aggregates results into a
``Dict[str, PSMList]`` keyed on canonical ``Peptidoform.modified_sequence``.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import polars as pl
from psm_utils import PSMList

from mspecviz.common.task_dispatch import run_tasks
from mspecviz.search.idxml_workers import _scan_single_idxml_peptides, _search_single_idxml

logger = logging.getLogger(__name__)


def _build_mzml_lookup(mzml_paths: Optional[List[str]]) -> Dict[str, str]:
    """Map ``stem.lower() → str(path)`` (matches ``psm_utils.PSM.run`` keying)."""
    if not mzml_paths:
        return {}
    return {Path(path).stem.lower(): str(path) for path in mzml_paths}


def search_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[callable] = None,
    error_callback: Optional[callable] = None,
    workers: int = 1,
) -> Dict[str, PSMList]:
    """Search target peptides across idXML files; return ``Dict[variant_key, PSMList]``.

    Each PSM carries mspecviz-specific resolution metadata (``idxml_file``,
    ``idxml_resolved_path``, ``mzml_resolved_path``) in ``PSM.metadata``,
    populated by :func:`mspecviz.search.idxml_workers._search_single_idxml`.
    """
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    if not target_set:
        if error_callback:
            error_callback("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    mzml_by_stem = _build_mzml_lookup(mzml_paths)
    tasks = [(path, Path(path).name, target_set, mzml_by_stem) for path in idxml_paths]
    total_files = len(idxml_paths)

    if workers > 1:
        logger.info(f"Parallel idXML search: {total_files} files, {workers} workers")

    variant_groups: Dict[str, PSMList] = {}
    total_psms = 0
    completed = 0
    for task, ok, result in run_tasks(tasks=tasks, worker=_search_single_idxml, workers=workers):
        completed += 1
        if not ok:
            message = f"Worker crashed processing {Path(task[0]).name}: {result}"
            if error_callback:
                error_callback(message)
            logger.error(message)
        else:
            _, matches, errors = result
            for variant_key, psm in matches:
                variant_groups.setdefault(variant_key, PSMList(psm_list=[])).append(psm)
                total_psms += 1
            for message in errors:
                if error_callback:
                    error_callback(message)
                logger.error(message)
        if progress_callback:
            progress_callback(completed, total_files, "Searching idXML files")

    if progress_callback:
        progress_callback(total_files, total_files, "Search complete")
    logger.info(f"Found {total_psms} PSMs across {len(variant_groups)} variants")
    return variant_groups


def get_unique_peptides_from_idxmls(
    idxml_paths: List[str],
    limit: int = 100,
    progress_callback: Optional[callable] = None,
    workers: int = 1,
) -> List[str]:
    """Return up to *limit* unique base sequences across *idxml_paths* (sorted)."""
    unique_peptides: Set[str] = set()
    total_files = len(idxml_paths)

    if workers > 1:
        logger.info(f"Parallel peptide scan: {total_files} files, {workers} workers")

    completed = 0
    # Single-arg tasks: run_tasks auto-wraps each ``idxml_path`` into a 1-tuple.
    for _task, ok, result in run_tasks(
        tasks=idxml_paths, worker=_scan_single_idxml_peptides, workers=workers,
    ):
        completed += 1
        if ok:
            _, file_peptides = result
            unique_peptides.update(file_peptides)
        else:
            logger.warning(f"Worker error scanning peptides: {result}")
        if progress_callback:
            progress_callback(completed, total_files, "Scanning idXML files")
        if len(unique_peptides) >= limit:
            break  # cancels remaining futures via run_tasks' try/finally

    if progress_callback:
        progress_callback(total_files, total_files, "Scan complete")
    return sorted(unique_peptides)[:limit]


_SUMMARY_SCHEMA = {
    "id": pl.Int32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}


def summarize_search_results(variant_groups: Dict[str, PSMList]) -> pl.LazyFrame:
    """Per-variant summary LazyFrame for the openms_insight Table component.

    Aggregation runs in Polars (``group_by``/``agg``); per-PSM rows are built
    once in Python from the PSM metadata and the variant_key (which is the
    canonical ``Peptidoform.modified_sequence``).
    """
    rows = []
    for variant_key, psm_list in variant_groups.items():
        if not psm_list:
            continue
        base_sequence = psm_list[0].peptidoform.sequence.upper()
        is_modified = variant_key != base_sequence
        for psm in psm_list:
            rows.append({
                "variant_key": variant_key,
                "base_sequence": base_sequence,
                "modification": variant_key if is_modified else "(unmodified)",
                "display_label": variant_key,
                "idxml_file": psm.metadata.get("idxml_file", ""),
                "mzml_resolved": bool(psm.metadata.get("mzml_resolved_path")),
            })

    if not rows:
        return pl.DataFrame(schema=_SUMMARY_SCHEMA).lazy()

    return (
        pl.DataFrame(rows).lazy()
        .group_by(["base_sequence", "modification", "display_label"], maintain_order=True)
        .agg([
            pl.len().alias("psm_count"),
            pl.col("idxml_file").n_unique().alias("idxml_files"),
            pl.col("mzml_resolved").sum().alias("mzml_resolved"),
            pl.col("mzml_resolved").not_().sum().alias("unresolved"),
        ])
        .with_row_index("id")
        .with_columns([pl.col(c).cast(dtype) for c, dtype in _SUMMARY_SCHEMA.items()])
    )
