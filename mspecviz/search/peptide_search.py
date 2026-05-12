"""High-level peptide-search orchestration.

Coordinates parallel/sequential execution of the per-file workers in
:mod:`mspecviz.search.idxml_workers` and aggregates results into a
``Dict[str, PSMList]`` keyed on canonical ``Peptidoform.modified_sequence``.
"""
import logging
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set

import polars as pl
from psm_utils import PSMList

from mspecviz.common.task_dispatch import run_tasks
from mspecviz.search.idxml_workers import _search_single_idxml, count_psms_per_peptide

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
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> Dict[str, PSMList]:
    """Search target peptides across idXML files; return ``Dict[variant_key, PSMList]``.

    Each PSM carries mspecviz-specific resolution metadata (``idxml_file``,
    ``idxml_resolved_path``, ``mzml_resolved_path``) in ``PSM.metadata``,
    populated by :func:`mspecviz.search.idxml_workers._search_single_idxml`.
    """
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)
    if not target_set:
        on_error("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    mzml_by_stem = _build_mzml_lookup(mzml_paths)
    tasks = [(path, Path(path).name, target_set, mzml_by_stem) for path in idxml_paths]
    total = len(idxml_paths)
    variant_groups: Dict[str, PSMList] = {}

    for completed, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=_search_single_idxml, workers=workers), start=1,
    ):
        if ok:
            _, matches, errors = result
            for variant_key, psm in matches:
                variant_groups.setdefault(variant_key, PSMList(psm_list=[])).append(psm)
        else:
            errors = [f"Worker crashed processing {Path(task[0]).name}: {result}"]
        for message in errors:
            on_error(message)
            logger.error(message)
        progress(completed, total, "Searching idXML files")

    progress(total, total, "Search complete")
    total_psms = sum(len(g) for g in variant_groups.values())
    logger.info(f"Found {total_psms} PSMs across {len(variant_groups)} variants")
    return variant_groups


def count_psms_in_idxmls(
    idxml_paths: List[str],
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> "Counter[str]":
    """Merge per-file PSM counts; ``{base_sequence: total_psms}`` across *idxml_paths*."""
    merged: Counter = Counter()
    total = len(idxml_paths)
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Single-arg tasks: run_tasks auto-wraps each ``idxml_path`` into a 1-tuple.
    for completed, (task, ok, result) in enumerate(run_tasks(
        tasks=idxml_paths, worker=count_psms_per_peptide, workers=workers,
    ), start=1):
        if ok:
            merged.update(result)
        else:
            on_error(f"{Path(task).name}: {result}")
        progress(completed, total, "Counting PSMs per peptide")

    progress(total, total, "Count complete")
    return merged


_SUMMARY_SCHEMA = {
    "id": pl.UInt32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}


def summarize_search_results(variant_groups: Dict[str, PSMList]) -> pl.LazyFrame:
    """Per-variant summary LazyFrame for the openms_insight Table component."""
    rows = []
    for variant_key, psm_list in variant_groups.items():
        if not psm_list:
            continue
        base = psm_list[0].peptidoform.sequence.upper()
        modification = variant_key if variant_key != base else "(unmodified)"
        rows.extend({
            "base_sequence": base,
            "modification": modification,
            "display_label": variant_key,
            "idxml_file": psm.metadata.get("idxml_file", ""),
            "mzml_resolved": bool(psm.metadata.get("mzml_resolved_path")),
        } for psm in psm_list)

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