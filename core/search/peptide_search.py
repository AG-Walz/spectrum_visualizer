"""Peptide-search orchestration (called by ``cli.main`` and ``streamlit_app.views.search_load``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set

import polars as pl
from psm_utils import PSMList

from core.common.task_dispatch import run_tasks
from core.search.idxml_workers import search_single_idxml, count_psms_per_peptide

logger = logging.getLogger(__name__)


def search_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> Dict[str, PSMList]:
    """Search ``target_peptides`` across all idXMLs; return ``{variant_key: PSMList}`` with per-PSM resolution metadata stamped."""

    # Normalise inputs: dedupe targets, wire optional callbacks
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)
    if not target_set:
        on_error("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    # Build tasks: pack each idXML into the 4-tuple search_single_idxml expects
    mzml_by_stem = {Path(path).stem.lower(): str(path) for path in (mzml_paths or [])}
    tasks = [(path, Path(path).name, target_set, mzml_by_stem) for path in idxml_paths]
    total = len(idxml_paths)
    variant_groups: Dict[str, PSMList] = {}

    # Dispatch & aggregate: group matches by variant, surface errors, tick progress
    for completed, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=search_single_idxml, workers=workers), start=1,
    ):
        if ok:
            _, matches, errors = result
            for variant_key, psm in matches:
                variant_groups.setdefault(variant_key, PSMList(psm_list=[])).append(psm)
        else:
            errors = [f"Worker crashed processing {Path(task[0]).name}: {result}"]
        for message in errors:
            on_error(message)
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
    """Merge per-file PSM counts across ``idxml_paths`` into a ``Counter`` keyed by base sequence (drives the Setup peptide-preview)."""

    # Setup PSM counter, progress and error callbacks
    merged: Counter = Counter()
    total = len(idxml_paths)
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Dispatch & merge: tick progress, surface errors, fold per-file counts in
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

# One row per variant: counts + resolution metadata aggregated from PSMs
_SUMMARY_SCHEMA = {
    "id": pl.UInt32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}


def summarize_search_results(variant_groups: Dict[str, PSMList]) -> pl.LazyFrame:
    """Build the per-variant summary LazyFrame consumed by the openms_insight Table."""

    # Flatten variant groups: one row per PSM, carrying variant & resolution metadata
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

    # Empty input: hand back a typed-empty frame so downstream schemas stay stable
    if not rows:
        return pl.DataFrame(schema=_SUMMARY_SCHEMA).lazy()

    # Group & aggregate: per-variant counts, then enforce schema
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