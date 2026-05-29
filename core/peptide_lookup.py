"""Peptide-lookup orchestration layer, dispatching the tasks of looking up target peptides
 across idXML files and summarising results for the UI; It uses core.helpers.run_tasks
to parallelize file-wise workersacross multiple processes.
(called by cli.main and streamlit_app.views.load_annotate).
"""
import logging
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

import polars as pl
from psm_utils import PSMList
from psm_utils.io.idxml import IdXMLReader

from core.helpers import run_tasks

logger = logging.getLogger(__name__)


# ============================================================================
# idXML scanning workers (dispatched across processes by run_tasks)
# ============================================================================

def count_psms_per_peptide(psm_list: PSMList) -> Counter:
    """Count PSMs per base sequence in an already-loaded PSMList.
    Returns a Counter keyed by uppercased unmodified peptide sequence,
    e.g. Counter({"PEPTIDE": 42}). Used by the Setup page to populate the
    "Pick peptides from idXML" tab.
    """
    return Counter(
        psm.peptidoform.sequence.upper()
        for psm in psm_list
        if psm.peptidoform.sequence
    )


def _count_psms_in_idxml(idxml_path: str) -> Counter:
    """Worker dispatched per idXML by count_psms_in_idxmls: load PSMs into a
    PSMList, then delegate the counting to count_psms_per_peptide.
    On parse errors, returns counts for whatever was read before the error.
    """
    psms: list = []
    try:
        for psm in IdXMLReader(idxml_path):
            psms.append(psm)
    except Exception as error:
        logger.warning(f"_count_psms_in_idxml partial read for {Path(idxml_path).name}: {error}")
    return count_psms_per_peptide(PSMList(psm_list=psms))


def lookup_peptides_in_single_idxml(
    idxml_path: str, idxml_name: str,
    target_set: Set[str], mzml_by_stem: Dict[str, str],
) -> Tuple[str, list, list]:
    """Scan one idXML file for target peptides; return the matching PSMs
     with extracted metadata attached (name, matches, errors).
    """
    # Open psm_utils' idxml-reader, which yields PSMs one at a time
    matches: list = []
    try:
        reader = IdXMLReader(idxml_path)
    # Catch file-level errors, returning empty results with an error message instead of crashing the whole lookup
    except Exception as load_error:
        return idxml_name, matches, [f"No peptide identifications in {idxml_name}: {load_error}"]

    # Fill a list with all PSMs, to surface parse errors as partial results
    psms: list = []
    parse_error: str | None = None
    try:
        for psm in reader:
            psms.append(psm)
    except Exception as error:
        parse_error = str(error)

    # Filter & annotate: keep target hits, stamp resolution metadata on each PSM
    for psm in psms:
        # Skip if no sequence or not in target set
        if psm.peptidoform.sequence not in target_set:
            continue
        # Annotate with idXML file
        psm.metadata["idxml_file"] = idxml_name
        psm.metadata["idxml_resolved_path"] = idxml_path
        # Annotate with mzML path if possible (match by run name, case-insensitive)
        if psm.run:
            resolved = mzml_by_stem.get(psm.run.lower(), "")
            if resolved:
                psm.metadata["mzml_resolved_path"] = resolved
        # Append modified sequence and metadata-enriched PSM to results (later becomes the variant_key)
        matches.append((psm.peptidoform.modified_sequence, psm))

    # Surface partial-read errors so the orchestrator can warn the user
    errors = []
    if parse_error is not None:
        errors.append(
            f"{idxml_name}: stopped after {len(psms)} PSMs — parse error: {parse_error}. "
            f"Try _pout_filtered.idXML for complete results."
        )
    return idxml_name, matches, errors


# ============================================================================
# Orchestration layer (called from cli.main and streamlit_app)
# ============================================================================

def count_psms_in_idxmls(
    idxml_paths: List[str],
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> "Counter[str]":
    """Small dispatcher run before peptide lookup to count by base sequence (ignoring modifications)
    across idXMLs; The result populates the "Pick peptides from idXML" tab in Streamlit
    keyed by base sequence (drives the Setup peptide-preview).
    """
    # Setup PSM counter, progress and error callbacks
    merged: Counter = Counter()
    total = len(idxml_paths)
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Dispatch the PSM counter worker across idXML paths,  and surfacing errors via callback;
    for completed, (task, ok, result) in enumerate(run_tasks(
        tasks=idxml_paths, worker=_count_psms_in_idxml, workers=workers,
    ), start=1):
        # On success, merge the per-file counters into the global counter
        if ok:
            merged.update(result)
        # On failure, surface errors via callback
        else:
            on_error(f"{Path(task).name}: {result}")
        progress(completed, total, "Counting PSMs per peptide")

    progress(total, total, "Count complete")
    return merged


def lookup_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> Dict[str, PSMList]:
    """Look up target_peptides across all idXMLs; return {variant_key: PSMList}
    with per-PSM resolution metadata stamped.
    """
    # Normalise inputs: target peptides are stripped, uppercased, and deduped into a set
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)
    if not target_set:
        on_error("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    # Build tasks: pack each idXML into the 4-tuple (path, name, target_set, mzml_by_stem) lookup_peptides_in_single_idxml expects
    mzml_by_stem = {Path(path).stem.lower(): path for path in (mzml_paths or [])}
    tasks = [(path, Path(path).name, target_set, mzml_by_stem) for path in idxml_paths]
    total = len(idxml_paths)
    variant_groups: Dict[str, PSMList] = {}

    # Dispatch  task stream tuples (task, ok, result) from run_tasks
    for completed, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=lookup_peptides_in_single_idxml, workers=workers), start=1,
    ):
        # On success, fold matches into variant groups (PSMList grouped by modified sequence)
        if ok:
            _, matches, errors = result
            for variant_key, psm in matches:
                variant_groups.setdefault(variant_key, PSMList(psm_list=[])).append(psm)
        # On failure, surface errors via callback
        else:
            errors = [f"Worker crashed processing {Path(task[0]).name}: {result}"]
        for message in errors:
            on_error(message)
        progress(completed, total, "Looking up peptides in idXML files")

    progress(total, total, "Lookup complete")
    total_psms = sum(len(g) for g in variant_groups.values())
    logger.info(f"Found {total_psms} PSMs across {len(variant_groups)} variants")
    return variant_groups


# LazyFrame schema for the per-variant summary consumed by the openms_insight Table
_SUMMARY_SCHEMA = {
    "id": pl.UInt32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}

def summarize_lookup_results(variant_groups: Dict[str, PSMList]) -> pl.LazyFrame:
    """Turn the variant groups returned by lookup_peptides_in_idxmls into a summary LazyFrame with one row per variant for display
    """
    # Flatten variant groups into one row per PSM, carrying variant & resolution metadata
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

    # Return empty LazyFrame with correct schema if no results are found
    if not rows:
        return pl.DataFrame(schema=_SUMMARY_SCHEMA).lazy()

    # Aggregate groups by variant (base sequence + modification) and count PSMs,
    # idXML files, and resolution status for summary display in the GUI
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