"""Peptide-lookup layer: load idXML files once into a PSMList (parallel file I/O via
core.helpers.run_tasks), then scan that in-memory PSMList — count PSMs per peptide for
the Setup picker, and group target-peptide hits into variant_groups for the pipeline.
Loading is separated from scanning so the same files are never parsed twice.
(called by cli.main and streamlit_app.views.load_annotate / widgets.peptide_input).
"""
import logging
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np
import polars as pl
from psm_utils import PSMList
from psm_utils.io.idxml import IdXMLReader

from core.helpers import run_tasks

logger = logging.getLogger(__name__)


# ============================================================================
# idXML loading (parallel file I/O — each file is parsed exactly once)
# ============================================================================

def _load_single_idxml(idxml_path: str, idxml_name: str) -> Tuple[str, list, list]:
    """Worker dispatched per idXML by load_idxmls: read every PSM once and stamp
    idXML provenance onto it. Returns (name, psms, errors); on a parse error the
    PSMs read before the failure are still returned as a partial result.
    """
    # Open psm_utils' idxml-reader, which yields PSMs one at a time
    try:
        reader = IdXMLReader(idxml_path)
    # Catch file-level errors, returning empty results with an error message instead of crashing the whole load
    except Exception as load_error:
        return idxml_name, [], [f"No peptide identifications in {idxml_name}: {load_error}"]

    # Fill a list with all PSMs, to surface parse errors as partial results
    psms: list = []
    parse_error: str | None = None
    try:
        for psm in reader:
            psms.append(psm)
    except Exception as error:
        parse_error = str(error)

    # Stamp idXML provenance on every PSM. mzML resolution is deliberately NOT done
    # here — it depends on the user's current mzML selection and is applied later by
    # lookup_peptides, so this loaded PSMList stays reusable across mzML selections.
    for psm in psms:
        psm.metadata["idxml_file"] = idxml_name
        psm.metadata["idxml_resolved_path"] = idxml_path
        # PSMList-native migration (Step 1): mirror the same provenance into the
        # native PSM.provenance_data dict. Additive dual-write — no consumer reads
        # it yet; psm.metadata stays the source of truth for now.
        psm.provenance_data = {
            **(psm.provenance_data or {}),
            "idxml_file": idxml_name,
            "idxml_resolved_path": idxml_path,
        }

    # Surface partial-read errors so the loader can warn the user
    errors = []
    if parse_error is not None:
        errors.append(
            f"{idxml_name}: stopped after {len(psms)} PSMs — parse error: {parse_error}. "
            f"Try _pout_filtered.idXML for complete results."
        )
    return idxml_name, psms, errors


def load_idxmls(
    idxml_paths: List[str],
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
    workers: int = 1,
) -> PSMList:
    """Read every idXML once into a single PSMList, stamping idXML provenance per PSM.

    This is the ONLY place idXML files are parsed. Both peptide counting
    (count_psms_per_peptide) and target lookup (lookup_peptides) operate on the
    returned in-memory PSMList, so the same file is never read twice — callers
    cache the PSMList and reuse it across the Setup and Load & Annotate stages.
    """
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Build tasks: pack each idXML into the (path, name) tuple _load_single_idxml expects
    tasks = [(path, Path(path).name) for path in idxml_paths]
    total = len(idxml_paths)
    psms: list = []

    # Dispatch the loader worker across idXML paths, collecting PSMs and surfacing errors
    for completed, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=_load_single_idxml, workers=workers), start=1,
    ):
        # On success, accumulate this file's PSMs into the shared list
        if ok:
            _, file_psms, errors = result
            psms.extend(file_psms)
        # On failure, surface errors via callback
        else:
            errors = [f"Worker crashed processing {Path(task[0]).name}: {result}"]
        for message in errors:
            on_error(message)
        progress(completed, total, "Loading idXML files")

    progress(total, total, "Load complete")
    logger.info(f"Loaded {len(psms)} PSMs from {total} idXML files")
    return PSMList(psm_list=psms)


# ============================================================================
# In-memory scanning (operate on a loaded PSMList — no idXML is re-read)
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


def lookup_peptides(
    psm_list: PSMList,
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
) -> Dict[str, PSMList]:
    """Group the target-peptide PSMs from an already-loaded PSMList into
    {variant_key: PSMList}, resolving each hit's mzML path and stamping it.

    Pure in-memory work: no idXML is re-read (see load_idxmls). Target hits are
    selected with a PSMList boolean mask, then grouped by modified sequence.
    """
    # Normalise inputs: target peptides are stripped, uppercased, and deduped into a set
    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    if not target_set:
        logger.warning("No valid target peptides provided")
        return {}
    if not psm_list:
        return {}

    # Select target hits via a PSMList boolean mask (returns a new PSMList of matches)
    is_target = np.fromiter(
        ((psm.peptidoform.sequence or "").upper() in target_set for psm in psm_list),
        dtype=bool, count=len(psm_list),
    )
    hits = psm_list[is_target]

    # Resolve each hit's mzML path (depends on the current selection) and group by variant
    mzml_by_stem = {Path(path).stem.lower(): path for path in (mzml_paths or [])}
    variant_groups: Dict[str, PSMList] = {}
    for psm in hits:
        # Annotate with mzML path if possible (match by run name, case-insensitive)
        if psm.run:
            resolved = mzml_by_stem.get(psm.run.lower(), "")
            if resolved:
                psm.metadata["mzml_resolved_path"] = resolved
                psm.provenance_data = {**(psm.provenance_data or {}), "mzml_resolved_path": resolved}
        # Group under the modified sequence (the variant_key)
        variant_groups.setdefault(psm.peptidoform.modified_sequence, PSMList(psm_list=[])).append(psm)

    total_psms = sum(len(group) for group in variant_groups.values())
    logger.info(f"Found {total_psms} PSMs across {len(variant_groups)} variants")
    return variant_groups


# LazyFrame schema for the per-variant summary consumed by the openms_insight Table
_SUMMARY_SCHEMA = {
    "id": pl.UInt32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}

def summarize_lookup_results(variant_groups: Dict[str, PSMList]) -> pl.LazyFrame:
    """Turn the variant groups returned by lookup_peptides into a summary LazyFrame with one row per variant for display
    """
    # Flatten variant groups into one row per PSM, carrying variant & resolution metadata
    rows = []
    for variant_key, psm_list in variant_groups.items():
        if not psm_list:
            continue
        base = psm_list[0].peptidoform.sequence.upper()
        modification = variant_key if psm_list[0].peptidoform.is_modified else "(unmodified)"
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