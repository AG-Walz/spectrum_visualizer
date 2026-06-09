"""Peptide-lookup layer: load idXML files once into a PSMList (parallel file I/O via
core.helpers.run_tasks), then scan that in-memory PSMList — count PSMs per peptide for
the Setup picker, and group target-peptide hits into peptidoform_psms for the pipeline.
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

def _load_idxml_worker(idxml_path: str, idxml_name: str) -> Tuple[str, list, list]:
    """Per-file worker that load_idxmls maps over via run_tasks (one call per idXML,
    dispatched sequentially or across processes). Reads every PSM once and stamps
    idXML provenance onto it. Returns (name, psms, errors); on a parse error the
    PSMs read before the failure are still returned as a partial result.
    Not a public single-file loader — load_idxmls is the entry point for any file count.
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
    # idxml_source is the parent directory name — it disambiguates same-basename
    # idXMLs loaded from different source folders (e.g. illumina/ vs pac/) so their
    # PSMs stay distinguishable downstream and resolve to their own co-located mzML.
    idxml_source = Path(idxml_path).parent.name
    for psm in psms:
        psm.metadata["idxml_file"] = idxml_name
        psm.metadata["idxml_resolved_path"] = idxml_path
        psm.metadata["idxml_source"] = idxml_source
        # PSMList-native migration (Step 1): mirror the same provenance into the
        # native PSM.provenance_data dict. Additive dual-write — no consumer reads
        # it yet; psm.metadata stays the source of truth for now.
        psm.provenance_data = {
            **(psm.provenance_data or {}),
            "idxml_file": idxml_name,
            "idxml_resolved_path": idxml_path,
            "idxml_source": idxml_source,
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
    This is the ONLY place idXML files are parsed. Both peptide counting (count_psms_per_peptide)
    and target lookup (lookup_peptides) operate on the returned in-memory PSMList, so the same file
    is never read twice — callers cache it and reuse it across the Setup and Load & Annotate stages.
    """
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Build tasks: pack each idXML into the (path, name) tuple _load_idxml_worker expects
    tasks = [(path, Path(path).name) for path in idxml_paths]
    total = len(idxml_paths)
    psms: list = []

    # Dispatch the loader worker across idXML paths, collecting PSMs and surfacing errors
    for completed, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=_load_idxml_worker, workers=workers), start=1,
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


def _index_mzml_paths(
    mzml_paths: List[str],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, List[str]]]:
    """Index selected mzML paths for source-aware run resolution.

    Returns (by_dir, by_stem): by_dir maps a parent directory to its
    {stem_lower: path}; by_stem maps a stem to every path sharing it (so the
    fallback can tell a unique stem from one that's ambiguous across folders).
    """
    by_dir: Dict[str, Dict[str, str]] = {}
    by_stem: Dict[str, List[str]] = {}
    for path in mzml_paths:
        parsed = Path(path)
        stem = parsed.stem.lower()
        by_dir.setdefault(str(parsed.parent), {})[stem] = path
        by_stem.setdefault(stem, []).append(path)
    return by_dir, by_stem


def _resolve_mzml_path(
    psm,
    mzml_by_dir: Dict[str, Dict[str, str]],
    mzml_by_stem: Dict[str, List[str]],
) -> str:
    """Resolve a PSM's run to a selected mzML path, preferring the mzML co-located
    with the PSM's source idXML.

    Falls back to a global stem match only when the idXML's own directory has no
    match. If the stem is ambiguous across multiple directories and none is
    co-located, it refuses to guess (returns "" and warns) rather than risk
    loading another source's spectra.
    """
    run = psm.run.lower()
    source_dir = str(Path(psm.metadata.get("idxml_resolved_path", "")).parent)
    same_dir = mzml_by_dir.get(source_dir, {})
    if run in same_dir:
        return same_dir[run]
    candidates = mzml_by_stem.get(run, [])
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        logger.warning(
            "mzML run %r is ambiguous across %d directories and no mzML is co-located "
            "with idXML %r; leaving it unresolved. Co-locate the mzML in the same "
            "folder as its idXML to disambiguate.",
            psm.run, len(candidates), psm.metadata.get("idxml_file", ""),
        )
    return ""


def lookup_peptides(
    psm_list: PSMList,
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
) -> Dict[str, PSMList]:
    """Group the target-peptide PSMs from an already-loaded PSMList into
    {peptidoform: PSMList}, resolving and stamping each hit's mzML path.
    Pure in-memory work: no idXML is re-read (see load_idxmls). Hits are selected
    with a PSMList boolean mask, then grouped by modified sequence.
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

    # Resolve each hit's mzML path (depends on the current selection) and group by
    # peptidoform. Resolution is source-aware: an idXML and its mzML are paired within
    # the same source directory, so a PSM resolves to the mzML co-located with its own
    # idXML. This keeps same-basename files from different folders (illumina/ vs pac/)
    # pointing at their own spectra instead of one silently shadowing the other.
    mzml_by_dir, mzml_by_stem = _index_mzml_paths(mzml_paths or [])
    peptidoform_psms: Dict[str, PSMList] = {}
    for psm in hits:
        # Annotate with mzML path if possible (match by run name, case-insensitive)
        if psm.run:
            resolved = _resolve_mzml_path(psm, mzml_by_dir, mzml_by_stem)
            if resolved:
                psm.metadata["mzml_resolved_path"] = resolved
                psm.provenance_data = {**(psm.provenance_data or {}), "mzml_resolved_path": resolved}
        # Group under the modified sequence (the peptidoform)
        peptidoform_psms.setdefault(psm.peptidoform.modified_sequence, PSMList(psm_list=[])).append(psm)

    total_psms = sum(len(group) for group in peptidoform_psms.values())
    logger.info(f"Found {total_psms} PSMs across {len(peptidoform_psms)} peptidoforms")
    return peptidoform_psms


# LazyFrame schema for the per-peptidoform summary consumed by the openms_insight Table
_SUMMARY_SCHEMA = {
    "id": pl.UInt32, "base_sequence": pl.String, "modification": pl.String,
    "display_label": pl.String, "psm_count": pl.Int32,
    "idxml_files": pl.Int32, "sources": pl.Int32,
    "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
}

def summarize_lookup_results(peptidoform_psms: Dict[str, PSMList]) -> pl.LazyFrame:
    """Turn the peptidoform PSMs returned by lookup_peptides into a summary LazyFrame with one row per peptidoform for display
    """
    # Flatten peptidoform PSMs into one row per PSM, carrying peptidoform & resolution metadata
    rows = []
    for peptidoform, psm_list in peptidoform_psms.items():
        if not psm_list:
            continue
        base = psm_list[0].peptidoform.sequence.upper()
        modification = peptidoform if psm_list[0].peptidoform.is_modified else "(unmodified)"
        rows.extend({
            "base_sequence": base,
            "modification": modification,
            "display_label": peptidoform,
            "idxml_file": psm.metadata.get("idxml_file", ""),
            "source": psm.metadata.get("idxml_source", ""),
            "mzml_resolved": bool(psm.metadata.get("mzml_resolved_path")),
        } for psm in psm_list)

    # Return empty LazyFrame with correct schema if no results are found
    if not rows:
        return pl.DataFrame(schema=_SUMMARY_SCHEMA).lazy()

    # Aggregate groups by peptidoform (base sequence + modification) and count PSMs,
    # idXML files, and resolution status for summary display in the GUI
    return (
        pl.DataFrame(rows).lazy()
        .group_by(["base_sequence", "modification", "display_label"], maintain_order=True)
        .agg([
            pl.len().alias("psm_count"),
            # Count distinct (source, file) pairs so same-basename idXMLs from
            # different source folders (illumina/ vs pac/) aren't collapsed into one.
            pl.struct("source", "idxml_file").n_unique().alias("idxml_files"),
            pl.col("source").n_unique().alias("sources"),
            pl.col("mzml_resolved").sum().alias("mzml_resolved"),
            pl.col("mzml_resolved").not_().sum().alias("unresolved"),
        ])
        .with_row_index("id")
        .with_columns([pl.col(c).cast(dtype) for c, dtype in _SUMMARY_SCHEMA.items()])
    )