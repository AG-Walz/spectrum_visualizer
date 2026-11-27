"""Load & Annotate action layer — called by load_annotate_page.py (and the Setup peptide picker).

The analysis runs as a background job (streamlit_app.utils.job_registry), so it continues
while the user changes pages: index the idXMLs (only files not indexed yet), look up the
target peptides, then stream, annotate and compare their spectra, and write the results to
the workspace analysis cache. Every page reloads that cache when the job finishes. Previous
results are replaced only when a run succeeds.
"""
from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from pathlib import Path

import streamlit as st

from streamlit_app.utils.job_registry import AnalysisJob, JobFailed, job_key_for, start_job
from streamlit_app.utils.progress_tracker import progress_status
from streamlit_app.utils.state import IDXML_INDEX_DIR, get_param, save_params
from streamlit_app.utils.page_state import analysis_settings, write_analysis_cache
from core.defaults import DEFAULTS
from core.helpers import format_time
from core.dataframes import build_psm_table, build_visualization_frames
from core.pipeline import (
    build_annotation_kwargs,
    build_comparison_frame,
    build_koina_kwargs,
    run_streaming_annotation,
)
from core.idxml_index import ensure_idxml_index, load_indexed_psms, missing_from_index
from core.peptide_lookup import find_mzml_files, lookup_peptides

logger = logging.getLogger(__name__)


def index_idxmls(idxmls: list[str], workers: int) -> tuple[list, list[str]]:
    """Make sure every idXML has an up-to-date index; return (index_files, errors).
    Only files not indexed yet are parsed, under a live progress bar; when all are indexed
    this returns at once without drawing anything. Used by the Setup peptide picker.
    """
    errors: list[str] = []
    missing = missing_from_index(idxmls, IDXML_INDEX_DIR)
    if not missing:
        return ensure_idxml_index(idxmls, IDXML_INDEX_DIR, error_callback=errors.append), errors
    workers = max(1, min(workers, len(missing)))
    with progress_status(f"Indexing {len(missing)} idXML file(s) ({workers} workers)...") as tracker:
        index_files = ensure_idxml_index(
            idxmls, IDXML_INDEX_DIR, workers=workers,
            progress_callback=lambda done, total, _message: tracker.update(done, total, units="files"),
            error_callback=errors.append,
        )
        if errors:
            tracker.finalize(f"Indexed {len(missing)} idXML file(s) with {len(errors)} warnings",
                             state="error")
        else:
            tracker.finalize(f"Indexed {len(missing)} idXML file(s)")
    return index_files, errors


def lookup_from_index(index_files: list, peptides: list[str], mzmls: list[str]) -> dict:
    """Return {peptidoform: PSMList} for the target peptides, read from the idXML index.
    Only the target-peptide rows are unpickled; grouping and mzML resolution as in core.
    """
    return lookup_peptides(load_indexed_psms(index_files, peptides), peptides, mzml_paths=mzmls)


def add_missing_mzmls(params: dict, run_names: list[str], idxmls: list[str]) -> str:
    """Find mzML files for runs that no selected mzML matches, add the unambiguous ones to the
    selection (and params.json), and return a one-line report for the page.
    Searches the Setup mzML folder and, for each idXML, the folder above its own (pipeline
    layouts often keep idXML and mzML in sibling folders).
    """
    search_dirs = [get_param(params, "mzml_browser_dir", "")] + [
        str(Path(path).parent.parent) for path in idxmls]
    found, complete = find_mzml_files(run_names, [d for d in search_dirs if d and Path(d).is_dir()])
    selected = list(get_param(params, "selected_mzmls", []))
    added = [paths[0] for paths in found.values() if len(paths) == 1 and paths[0] not in selected]
    ambiguous = sorted(run for run, paths in found.items() if len(paths) > 1)
    missing = sorted(run for run, paths in found.items() if not paths)
    if added:
        st.session_state["selected_mzmls"] = selected + added
        save_params(params)

    report = [f"Added {len(added)} mzML file(s)."]
    if ambiguous:
        report.append(f"{len(ambiguous)} run(s) match several files; select the right one on "
                      f"the Setup page: {', '.join(ambiguous[:5])}.")
    if missing:
        report.append(f"No file found for {len(missing)} run(s): {', '.join(missing[:5])}.")
    if not complete:
        report.append("The search stopped early because the folders are very large; "
                      "point the Setup mzML browser at a smaller folder.")
    return " ".join(report)


def show_messages(messages: list[str], label: str) -> None:
    """Collapse warnings into one expander (first 20 shown) so they never flood the page."""
    if not messages:
        return
    with st.expander(f"{len(messages)} {label}", expanded=False):
        for message in messages[:20]:
            st.caption(message)
        if len(messages) > 20:
            st.caption(f"... and {len(messages) - 20} more")


def analysis_job_key() -> str:
    """Job-registry key of the active workspace (one analysis job per workspace)."""
    return job_key_for(st.session_state.workspace)


@dataclass(frozen=True)
class AnalysisInputs:
    """Everything the job thread needs, resolved in the script thread (the job thread must
    not read session state or widgets)."""
    workspace: Path
    idxmls: list[str]
    mzmls: list[str]
    peptides: list[str]
    workers: int
    reference_substring: str
    annotation_kwargs: dict
    koina_kwargs: dict
    settings: dict


def start_analysis(params: dict, idxmls: list[str], peptides: list[str]) -> bool:
    """Start the analysis as a background job; False when this workspace already runs one."""
    getter = lambda key, default: get_param(params, key, default)
    inputs = AnalysisInputs(
        workspace=Path(st.session_state.workspace),
        idxmls=list(idxmls),
        mzmls=list(get_param(params, "selected_mzmls", [])),
        peptides=list(peptides),
        workers=int(get_param(params, "workers", DEFAULTS.workers)),
        reference_substring=get_param(params, "reference_substring", ""),
        annotation_kwargs=build_annotation_kwargs(getter),
        koina_kwargs=build_koina_kwargs(getter),
        settings=analysis_settings(params),
    )
    return start_job(analysis_job_key(), lambda job: _run_analysis_job(job, inputs)) is not None


def _run_analysis_job(job: AnalysisJob, inputs: AnalysisInputs) -> str:
    """Job-thread body: index, look up, stream + annotate, then write the analysis cache.
    Reports progress through job.progress (which also raises JobCancelled on cancel) and
    returns the one-line summary. Raises JobFailed when there is nothing to store.
    """
    # Step 1: index the idXMLs not indexed yet, then look up the target peptides
    missing = missing_from_index(inputs.idxmls, IDXML_INDEX_DIR)
    job.progress("Indexing idXML files", 0, len(missing), "files")
    index_files = ensure_idxml_index(
        inputs.idxmls, IDXML_INDEX_DIR, workers=max(1, min(inputs.workers, len(missing) or 1)),
        progress_callback=lambda done, total, _message: job.progress(
            "Indexing idXML files", done, total, "files"),
        error_callback=job.warn,
    )
    job.progress("Looking up peptides")
    peptidoform_psms = lookup_from_index(index_files, inputs.peptides, inputs.mzmls)
    if not peptidoform_psms:
        raise JobFailed("None of the target peptides is in the selected idXML files. "
                        "The previous results (if any) are unchanged.")

    # Step 2: deep-copy peptidoforms (workers mutate metadata) and build the PSM table
    peptidoform_psms_copy = copy.deepcopy(peptidoform_psms)
    psm_frame, psm_lookup = build_psm_table(peptidoform_psms_copy)
    total_psms = sum(len(psm_list) for psm_list in peptidoform_psms_copy.values())
    # File phase total: the number of mzMLs the lookup resolved PSMs to
    total_files = max(1, len({
        psm.metadata.get("mzml_resolved_path")
        for psm_list in peptidoform_psms_copy.values() for psm in psm_list
        if psm.metadata.get("mzml_resolved_path")
    }))
    files_done = [0]

    def on_file_complete(name, loaded, not_found) -> None:
        files_done[0] += 1
        logger.info(f"Completed {Path(name).name}: {loaded} loaded, {not_found} not found")
        job.progress("Reading mzML files", min(files_done[0], total_files), total_files, "files")

    def on_error(error: str) -> None:
        logger.error(f"Streaming error: {error}")
        job.warn(error)

    # Step 3: stream spectra (reading files), then annotate and score each query
    job.progress("Reading mzML files", 0, total_files, "files")
    annotated_queries, reference_store, comparison_rows = run_streaming_annotation(
        peptidoform_psms=peptidoform_psms_copy,
        mzml_paths=inputs.mzmls,
        reference_substring=inputs.reference_substring,
        psm_lookup=psm_lookup,
        annotation_kwargs=inputs.annotation_kwargs,
        koina_kwargs=inputs.koina_kwargs,
        workers=inputs.workers,
        on_progress=lambda count: job.progress("Annotating spectra", count, total_psms, "spectra"),
        on_file_complete=on_file_complete,
        on_error=on_error,
    )
    if not annotated_queries:
        raise JobFailed("No spectra could be loaded. The previous results (if any) are unchanged.")

    # Step 4: build the visualization frames (query AND reference spectra, #94) and swap the
    # new results into the analysis cache. Last cancel point is here, before anything is written.
    job.progress("Saving results")
    frames = build_visualization_frames(
        annotated_queries, psm_frame=psm_frame, reference_store=reference_store)
    write_analysis_cache(inputs.workspace, {
        "psm_frame": psm_frame,
        "peaks_frame": frames["peaks"],
        "sequence_frame": frames["sequence"],
        "heatmap_frame": frames["heatmap"],
        "comparison_frame": build_comparison_frame(comparison_rows),
        "psm_lookup": psm_lookup,
        "annotated_queries": {
            psm_id: (query_psm, annotated_peaks_frame)
            for psm_id, query_psm, annotated_peaks_frame in annotated_queries
        },
        "reference_store": reference_store,
        "peptidoform_psms": peptidoform_psms,
    }, inputs.settings)
    return (f"Loaded {len(annotated_queries)} spectra and {len(comparison_rows)} "
            f"comparisons in {format_time(job.elapsed)}.")
