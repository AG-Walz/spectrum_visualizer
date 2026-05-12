"""
Search & Load pipeline helpers — peptide search, spectrum streaming, annotation.

These are the long-running operations triggered from ``page.py`` buttons.
They take ``params`` explicitly (not from module state) so the split remains
testable and re-runnable without reloading the view module.
"""
from __future__ import annotations

import copy
import logging
from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_app.utils.progress_tracker import progress_status
from streamlit_app.utils.state import get_param, get_state_key, save_params
from streamlit_app.utils.page_state import clear_analysis_cache, save_analysis_cache
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.helpers import format_time
from mspecviz.dataframes import build_psm_table, build_visualization_frames
from mspecviz.pipeline import (
    build_annotation_kwargs,
    build_koina_kwargs,
    run_streaming_annotation,
)
from mspecviz.search import search_peptides_in_idxmls

logger = logging.getLogger(__name__)


def _show_search_results_summary(results, error_messages, tracker, params: dict) -> None:
    """Finalize the status block with PSM/variant counts and any warnings."""
    # Empty path: mark the status block as an error and bail
    if not results:
        tracker.finalize("No matches found", state="error")
        return

    # Success path: report PSM/variant counts and persist params
    total_matches = sum(len(psm_list) for psm_list in results.values())
    tracker.finalize(f"Found {total_matches} PSMs across {len(results)} variants")
    save_params(params)

    # Warnings: surface unresolved-mzML count + collapse per-file errors
    unresolved = sum(
        1 for psm_list in results.values()
        for psm in psm_list
        if not psm.metadata.get("mzml_resolved_path")
    )
    if unresolved > 0:
        st.warning(f"{unresolved} PSMs could not be mapped to mzML files")
    if error_messages:
        with st.expander(f"{len(error_messages)} warnings", expanded=False):
            for message in error_messages[:10]:
                st.caption(message)


def run_peptide_search(params: dict, idxmls: list, peptides: list) -> None:
    """Execute peptide search."""
    mzmls = get_param(params, "selected_mzmls", [])
    error_messages: list[str] = []

    with progress_status("Searching peptides...") as tracker:
        results = search_peptides_in_idxmls(
            idxmls, peptides, mzml_paths=mzmls,
            progress_callback=lambda current, total, message: tracker.update(
                current, total, units="files", message=message,
            ),
            error_callback=error_messages.append,
            workers=get_param(params, "workers", DEFAULTS.workers),
        )
        st.session_state[get_state_key("variant_groups")] = results
        _show_search_results_summary(results, error_messages, tracker, params)


def _params_getter(params: dict):
    """Adapter so shared kwargs builders can read from Streamlit ``params``.

    Streamlit stores some keys under different names (e.g. ``selected_mzmls``
    in session vs ``mzml_files`` in run_settings); for the annotation/Koina
    kwargs the names already match.
    """
    return lambda key, default: get_param(params, key, default)


def _store_analysis_results(psm_frame, psm_lookup, annotated_queries,
                            reference_store, comparison_rows) -> None:
    """Persist all analysis results to session state and disk cache."""
    # Build visualization frames (peaks / sequence / heatmap) from annotations
    frames = build_visualization_frames(annotated_queries)

    # Stash PSM table + visualization frames in session state
    st.session_state[get_state_key("psm_frame")] = psm_frame
    st.session_state[get_state_key("psm_lookup")] = psm_lookup
    st.session_state[get_state_key("peaks_frame")] = frames["peaks"]
    st.session_state[get_state_key("sequence_frame")] = frames["sequence"]
    st.session_state[get_state_key("heatmap_frame")] = frames["heatmap"]

    # Stash annotation payloads + comparison rows for the Explore page
    st.session_state[get_state_key("annotated_queries")] = {
        psm_id: (spec, df) for psm_id, spec, df in annotated_queries
    }
    st.session_state[get_state_key("reference_store")] = reference_store
    st.session_state[get_state_key("comparison_frame")] = (
        pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()
    )

    # Mark data ready, bump cache version, and persist to disk
    st.session_state[get_state_key("data_ready")] = True
    st.session_state[get_state_key("cache_version")] = \
        st.session_state.get(get_state_key("cache_version"), 0) + 1

    save_analysis_cache()


def _stream_and_annotate(
    params: dict, variant_groups_copy, mzmls,
    psm_lookup, annotation_kwargs, koina_kwargs,
    total_psms: int, tracker,
) -> tuple[list, dict, list]:
    """Run the streaming-annotation orchestrator with Streamlit progress + error wiring.

    Returns ``(annotated_queries, reference_store, comparison_rows)``.
    """
    def on_progress(loaded_count: int) -> None:
        if loaded_count % 10 == 0 or loaded_count == total_psms:
            tracker.update(loaded_count, total_psms, units="spectra")

    def on_file_start(path, is_reference):
        kind = "reference" if is_reference else "query"
        tracker.info(f"Processing {kind}: {Path(path).name}")

    return run_streaming_annotation(
        variant_groups=variant_groups_copy,
        mzml_paths=mzmls,
        reference_substring=get_param(params, "reference_substring", ""),
        psm_lookup=psm_lookup,
        annotation_kwargs=annotation_kwargs,
        koina_kwargs=koina_kwargs,
        workers=get_param(params, "workers", DEFAULTS.workers),
        on_progress=on_progress,
        on_file_start=on_file_start,
        on_file_complete=lambda path, loaded, total: logger.info(
            f"Completed {Path(path).name}: {loaded}/{total}"
        ),
        on_error=lambda error: logger.error(f"Streaming error: {error}"),
    )


def _validate_load_inputs(params: dict) -> tuple | None:
    """Pre-flight: ensure search ran and mzMLs were picked. Returns ``(variant_groups, mzmls)``."""
    variant_groups = st.session_state.get(get_state_key("variant_groups"))
    if not variant_groups:
        st.error("No search results. Run peptide search first.")
        return None
    mzmls = get_param(params, "selected_mzmls", [])
    if not mzmls:
        st.error("No mzML files selected.")
        return None
    return variant_groups, mzmls


def _finalize_load(
    psm_frame, psm_lookup, annotated_queries, reference_store, comparison_rows, tracker,
) -> None:
    """Persist results + update the status block with the run summary."""
    tracker.info(f"Building DataFrames from {len(annotated_queries)} spectra...")
    _store_analysis_results(psm_frame, psm_lookup, annotated_queries,
                            reference_store, comparison_rows)
    tracker.finalize(
        f"Loaded {len(annotated_queries)} spectra, "
        f"{len(comparison_rows)} comparisons in {format_time(tracker.elapsed)}"
    )


def run_search_and_load(params: dict) -> None:
    """Stream spectra, annotate, compute metrics, build DataFrames."""
    # Reset any previous analysis cache + validate prerequisites
    clear_analysis_cache()

    inputs = _validate_load_inputs(params)
    if inputs is None:
        return
    variant_groups, mzmls = inputs

    # Prepare inputs: snapshot variants, build PSM table, gather kwargs
    variant_groups_copy = copy.deepcopy(variant_groups)
    psm_frame, psm_lookup = build_psm_table(variant_groups_copy)
    getter = _params_getter(params)
    annotation_kwargs = build_annotation_kwargs(getter)
    koina_kwargs = build_koina_kwargs(getter)

    total_psms = sum(len(plist) for plist in variant_groups_copy.values())

    # Stream + annotate inside a progress status; persist results when non-empty
    with progress_status(f"Loading {total_psms} spectra...") as tracker:
        annotated_queries, reference_store, comparison_rows = _stream_and_annotate(
            params, variant_groups_copy, mzmls,
            psm_lookup, annotation_kwargs, koina_kwargs,
            total_psms, tracker,
        )
        if not annotated_queries:
            tracker.finalize("No spectra could be loaded or annotated", state="error")
            return
        _finalize_load(
            psm_frame, psm_lookup, annotated_queries, reference_store, comparison_rows, tracker,
        )

    st.rerun()
