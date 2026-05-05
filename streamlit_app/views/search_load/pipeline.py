"""
Search & Load pipeline helpers — peptide search, spectrum streaming, annotation.

These are the long-running operations triggered from ``page.py`` buttons.
They take ``params`` explicitly (not from module state) so the split remains
testable and re-runnable without reloading the view module.
"""
from __future__ import annotations

import copy
import logging
import time as _time
from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_app.utils.state import get_param, get_state_key, save_params
from streamlit_app.utils.page_state import clear_analysis_cache, save_analysis_cache
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.text_utils import format_time
from mspecviz.dataframes import (
    build_heatmap_dataframe,
    build_peaks_dataframe,
    build_psm_table,
    build_sequence_dataframe,
)
from mspecviz.pipeline import process_query_result
from mspecviz.search import search_peptides_in_idxmls
from mspecviz.spectra import load_spectra_streaming

logger = logging.getLogger(__name__)


def _make_search_progress_callback(start_time: float, progress_bar, status_text):
    """Build the per-file progress callback for the peptide search stage."""
    def progress_callback(current, total, message):
        if total <= 0:
            return
        fraction = current / total
        progress_bar.progress(fraction, message)
        elapsed = _time.time() - start_time
        if fraction > 0 and current < total:
            remaining = (elapsed / fraction) - elapsed
            status_text.caption(
                f"{current}/{total} files  "
                f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)"
            )
        elif current >= total:
            status_text.caption(f"Done in {format_time(elapsed)}")
    return progress_callback


def _show_search_results_summary(results, error_messages, status, params: dict) -> None:
    """Update the status block with PSM/variant counts and any warnings."""
    if not results:
        status.update(label="No matches found", state="error")
        return

    total_matches = sum(len(variant_group.matches) for variant_group in results.values())
    status.update(
        label=f"Found {total_matches} PSMs across {len(results)} variants",
        state="complete",
    )
    save_params(params)

    unresolved = sum(
        1 for variant_group in results.values()
        for m in variant_group.matches if not m.mzml_filename
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

    with st.status("Searching peptides...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        progress_callback = _make_search_progress_callback(
            _time.time(), progress_bar, status_text,
        )

        results = search_peptides_in_idxmls(
            idxmls, peptides, mzml_paths=mzmls,
            progress_callback=progress_callback,
            error_callback=error_messages.append,
            workers=get_param(params, "workers", DEFAULTS.workers),
        )
        st.session_state[get_state_key("variant_groups")] = results
        progress_bar.progress(1.0, "Search complete!")

        _show_search_results_summary(results, error_messages, status, params)


def _build_annotation_kwargs(params: dict) -> dict:
    """Collect annotation settings from params into a kwargs dict."""
    return dict(
        ion_types=get_param(params, "ion_types", DEFAULTS.ion_types),
        mz_tolerance=get_param(params, "mz_tolerance", DEFAULTS.mz_tolerance),
        normalize=get_param(params, "intensity_mode", DEFAULTS.intensity_mode) == "Relative",
        add_neutral_losses=get_param(params, "add_neutral_losses", DEFAULTS.add_neutral_losses),
    )


def _store_analysis_results(psm_frame, psm_lookup, annotated_queries,
                            reference_store, comparison_rows) -> None:
    """Persist all analysis results to session state and disk cache."""
    peaks_frame = build_peaks_dataframe([(psm_id, df) for psm_id, _, df in annotated_queries])
    sequence_frame = build_sequence_dataframe([(psm_id, s) for psm_id, s, _ in annotated_queries])
    heatmap_frame = build_heatmap_dataframe([s for _, s, _ in annotated_queries])

    st.session_state[get_state_key("psm_frame")] = psm_frame
    st.session_state[get_state_key("psm_lookup")] = psm_lookup
    st.session_state[get_state_key("peaks_frame")] = peaks_frame
    st.session_state[get_state_key("sequence_frame")] = sequence_frame
    st.session_state[get_state_key("heatmap_frame")] = heatmap_frame

    st.session_state[get_state_key("annotated_queries")] = {
        psm_id: (spec, df) for psm_id, spec, df in annotated_queries
    }
    st.session_state[get_state_key("reference_store")] = reference_store
    st.session_state[get_state_key("comparison_frame")] = (
        pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()
    )

    st.session_state[get_state_key("data_ready")] = True
    st.session_state[get_state_key("cache_version")] = \
        st.session_state.get(get_state_key("cache_version"), 0) + 1

    save_analysis_cache()


def _koina_kwargs(params: dict) -> dict:
    """Collect Koina settings from params for ``process_query_result``."""
    return dict(
        koina_enabled=get_param(params, "koina_enabled", DEFAULTS.koina_enabled),
        koina_model=get_param(params, "koina_model", DEFAULTS.koina_model),
        koina_collision_energy=get_param(
            params, "koina_collision_energy", DEFAULTS.koina_collision_energy,
        ),
        koina_instrument_type=get_param(
            params, "koina_instrument_type", DEFAULTS.koina_instrument_type,
        ),
    )


def _make_query_ready_callback(
    psm_lookup, annotation_kwargs, reference_store, koina_kwargs,
    annotated_queries: list, comparison_rows: list,
    total_psms: int, start_time: float, progress_bar, status_text,
):
    """Return ``on_query_ready`` callback that buffers results and updates UI."""
    loaded_count = [0]

    def on_query_ready(result):
        processed = process_query_result(
            result, psm_lookup, annotation_kwargs, reference_store, **koina_kwargs,
        )
        if processed is None:
            return
        psm_id, query_spectrum, query_frame, comparison_rows_for_psm = processed
        annotated_queries.append((psm_id, query_spectrum, query_frame))
        comparison_rows.extend(comparison_rows_for_psm)

        loaded_count[0] += 1
        if loaded_count[0] % 10 == 0 or loaded_count[0] == total_psms:
            fraction = min(loaded_count[0] / max(total_psms, 1), 1.0)
            elapsed = _time.time() - start_time
            progress_bar.progress(fraction)
            status_text.caption(
                f"{loaded_count[0]}/{total_psms} spectra "
                f"({format_time(elapsed)} elapsed)"
            )

    return on_query_ready


def _stream_and_annotate(
    params: dict, variant_groups_copy, mzmls,
    psm_lookup, annotation_kwargs,
    annotated_queries: list, reference_store: dict, comparison_rows: list,
    total_psms: int, status_text, progress_bar,
) -> None:
    """Run ``load_spectra_streaming`` with progress + error wiring."""
    start_time = _time.time()
    on_query_ready = _make_query_ready_callback(
        psm_lookup, annotation_kwargs, reference_store, _koina_kwargs(params),
        annotated_queries, comparison_rows,
        total_psms, start_time, progress_bar, status_text,
    )

    def on_file_start(path, is_reference):
        kind = "reference" if is_reference else "query"
        status_text.caption(f"Processing {kind}: {Path(path).name}")

    load_spectra_streaming(
        variant_groups=variant_groups_copy,
        mzml_paths=[Path(m) for m in mzmls],
        reference_substring=get_param(params, "reference_substring", ""),
        on_reference_loaded=lambda vk, c: None,
        on_query_ready=on_query_ready,
        on_file_start=on_file_start,
        on_file_complete=lambda p, loaded, total: logger.info(
            f"Completed {Path(p).name}: {loaded}/{total}"
        ),
        error_callback=lambda e: logger.error(f"Streaming error: {e}"),
        workers=get_param(params, "workers", DEFAULTS.workers),
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
    psm_frame, psm_lookup, annotated_queries, reference_store, comparison_rows,
    elapsed: float, status, status_text,
) -> None:
    """Persist results + update the status block with the run summary."""
    status_text.caption(f"Building DataFrames from {len(annotated_queries)} spectra...")
    _store_analysis_results(psm_frame, psm_lookup, annotated_queries,
                            reference_store, comparison_rows)
    status.update(
        label=f"Loaded {len(annotated_queries)} spectra, "
              f"{len(comparison_rows)} comparisons in {format_time(elapsed)}",
        state="complete",
    )


def run_search_and_load(params: dict) -> None:
    """Stream spectra, annotate, compute metrics, build DataFrames."""
    clear_analysis_cache()

    inputs = _validate_load_inputs(params)
    if inputs is None:
        return
    variant_groups, mzmls = inputs

    variant_groups_copy = copy.deepcopy(variant_groups)
    psm_frame, psm_lookup = build_psm_table(variant_groups_copy)
    annotation_kwargs = _build_annotation_kwargs(params)

    annotated_queries: list = []
    reference_store: dict = {}
    comparison_rows: list = []
    total_psms = sum(len(group.matches) for group in variant_groups_copy.values())

    with st.status(f"Loading {total_psms} spectra...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        load_start = _time.time()

        _stream_and_annotate(
            params, variant_groups_copy, mzmls,
            psm_lookup, annotation_kwargs,
            annotated_queries, reference_store, comparison_rows,
            total_psms, status_text, progress_bar,
        )
        progress_bar.progress(1.0)

        if not annotated_queries:
            status.update(label="No spectra could be loaded or annotated", state="error")
            return

        _finalize_load(
            psm_frame, psm_lookup, annotated_queries, reference_store, comparison_rows,
            _time.time() - load_start, status, status_text,
        )

    st.rerun()
