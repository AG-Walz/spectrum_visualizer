"""
Search & Load Page — Peptide search, spectrum streaming, annotation, metrics.

Long-running operations use st.status with progress bars.
The results display section is an @st.fragment.
"""
import copy
import logging
import time as _time
from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_app.utils.state import page_setup, save_params
from mspecviz.insight_dataframes import (
    build_psm_table, build_peaks_dataframe,
    build_sequence_dataframe, build_heatmap_dataframe,
)
from mspecviz.peptide_search import search_peptides_in_idxmls, summarize_search_results
from mspecviz.spectrum_loader import load_spectra_streaming
from mspecviz.query_processing import process_query_result
from streamlit_app.views.page_state import (
    STATE_PREFIX, get_state_key, get_param,
    restore_analysis_if_needed, save_analysis_cache, clear_analysis_cache,
)
from mspecviz.common.text_utils import format_time, parse_peptide_list

params = page_setup()
logger = logging.getLogger(__name__)


def run_peptide_search(idxmls: list, peptides: list):
    """Execute peptide search."""
    mzmls = get_param(params, 'selected-mzmls', [])
    error_messages = []

    with st.status("Searching peptides...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        _search_start = _time.time()

        def progress_callback(current, total, message):
            if total > 0:
                frac = current / total
                progress_bar.progress(frac, message)
                elapsed = _time.time() - _search_start
                if frac > 0 and current < total:
                    remaining = (elapsed / frac) - elapsed
                    status_text.caption(f"{current}/{total} files  "
                                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)")
                elif current >= total:
                    status_text.caption(f"Done in {format_time(elapsed)}")

        def error_callback(error_msg: str):
            error_messages.append(error_msg)

        workers = get_param(params, 'workers', 4)
        results = search_peptides_in_idxmls(
            idxmls, peptides, mzml_paths=mzmls,
            progress_callback=progress_callback, error_callback=error_callback,
            workers=workers,
        )
        st.session_state[get_state_key('variant_groups')] = results
        progress_bar.progress(1.0, "Search complete!")

        if results:
            total_matches = sum(len(variant_group.matches) for variant_group in results.values())
            status.update(label=f"Found {total_matches} PSMs across {len(results)} variants",
                          state="complete")

            save_params(params)

            unresolved = sum(1 for variant_group in results.values() for m in variant_group.matches if not m.mzml_filename)
            if unresolved > 0:
                st.warning(f"{unresolved} PSMs could not be mapped to mzML files")
            if error_messages:
                with st.expander(f"{len(error_messages)} warnings", expanded=False):
                    for msg in error_messages[:10]:
                        st.caption(msg)
        else:
            status.update(label="No matches found", state="error")


def _build_annotation_kwargs(params) -> dict:
    """Collect annotation settings from params into a kwargs dict."""
    return dict(
        ion_types=get_param(params, 'ion-types', ['b', 'y']),
        mz_tolerance=get_param(params, 'mz-tolerance', 0.02),
        normalize=get_param(params, 'intensity-mode', 'Relative') == 'Relative',
        add_losses=get_param(params, 'add-losses', False),
        add_first_prefix_ion=get_param(params, 'add-first-prefix-ion', True),
        remove_precursor_peaks=get_param(params, 'remove-precursor-peaks', False),
    )


def _store_analysis_results(psm_lf, psm_lookup, annotated_queries,
                            reference_store, comparison_rows):
    """Persist all analysis results to session state and disk cache."""
    peaks_lf = build_peaks_dataframe([(pid, df) for pid, _, df in annotated_queries])
    seq_lf = build_sequence_dataframe([(pid, s) for pid, s, _ in annotated_queries])
    heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated_queries])

    st.session_state[get_state_key('psm_lf')] = psm_lf
    st.session_state[get_state_key('psm_lookup')] = psm_lookup
    st.session_state[get_state_key('peaks_lf')] = peaks_lf
    st.session_state[get_state_key('seq_lf')] = seq_lf
    st.session_state[get_state_key('heatmap_lf')] = heatmap_lf

    st.session_state[get_state_key('annotated_queries')] = {
        pid: (spec, df) for pid, spec, df in annotated_queries
    }
    st.session_state[get_state_key('reference_store')] = reference_store
    st.session_state[get_state_key('comparison_df')] = (
        pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()
    )

    st.session_state[get_state_key('data_ready')] = True
    st.session_state[get_state_key('cache_version')] = \
        st.session_state.get(get_state_key('cache_version'), 0) + 1

    save_analysis_cache()


def run_search_and_load():
    """Stream spectra, annotate, compute metrics, build DataFrames."""
    clear_analysis_cache()

    variant_groups = st.session_state.get(get_state_key('variant_groups'))
    if not variant_groups:
        st.error("No search results. Run peptide search first.")
        return

    mzmls = get_param(params, 'selected-mzmls', [])
    if not mzmls:
        st.error("No mzML files selected.")
        return

    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)
    annotation_kwargs = _build_annotation_kwargs(params)

    annotated_queries = []
    reference_store = {}
    comparison_rows = []
    total_psms = sum(len(vg.matches) for vg in vg_copy.values())

    with st.status(f"Loading {total_psms} spectra...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        _load_start = _time.time()
        loaded_count = [0]

        def on_query_ready(result):
            processed = process_query_result(
                result, psm_lookup, annotation_kwargs, reference_store,
                koina_enabled=get_param(params, 'koina-enabled', False),
                koina_model=get_param(params, "koina-model", "ms2pip_timsTOF2024"),
                koina_ce=get_param(params, "koina-ce", 35),
                koina_instrument=get_param(params, "koina-instrument", "TIMSTOF"),
            )
            if processed is None:
                return
            psm_id, query_spec, query_df, comp_rows = processed
            annotated_queries.append((psm_id, query_spec, query_df))
            comparison_rows.extend(comp_rows)

            loaded_count[0] += 1
            if loaded_count[0] % 10 == 0 or loaded_count[0] == total_psms:
                frac = min(loaded_count[0] / max(total_psms, 1), 1.0)
                elapsed = _time.time() - _load_start
                progress_bar.progress(frac)
                status_text.caption(f"{loaded_count[0]}/{total_psms} spectra "
                                    f"({format_time(elapsed)} elapsed)")

        def on_file_start(path, is_reference):
            kind = "reference" if is_reference else "query"
            status_text.caption(f"Processing {kind}: {Path(path).name}")

        load_spectra_streaming(
            variant_groups=vg_copy,
            mzml_paths=[Path(m) for m in mzmls],
            reference_substring=get_param(params, 'reference-substring', ''),
            on_reference_loaded=lambda vk, c: None,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=lambda p, l, t: logger.info(f"Completed {Path(p).name}: {l}/{t}"),
            error_callback=lambda e: logger.error(f"Streaming error: {e}"),
            workers=get_param(params, 'workers', 4),
        )

        progress_bar.progress(1.0)
        elapsed = _time.time() - _load_start

        if not annotated_queries:
            status.update(label="No spectra could be loaded or annotated", state="error")
            return

        status_text.caption(f"Building DataFrames from {len(annotated_queries)} spectra...")
        _store_analysis_results(psm_lf, psm_lookup, annotated_queries,
                                reference_store, comparison_rows)

        status.update(
            label=f"Loaded {len(annotated_queries)} spectra, {len(comparison_rows)} comparisons "
                  f"in {format_time(elapsed)}",
            state="complete",
        )

    st.rerun()


def display_search_results(results: dict):
    """Display search results table."""
    summary_lf = summarize_search_results(results)
    summary_df = summary_lf.collect().to_pandas()
    if not summary_df.empty:
        display_df = summary_df.rename(columns={
            'base_sequence': 'Peptide', 'modification': 'Modification',
            'psm_count': 'PSMs', 'idxml_files': 'idXML Files',
            'mzml_resolved': 'mzML Resolved', 'unresolved': 'Unresolved'
        })
        display_df = display_df.drop(columns=['display_label'], errors='ignore')
        st.dataframe(display_df, width="stretch", hide_index=True)


@st.fragment
def search_results_fragment():
    """Display search results and load controls as a fragment."""
    variant_groups = st.session_state.get(get_state_key('variant_groups'))
    if variant_groups:
        st.subheader("Search Results")
        display_search_results(variant_groups)

        total_matches = sum(len(variant_group.matches) for variant_group in variant_groups.values())
        st.success(f"Found {total_matches} PSMs across {len(variant_groups)} peptide variants. "
                   f"Click **Search & Load Spectra** below to load and annotate them.")
        mzmls = get_param(params, 'selected-mzmls', [])
        if mzmls:
            st.subheader("Load Spectra & Compute Metrics")
            reference_substring = get_param(params, 'reference-substring', '')
            if reference_substring:
                st.caption(f"Will pair query spectra with references (substring: '{reference_substring}') "
                           f"and compute 7 comparison metrics per pair.")
            else:
                koina_enabled = get_param(params, 'koina-enabled', False)
                if koina_enabled:
                    st.caption("No reference configured -- will use Koina AI prediction as fallback.")
                else:
                    st.caption("No reference configured -- will annotate individual spectra.")
            if st.button("Search & Load Spectra", key=f"{STATE_PREFIX}_load_btn", type="primary"):
                run_search_and_load()

            if st.session_state.get(get_state_key('data_ready')):
                comparison_df = st.session_state.get(get_state_key('comparison_df'))
                if comparison_df is not None and not comparison_df.empty:
                    st.success(f"Data loaded with {len(comparison_df)} query/reference comparisons. "
                               f"Go to **Explore** page.")
                else:
                    st.success("Data loaded. Go to **Explore** page.")
        else:
            st.warning("Please select mzML files first (Setup page)")
    else:
        st.caption("Search results will appear here after running the search.")


def render():
    restore_analysis_if_needed()
    st.title("Search & Load")
    st.caption("Search for peptides in idXML files and load spectra for analysis")

    idxmls = get_param(params, 'selected-idxmls', [])
    peptide_text = get_param(params, 'target-peptides', '')
    peptides = parse_peptide_list(peptide_text)

    if not idxmls:
        st.warning("Please select idXML files first (Setup page)")
        return
    if not peptides:
        st.warning("Please enter target peptides first (Setup page)")
        return

    if st.button("Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(idxmls, peptides)

    search_results_fragment()

    save_params(params)


render()
