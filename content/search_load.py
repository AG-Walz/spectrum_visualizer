"""
Search & Load Page — Peptide search, spectrum streaming, annotation, metrics.

Long-running operations use st.status with progress bars.
The results display section is an @st.fragment.
"""
import copy
import logging
import time as _time
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from psm_utils import Peptidoform as _PF

from src.common.common import page_setup, save_params
from src.koina_prediction import predict_koina_spectrum, KOINA_AVAILABLE
from src.parallel_processor import compute_ion_match_statistics
from src.parquet_builder import (
    build_psm_table, build_peaks_dataframe,
    build_sequence_dataframe, build_heatmap_dataframe,
)
from src.peak_annotation import prepare_spectrum_for_plot, prepare_koina_spectrum_for_plot
from src.peptide_search import search_peptides_in_idxmls, summarize_search_results
from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key, LoadedSpectrum
from src.spectrum_metrics import compute_all_metrics, format_metrics_text, METRIC_SHORT_TO_DISPLAY
from content._shared import (
    STATE_PREFIX, get_state_key,
    restore_analysis_if_needed, save_analysis_cache, clear_analysis_cache,
)
from utils.helpers import format_time, parse_peptide_list

params = page_setup()
logger = logging.getLogger(__name__)


def _p(key, default=None):
    """Get a param value: session_state first (if widget present), then params.json."""
    return st.session_state.get(key, params.get(key, default))


# ---------------------------------------------------------------------------
# Peptide search
# ---------------------------------------------------------------------------

def run_peptide_search(idxmls: list, peptides: list):
    """Execute peptide search."""
    mzmls = _p('selected-mzmls', [])
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

        try:
            workers = _p('workers', 4)
            results = search_peptides_in_idxmls(
                idxmls, peptides, mzml_paths=mzmls,
                progress_callback=progress_callback, error_callback=error_callback,
                workers=workers,
            )
            st.session_state[get_state_key('variant_groups')] = results
            progress_bar.progress(1.0, "Search complete!")

            if results:
                total_matches = sum(len(vg.matches) for vg in results.values())
                status.update(label=f"Found {total_matches} PSMs across {len(results)} variants",
                              state="complete")

                save_params(params)

                unresolved = sum(1 for vg in results.values() for m in vg.matches if not m.mzml_filename)
                if unresolved > 0:
                    st.warning(f"{unresolved} PSMs could not be mapped to mzML files")
                if error_messages:
                    with st.expander(f"{len(error_messages)} warnings", expanded=False):
                        for msg in error_messages[:10]:
                            st.caption(msg)
            else:
                status.update(label="No matches found", state="error")
        except Exception as e:
            status.update(label=f"Search failed: {e}", state="error")
            logger.exception("Peptide search error")


# ---------------------------------------------------------------------------
# Spectrum loading + annotation + metrics
# ---------------------------------------------------------------------------

def _get_koina_prediction(peptidoform: str, charge: int, annotation_kwargs: dict):
    """Generate Koina AI-predicted spectrum as a fallback reference."""
    try:
        if not KOINA_AVAILABLE:
            return None

        model_id = _p("koina-model", "ms2pip_timsTOF2024")
        ce = _p("koina-ce", 35)
        instrument = _p("koina-instrument", "TIMSTOF")

        pred_df = predict_koina_spectrum(
            peptidoform_str=peptidoform,
            charge=charge,
            model_id=model_id,
            collision_energy=ce,
            instrument_type=instrument,
        )

        if pred_df is None or pred_df.empty:
            return None

        try:
            base_seq = _PF(peptidoform.split("/")[0]).sequence
        except Exception:
            base_seq = peptidoform.split("/")[0]

        ref_spec = LoadedSpectrum(
            spectrum_id="koina_prediction",
            mz=pred_df["mz"].values.astype(np.float64),
            intensity=pred_df["intensity"].values.astype(np.float64),
            precursor_mz=0.0,
            precursor_charge=charge,
            peptidoform=peptidoform,
            base_sequence=base_seq,
            rt=0.0,
            source_file="Koina",
            is_reference=True,
        )

        # Use Koina's native AI-predicted annotations directly instead of
        # re-matching against PyOpenMS theoretical fragments
        ref_df = prepare_koina_spectrum_for_plot(pred_df, **annotation_kwargs)

        return (ref_spec, ref_df)

    except Exception as e:
        logger.warning(f"Koina prediction failed for {peptidoform}: {e}")
        return None


def run_search_and_load():
    """Stream spectra, capture references, annotate, compute metrics, build DataFrames."""
    clear_analysis_cache()

    variant_groups = st.session_state.get(get_state_key('variant_groups'))
    if not variant_groups:
        st.error("No search results. Run peptide search first.")
        return

    mzmls = _p('selected-mzmls', [])
    if not mzmls:
        st.error("No mzML files selected.")
        return

    ref_substring = _p('reference-substring', '')
    ion_types = _p('ion-types', ['b', 'y'])
    mz_tolerance = _p('mz-tolerance', 0.02)
    normalize = _p('intensity-mode', 'Relative') == 'Relative'
    workers = _p('workers', 4)

    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)

    annotated_queries = []
    reference_store = {}
    comparison_rows = []

    total_psms = sum(len(vg.matches) for vg in vg_copy.values())

    add_losses = _p('add-losses', False)
    add_first_prefix_ion = _p('add-first-prefix-ion', True)
    remove_precursor_peaks = _p('remove-precursor-peaks', False)

    annotation_kwargs = dict(
        ion_types=ion_types,
        mz_tolerance=mz_tolerance,
        normalize=normalize,
        add_losses=add_losses,
        add_first_prefix_ion=add_first_prefix_ion,
        remove_precursor_peaks=remove_precursor_peaks,
    )

    with st.status(f"Loading {total_psms} spectra...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        _load_start = _time.time()
        loaded_count = [0]

        def on_reference_loaded(variant_key, count):
            pass

        def on_query_ready(result):
            key = spectrum_lookup_key(result.query_spectrum)
            psm_id = psm_lookup.get(key)
            if psm_id is None:
                return

            try:
                query_df = prepare_spectrum_for_plot(
                    result.query_spectrum.mz,
                    result.query_spectrum.intensity,
                    result.query_spectrum.peptidoform,
                    result.query_spectrum.precursor_charge,
                    **annotation_kwargs,
                )
                annotated_queries.append((psm_id, result.query_spectrum, query_df))
            except Exception as e:
                logger.warning(f"Annotation failed psm_id={psm_id}: {e}")
                return

            vk = result.variant_key
            if vk not in reference_store and result.reference_spectra:
                ref_list = []
                for ref_spec in result.reference_spectra:
                    try:
                        ref_df = prepare_spectrum_for_plot(
                            ref_spec.mz, ref_spec.intensity,
                            ref_spec.peptidoform, ref_spec.precursor_charge,
                            **annotation_kwargs,
                        )
                        ref_list.append((ref_spec, ref_df))
                    except Exception as e:
                        logger.warning(f"Reference annotation failed: {e}")
                reference_store[vk] = ref_list

            refs = reference_store.get(vk, [])
            if not refs and _p('koina-enabled', False):
                koina_ref = _get_koina_prediction(
                    result.query_spectrum.peptidoform,
                    result.query_spectrum.precursor_charge,
                    annotation_kwargs,
                )
                if koina_ref is not None:
                    reference_store[vk] = [koina_ref]
                    refs = reference_store[vk]

            for ref_spec, ref_df in refs:
                try:
                    metrics = compute_all_metrics(query_df, ref_df)
                    row = {
                        "psm_id": psm_id,
                        "Peptide": result.base_sequence,
                        "Modification": result.modification_key,
                        "Charge": result.query_spectrum.precursor_charge,
                        "Query": Path(result.query_spectrum.source_file).stem,
                        "Reference": Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina",
                        "query_spectrum_id": result.query_spectrum.spectrum_id,
                        "ref_spectrum_id": ref_spec.spectrum_id,
                        "query_rt": result.query_spectrum.rt,
                        "ref_rt": ref_spec.rt,
                    }
                    for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
                        row[display_name] = metrics.get(short_key, float("nan"))
                    for count_key in ("total_q", "total_r", "ann_q", "ann_r", "shared"):
                        row[count_key] = metrics.get(count_key, 0)

                    try:
                        ims_stats = compute_ion_match_statistics(
                            query_df, result.query_spectrum.peptidoform,
                            result.query_spectrum.precursor_charge,
                            result.query_spectrum.precursor_mz,
                        )
                        row.update(ims_stats)
                    except Exception:
                        pass

                    comparison_rows.append(row)
                except Exception as e:
                    logger.warning(f"Metric computation failed: {e}")

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

        def on_file_complete(path, loaded, total):
            logger.info(f"Completed {Path(path).name}: {loaded}/{total} spectra")

        def error_callback(error):
            logger.error(f"Streaming error: {error}")

        load_spectra_streaming(
            variant_groups=vg_copy,
            mzml_paths=[Path(m) for m in mzmls],
            reference_substring=ref_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=on_file_complete,
            error_callback=error_callback,
            workers=workers,
        )

        progress_bar.progress(1.0)
        elapsed = _time.time() - _load_start

        if not annotated_queries:
            status.update(label="No spectra could be loaded or annotated", state="error")
            return

        status_text.caption(f"Building DataFrames from {len(annotated_queries)} spectra...")

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

        if comparison_rows:
            comparison_df = pd.DataFrame(comparison_rows)
        else:
            comparison_df = pd.DataFrame()
        st.session_state[get_state_key('comparison_df')] = comparison_df

        st.session_state[get_state_key('data_ready')] = True
        st.session_state[get_state_key('cache_version')] = \
            st.session_state.get(get_state_key('cache_version'), 0) + 1

        save_analysis_cache()

        n_comparisons = len(comparison_rows)
        status.update(
            label=f"Loaded {len(annotated_queries)} spectra, {n_comparisons} comparisons "
                  f"in {format_time(elapsed)}",
            state="complete",
        )

    st.rerun()


# ---------------------------------------------------------------------------
# Search results display
# ---------------------------------------------------------------------------

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

        total_matches = sum(len(vg.matches) for vg in variant_groups.values())
        st.success(f"Found {total_matches} PSMs across {len(variant_groups)} peptide variants. "
                   f"Click **Search & Load Spectra** below to load and annotate them.")
        mzmls = _p('selected-mzmls', [])
        if mzmls:
            st.subheader("Load Spectra & Compute Metrics")
            ref_sub = _p('reference-substring', '')
            if ref_sub:
                st.caption(f"Will pair query spectra with references (substring: '{ref_sub}') "
                           f"and compute 7 comparison metrics per pair.")
            else:
                koina_on = _p('koina-enabled', False)
                if koina_on:
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


# ---------------------------------------------------------------------------
# Page render
# ---------------------------------------------------------------------------

def render():
    restore_analysis_if_needed()
    st.title("Search & Load")
    st.caption("Search for peptides in idXML files and load spectra for analysis")

    idxmls = _p('selected-idxmls', [])
    peptide_text = _p('target-peptides', '')
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
