"""Per-query workflow layer: annotate one query spectrum, resolve its references,
score every (query, reference) pair, and package the results either as comparison
rows or as ComparisonTask plot jobs.

`process_query_result` is the per-query step of the streaming pipeline; it runs
through the annotation phases (annotate query, cache experimental references,
optional Koina AI fallback) and the scoring via core.comparison_metrics.
`run_streaming_annotation` is the top-level load-stage driver: it streams every
mzML and feeds each arriving query to `process_query_result`, buffering the
outcomes for downstream consumption.
`build_plot_tasks` converts already-annotated queries into one plot task per pair,
pre-computing the metrics overlay so the parallel plot workers never recompute it.

run_streaming_annotation called by cli.main and streamlit_app.views.load_annotate;
build_plot_tasks called by core.export.run_export;
build_psm_to_peptidoform_mapping and format_spectrum_labels also used by
streamlit_app.views.explore.mirror_plots.
"""
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd
from psm_utils import Peptidoform, PSM

from core.comparison_metrics import (
    METRIC_KEY_TO_COLUMN_NAME,
    compute_all_metrics,
    format_metrics_text,
    rt_delta,
)
from core.peak_annotation import (
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.defaults import DEFAULTS
from core.plotting.mirror import group_mz_axis_max
from core.plotting.plot_worker import ComparisonTask
from core.spectra import load_spectra_streaming
from core.spectra.koina_prediction import predict_koina_spectrum
from core.spectra.streaming_types import psm_lookup_key, source_label


# ============================================================================
# Comparison-summary columns — single source of truth for the CSV/XLSX export
# ============================================================================
# The comparison summary (one row per query/reference pair) is assembled as plain
# dicts in `process_query_result`. These constants name and order its columns in one
# place instead of leaving the layout implied by dict-insertion order; the ordering
# is identity -> metrics (from METRIC_KEY_TO_COLUMN_NAME) -> counts. `build_comparison_frame`
# reindexes to this order, so adding a column here is what makes it appear in the export.

COMPARISON_IDENTITY_COLUMNS = (
    "psm_id", "Peptide", "Modification", "Charge", "Query", "Reference",
    "query_spectrum_id", "ref_spectrum_id", "query_rt", "ref_rt", "rt_delta",
)
# Retention-time columns surfaced in the Explore summaries (#113). Kept apart from the
# rest of the identity block because the summary table and the summary statistics both
# format them as seconds, unlike the string/int identity columns beside them.
COMPARISON_RT_COLUMNS = ("query_rt", "ref_rt", "rt_delta")
COMPARISON_COUNT_COLUMNS = (
    "query_peak_count", "reference_peak_count",
    "query_annotation_count", "reference_annotation_count", "shared",
)
COMPARISON_SUMMARY_COLUMNS = (
    *COMPARISON_IDENTITY_COLUMNS,
    *METRIC_KEY_TO_COLUMN_NAME.values(),
    *COMPARISON_COUNT_COLUMNS,
)


def build_comparison_frame(comparison_rows: list[dict]) -> pd.DataFrame:
    """Package comparison rows into a DataFrame with the canonical column order.
    Empty input yields an empty, column-less frame — the "no metrics" placeholder
    both front ends and the analysis cache already expect.
    """
    if not comparison_rows:
        return pd.DataFrame()
    return pd.DataFrame(comparison_rows, columns=list(COMPARISON_SUMMARY_COLUMNS))


# ============================================================================
# Per-query annotation phases (annotate query, cache references, Koina fallback)
# ============================================================================

def get_koina_prediction(
    peptidoform, charge, annotation_kwargs,
    koina_model: str = DEFAULTS.koina_model,
    koina_collision_energy: int = DEFAULTS.koina_collision_energy,
    koina_instrument_type: str = DEFAULTS.koina_instrument_type,
):
    """Predict a spectrum with Koina and wrap it as an annotated reference;
    return (spectrum, frame) or None.
    """
    # Back the Koina reference with a synthetic PSM carrier built up front, so predicted
    # references are PSMs like experimental ones (no idXML provenance, run="Koina:<model>").
    # Charge is embedded in the peptidoform so get_precursor_charge() resolves it;
    # predict_koina_spectrum reads the peptidoform + charge straight off this PSM.
    reference_psm = PSM(
        peptidoform=Peptidoform(f"{peptidoform}/{charge}"),
        spectrum_id="koina_prediction",
        run=f"Koina:{koina_model}",
        retention_time=0.0,
    )

    # Predict: call Koina with model + CE + instrument, bail on empty frame
    predicted_frame = predict_koina_spectrum(
        reference_psm,
        model_id=koina_model,
        collision_energy=koina_collision_energy,
        instrument_type=koina_instrument_type,
    )
    if predicted_frame is None or predicted_frame.empty:
        return None

    # Attach the predicted peaks and the theoretical precursor m/z to the carrier
    # (no measured precursor exists for a predicted spectrum).
    reference_psm.spectrum = {
        "mz": predicted_frame["mz"].to_numpy(dtype=np.float64),
        "intensity": predicted_frame["intensity"].to_numpy(dtype=np.float64),
    }
    reference_psm.precursor_mz = reference_psm.peptidoform.theoretical_mz

    # Annotate: reuse Koina's native AI-predicted ion labels directly
    predicted_reference_frame = prepare_koina_spectrum_for_plot(predicted_frame, **annotation_kwargs)
    return (reference_psm, predicted_reference_frame)


def _annotate_query_spectrum(streaming_result, annotation_kwargs):
    """Match the query spectrum's peaks against theoretical fragment ions; return the annotated frame.
    The query spectrum is the carrier PSM (peaks on psm.spectrum).
    """
    psm = streaming_result.query_psm
    return prepare_spectrum_for_plot(
        psm.spectrum["mz"],
        psm.spectrum["intensity"],
        psm.peptidoform.modified_sequence,
        psm.get_precursor_charge() or 0,
        precursor_mz=psm.precursor_mz or 0.0,
        **annotation_kwargs,
    )


def _cache_reference_spectra(streaming_result, reference_store, annotation_kwargs):
    """Annotate this peptidoform's reference spectra and cache them in reference_store."""
    peptidoform = streaming_result.peptidoform
    # Guard: skip if this peptidoform is already cached or has no references to annotate
    if peptidoform in reference_store or not streaming_result.reference_psms:
        return
    # Annotate each experimental reference against its own theoretical fragments;
    # each reference is itself the carrier PSM.
    annotated_references = []
    for reference_psm in streaming_result.reference_psms:
        reference_frame = prepare_spectrum_for_plot(
            reference_psm.spectrum["mz"], reference_psm.spectrum["intensity"],
            reference_psm.peptidoform.modified_sequence, reference_psm.get_precursor_charge() or 0,
            precursor_mz=reference_psm.precursor_mz or 0.0,
            **annotation_kwargs,
        )
        annotated_references.append((reference_psm, reference_frame))
    reference_store[peptidoform] = annotated_references


def _try_koina_fallback(streaming_result, reference_store, annotation_kwargs,
                        koina_model, koina_collision_energy, koina_instrument_type):
    """When no experimental reference exists for this peptidoform, request a Koina prediction; return the resolved list.
    """
    peptidoform = streaming_result.peptidoform
    # Prefer experimental references: only fall back to Koina when the cache is empty
    references = reference_store.get(peptidoform, [])
    if references:
        return references
    # Request an AI-predicted reference for this peptidoform + charge
    koina_ref = get_koina_prediction(
        streaming_result.query_psm.peptidoform.modified_sequence,
        streaming_result.query_psm.get_precursor_charge() or 0,
        annotation_kwargs,
        koina_model=koina_model,
        koina_collision_energy=koina_collision_energy,
        koina_instrument_type=koina_instrument_type,
    )
    # On success, cache the synthetic reference; on failure, return the empty list
    if koina_ref is not None:
        reference_store[peptidoform] = [koina_ref]
        return reference_store[peptidoform]
    return references


# ============================================================================
# Label/lookup helpers (also imported by streamlit_app.views.explore.mirror_plots)
# ============================================================================

def format_spectrum_labels(query_psm, reference_psm=None) -> tuple[str, str]:
    """Build the human-readable plot labels for a (query, reference) pair.
    Returns (query_label, reference_label), each the source file name with its RT
    appended, e.g. "Query: run01.mzML (RT: 42.3)". The reference label is ""
    when reference_psm is None.
    """
    query_label = f"Query: {source_label(query_psm)} (RT: {query_psm.retention_time:.1f})"
    # Guard clause: an unpaired query has no reference label to build
    if reference_psm is None:
        return query_label, ""
    reference_label = f"Ref: {source_label(reference_psm)} (RT: {reference_psm.retention_time:.1f})"
    return query_label, reference_label


def build_psm_to_peptidoform_mapping(
    peptidoform_psms: dict, annotated_queries: list[tuple],
) -> dict[str, str]:
    """Map each PSM id to its peptidoform key; returns {psm_id: peptidoform}.
    annotated_queries is the canonical [(psm_id, query_psm, query_frame), ...] shape produced by
    run_streaming_annotation; callers holding another representation reshape to it first. A PSM
    whose spectrum key has no annotated query is simply omitted from the result.
    """
    # Invert to a per-PSM-key lookup so each peptidoform PSM can find its PSM id.
    # The key is the full psm_lookup_key ("spectrum_id|peptidoform"), not a bare spectrum id.
    lookup_key_to_psm_id = {
        psm_lookup_key(query_psm): psm_id
        for psm_id, query_psm, _query_frame in annotated_queries
    }

    # Walk peptidoforms: rebuild each PSM's lookup key and resolve it to a PSM id
    psm_id_to_peptidoform: dict[str, str] = {}
    for peptidoform, psm_list in peptidoform_psms.items():
        for psm in psm_list:
            psm_id = lookup_key_to_psm_id.get(psm_lookup_key(psm))
            if psm_id is not None:
                psm_id_to_peptidoform[psm_id] = peptidoform
    return psm_id_to_peptidoform


# ============================================================================
# Per-query orchestration
# ============================================================================

def process_query_result(streaming_result, psm_lookup, annotation_kwargs,
                         reference_store, koina_enabled=False,
                         koina_model: str = DEFAULTS.koina_model,
                         koina_collision_energy: int = DEFAULTS.koina_collision_energy,
                         koina_instrument_type: str = DEFAULTS.koina_instrument_type):
    """Annotate one query spectrum, resolve its references, and score every pair.
    References come from the experimental reference_store cache, with an optional Koina AI fallback
    when no experimental reference exists. Returns (psm_id, query_psm, query_frame, comparison_rows) —
    one metrics dict per (query, reference) pair — or None when the query has no matching PSM.
    """
    # PSM lookup: skip spectra that have no matching PSM in the lookup results
    psm_id = psm_lookup.get(psm_lookup_key(streaming_result.query_psm))
    if psm_id is None:
        return None

    # Annotate query + warm reference cache for this peptidoform (once)
    query_frame = _annotate_query_spectrum(streaming_result, annotation_kwargs)
    _cache_reference_spectra(streaming_result, reference_store, annotation_kwargs)

    # Resolve references: experimental cache first, optional Koina AI fallback
    if koina_enabled:
        references = _try_koina_fallback(
            streaming_result, reference_store, annotation_kwargs,
            koina_model, koina_collision_energy, koina_instrument_type,
        )
    else:
        references = reference_store.get(streaming_result.peptidoform, [])

    # Compute metrics: one comparison row per (query, reference) pair. The query
    # spectrum and each reference are carrier PSMs read directly.
    query_psm = streaming_result.query_psm
    base_sequence = query_psm.peptidoform.sequence.upper()
    # Modification column is blank for the unmodified peptidoform
    modification = streaming_result.peptidoform if query_psm.peptidoform.is_modified else ""
    comparison_rows = []
    for reference_psm, reference_frame in references:
        metrics = compute_all_metrics(query_frame, reference_frame)
        # Identity columns: per-PSM metadata + query/reference labels and RTs
        row: dict = {
            "psm_id": psm_id,
            "Peptide": base_sequence,
            "Modification": modification,
            "Charge": query_psm.get_precursor_charge() or 0,
            "Query": source_label(query_psm),
            "Reference": source_label(reference_psm),
            "query_spectrum_id": query_psm.spectrum_id,
            "ref_spectrum_id": reference_psm.spectrum_id,
            "query_rt": query_psm.retention_time,
            "ref_rt": reference_psm.retention_time,
            "rt_delta": rt_delta(query_psm.retention_time, reference_psm.retention_time),
        }
        # Metric columns: metric key -> spreadsheet column name, NaN for missing values
        for metric_key, column_name in METRIC_KEY_TO_COLUMN_NAME.items():
            row[column_name] = metrics.get(metric_key, float("nan"))
        # Count columns: peak/annotation tallies on both sides plus shared
        for count_key in COMPARISON_COUNT_COLUMNS:
            row[count_key] = metrics.get(count_key, 0)
        comparison_rows.append(row)
    return (psm_id, query_psm, query_frame, comparison_rows)


# ============================================================================
# Streaming-annotation orchestrator (load stage)
# ============================================================================

def make_query_ready_callback(
    *,
    psm_lookup: dict,
    annotation_kwargs: dict,
    reference_store: dict,
    koina_kwargs: dict,
    annotated_queries: list,
    comparison_rows: list,
    on_progress: Callable[[int], None],
):
    """Build the per-PSM callback fed to load_spectra_streaming as on_query_ready.
    The returned closure annotates each arriving query, appends results to the
    caller-owned annotated_queries and comparison_rows buffers, and ticks
    on_progress with the running loaded-PSM count.
    """
    # Mutable-cell closure idiom: a one-element list whose single slot the nested
    # closure can mutate (a plain int would be read-only from inside on_query_ready).
    loaded_count = [0]

    def on_query_ready(streaming_result):
        """Annotate one streamed query spectrum and record it into the shared buffers.
        Invoked once per query by load_spectra_streaming; delegates to
        process_query_result. A None result (no matching PSM) is skipped; otherwise the
        annotated query and its metric rows are buffered and the progress counter advances.
        """
        # Annotate the query, resolve its references, and score every pair.
        annotation_outcome = process_query_result(
            streaming_result, psm_lookup, annotation_kwargs, reference_store, **koina_kwargs,
        )
        # Skip spectra with no matching PSM in the lookup results.
        if annotation_outcome is None:
            return
        # Buffer the annotated query and its per-pair metric rows for the caller.
        psm_id, query_psm, query_frame, rows = annotation_outcome
        annotated_queries.append((psm_id, query_psm, query_frame))
        comparison_rows.extend(rows)
        # Advance the running count and report progress to the frontend.
        loaded_count[0] += 1
        on_progress(loaded_count[0])

    return on_query_ready


def run_streaming_annotation(
    *,
    peptidoform_psms: dict,
    mzml_paths,
    reference_substring: str,
    psm_lookup: dict,
    annotation_kwargs: dict,
    koina_kwargs: dict,
    workers: int,
    on_progress: Optional[Callable[[int], None]] = None,
    on_file_start: Optional[Callable] = None,
    on_file_complete: Optional[Callable] = None,
    on_error: Optional[Callable[[str], None]] = None,
) -> tuple[list, dict, list]:
    """Stream-load and annotate every PSM, returning the three load-stage buffers.
    Wires a per-PSM callback into load_spectra_streaming and runs the full streaming load; the
    optional on_* callbacks track progress, file lifecycle, and errors (each None falls back to a
    no-op). Returns (annotated_queries, reference_store, comparison_rows): the
    (psm_id, query_psm, query_frame) tuples, the reference cache keyed by peptidoform, and the metric dicts.
    """
    # Init buffers: shared accumulators the on_query_ready callback writes into;
    # passing them by reference lets the callback mutate them as spectra arrive.
    annotated_queries: list = []
    reference_store: dict = {}
    comparison_rows: list = []
    # No-op stand-in for any optional callback the caller did not supply.
    noop = lambda *_: None

    # Build callback: per-PSM closure that annotates and appends to buffers
    on_query_ready = make_query_ready_callback(
        psm_lookup=psm_lookup,
        annotation_kwargs=annotation_kwargs,
        reference_store=reference_store,
        koina_kwargs=koina_kwargs,
        annotated_queries=annotated_queries,
        comparison_rows=comparison_rows,
        on_progress=on_progress or noop,
    )

    # Dispatch: stream every mzML, fanning lifecycle callbacks to the frontend
    load_spectra_streaming(
        peptidoform_psms=peptidoform_psms,
        mzml_paths=[Path(path) for path in mzml_paths],
        reference_substring=reference_substring,
        on_reference_loaded=noop,
        on_query_ready=on_query_ready,
        on_file_start=on_file_start or noop,
        on_file_complete=on_file_complete or noop,
        error_callback=on_error or noop,
        workers=workers,
    )

    return annotated_queries, reference_store, comparison_rows


# ============================================================================
# Plot-task builder (export stage)
# ============================================================================

def peptidoform_mz_axis_maxima(annotated_queries, reference_store,
                               psm_id_to_peptidoform) -> dict[str, float]:
    """Map each peptidoform to the one x-axis upper bound all its plots share (#114).
    A group's bound covers every spectrum it can be plotted from — all of its query panels at
    any charge plus all of its references — so the axis stays put whichever pair is drawn, and
    the peptide's plots can be read side by side (in the Explore view or the export gallery)
    without silently rescaling. annotated_queries is the canonical [(psm_id, psm, frame)] list.
    """
    frames_by_peptidoform: dict[str, list] = {}
    for psm_id, _query_psm, query_frame in annotated_queries:
        peptidoform = psm_id_to_peptidoform.get(psm_id)
        if peptidoform:
            frames_by_peptidoform.setdefault(peptidoform, []).append(query_frame)
    for peptidoform, references in reference_store.items():
        frames_by_peptidoform.setdefault(peptidoform, []).extend(
            reference_frame for _reference_psm, reference_frame in references)
    return {peptidoform: group_mz_axis_max(frames)
            for peptidoform, frames in frames_by_peptidoform.items()}


def build_plot_tasks(annotated_queries, reference_store, peptidoform_psms,
                     sanitize_fn, format_dirs=None, figure_size=(12, 8),
                     export_resolution=300, html_dir="", plot_unpaired_queries=False,
                     max_peak_labels=0):
    """Build the list of ComparisonTask plot jobs for a set of annotated queries.
    Emits one "mirror" task per (query, reference) pair; a query with no reference yields a
    "single" task only when plot_unpaired_queries is set. The metrics overlay text is pre-computed
    here so the parallel plot workers never recompute it. Every task of one peptidoform also carries
    that group's shared x-axis bound. Returns [ComparisonTask, ...] ordered by psm_id.
    """
    # Setup: normalise format_dirs, derive export_formats, build PSM->peptidoform lookup
    if format_dirs is None:
        format_dirs = {}
    export_formats = [*format_dirs] + (["html"] if html_dir else [])
    psm_id_to_peptidoform = build_psm_to_peptidoform_mapping(peptidoform_psms, annotated_queries)
    png_dir = format_dirs.get("png", "")
    svg_dir = format_dirs.get("svg", "")
    pdf_dir = format_dirs.get("pdf", "")
    axis_max_by_peptidoform = peptidoform_mz_axis_maxima(
        annotated_queries, reference_store, psm_id_to_peptidoform)

    # Build tasks: iterate PSMs in stable order, mirror when refs exist else single
    tasks: list[ComparisonTask] = []
    task_id = 0
    for psm_id, query_psm, query_frame in sorted(annotated_queries, key=lambda x: x[0]):
        # Identity fields off the carrier PSM; source_label still reads the spectrum.
        peptidoform = query_psm.peptidoform.modified_sequence
        charge = query_psm.get_precursor_charge() or 0
        peptidoform = psm_id_to_peptidoform.get(psm_id)
        references = reference_store.get(peptidoform, []) if peptidoform else []
        folder_name = sanitize_fn(peptidoform)
        plot_title = f"{peptidoform} +{charge}"
        query_src = source_label(query_psm)

        if references:
            # Mirror tasks: one per reference, with metrics + HTML output
            for reference_psm, reference_frame in references:
                task_id += 1
                metrics_text = format_metrics_text(compute_all_metrics(query_frame, reference_frame))
                query_label, reference_label = format_spectrum_labels(query_psm, reference_psm)
                # Include both spectrum ids so every (query PSM, reference PSM) pair gets its
                # own file. source_label is only the run name, so multiple PSMs of one peptide
                # in the same file pair would otherwise collide onto one PNG and overwrite each
                # other — leaving fewer PNGs than comparison rows / boxplot points.
                safe_filename = sanitize_fn(
                    f"{query_src}_{query_psm.spectrum_id}"
                    f"_vs_{source_label(reference_psm)}_{reference_psm.spectrum_id}"
                )
                html_path = (
                    str(Path(html_dir) / folder_name / f"{safe_filename}.html")
                    if html_dir else ""
                )
                tasks.append(ComparisonTask(
                    task_id=task_id, task_type="mirror",
                    query_frame=query_frame, reference_frame=reference_frame,
                    query_label=query_label, reference_label=reference_label,
                    plot_title=plot_title, folder_name=folder_name, safe_filename=safe_filename,
                    png_dir=png_dir, svg_dir=svg_dir, pdf_dir=pdf_dir,
                    export_formats=export_formats,
                    figure_size=figure_size, export_resolution=export_resolution,
                    max_peak_labels=max_peak_labels,
                    mz_axis_max=axis_max_by_peptidoform.get(peptidoform, 0.0),
                    html_plot_path=html_path, metrics_text=metrics_text,
                    peptide=peptidoform, charge=charge,
                    query_rt=query_psm.retention_time, reference_rt=reference_psm.retention_time,
                ))
        elif plot_unpaired_queries:
            # Single task: unpaired query, no metrics overlay, no HTML
            task_id += 1
            tasks.append(ComparisonTask(
                task_id=task_id, task_type="single",
                query_frame=query_frame, reference_frame=None,
                query_label=f"{query_src} (RT: {query_psm.retention_time:.1f})", reference_label="",
                plot_title=plot_title, folder_name=folder_name,
                safe_filename=sanitize_fn(f"{query_src}_{query_psm.spectrum_id}"),
                png_dir=png_dir, svg_dir=svg_dir, pdf_dir=pdf_dir,
                export_formats=[f for f in export_formats if f != "html"],
                figure_size=figure_size, export_resolution=export_resolution,
                max_peak_labels=max_peak_labels,
                mz_axis_max=axis_max_by_peptidoform.get(peptidoform, 0.0),
                html_plot_path="", metrics_text="",
                peptide=peptidoform, charge=charge,
                query_rt=query_psm.retention_time,
            ))
    return tasks
