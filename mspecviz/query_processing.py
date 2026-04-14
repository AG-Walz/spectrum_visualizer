"""
Shared processing logic for the MSpecViz pipeline.

Used by both ``cli/main.py`` (terminal mode) and the Streamlit app in
``streamlit_app/views/`` (interactive mode).  Keeps all duplicated
spectrum-processing helpers in a single authoritative location.
"""
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from psm_utils import Peptidoform

from mspecviz.koina_prediction import predict_koina_spectrum, KOINA_AVAILABLE
from mspecviz.peak_annotation import prepare_spectrum_for_plot, prepare_koina_spectrum_for_plot
from mspecviz.plot_worker import ComparisonTask
from mspecviz.spectrum_loader import spectrum_lookup_key, LoadedSpectrum
from mspecviz.spectrum_metrics import compute_all_metrics, format_metrics_text, METRIC_SHORT_TO_DISPLAY

logger = logging.getLogger("mspecviz.pipeline")


def get_source_label(spec) -> str:
    """Extract filename stem from spectrum, or 'Koina' if no source file."""
    return Path(spec.source_file).stem if spec.source_file else "Koina"


def format_spectrum_labels(query_spec, ref_spec=None):
    """Format query/reference spectrum labels with RT.

    Returns ``(query_label, ref_label)`` tuple.  *ref_label* is ``""``
    when *ref_spec* is ``None``.
    """
    query_label = f"Query: {get_source_label(query_spec)} (RT: {query_spec.rt:.1f})"
    if ref_spec is None:
        return query_label, ""
    ref_label = f"Ref: {get_source_label(ref_spec)} (RT: {ref_spec.rt:.1f})"
    return query_label, ref_label


def build_psm_to_variant_mapping(variant_groups, annotated_queries):
    """Build mapping from PSM ID to variant key.

    *annotated_queries* can be:

    * ``dict``  -- ``{psm_id: (spec, df)}``   (Streamlit format)
    * ``list``  -- ``[(psm_id, spec, df)]``    (CLI format)
    """
    psm_id_to_variant_key: Dict[str, str] = {}

    # Normalise to iterable of (psm_id, spec)
    if isinstance(annotated_queries, dict):
        items = [(psm_id, spec) for psm_id, (spec, _) in annotated_queries.items()]
    else:
        items = [(psm_id, spec) for psm_id, spec, _ in annotated_queries]

    for variant_key, variant_group in variant_groups.items():
        for match in variant_group.matches:
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for psm_id, spec in items:
                if spectrum_lookup_key(spec) == key:
                    psm_id_to_variant_key[psm_id] = variant_key
                    break
    return psm_id_to_variant_key


def get_koina_prediction(peptidoform, charge, annotation_kwargs,
                         koina_model="ms2pip_timsTOF2024",
                         koina_ce=35, koina_instrument="TIMSTOF"):
    """Get Koina AI prediction and annotate as reference spectrum.

    Returns ``(LoadedSpectrum, annotated_df)`` tuple or ``None``.
    """
    if not KOINA_AVAILABLE:
        return None

    pred_df = predict_koina_spectrum(
        peptidoform_str=peptidoform,
        charge=charge,
        model_id=koina_model,
        collision_energy=koina_ce,
        instrument_type=koina_instrument,
    )
    if pred_df is None or pred_df.empty:
        return None

    mz_arr = pred_df["mz"].values.astype(np.float64)
    int_arr = pred_df["intensity"].values.astype(np.float64)

    base_seq = Peptidoform(peptidoform.split("/")[0]).sequence

    ref_spec = LoadedSpectrum(
        spectrum_id="koina_prediction",
        mz=mz_arr,
        intensity=int_arr,
        precursor_mz=0.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence=base_seq,
        rt=0.0,
        source_file=f"Koina:{koina_model}",
        is_reference=True,
    )

    # Use Koina's native AI-predicted annotations directly
    ref_df = prepare_koina_spectrum_for_plot(pred_df, **annotation_kwargs)
    return (ref_spec, ref_df)


def build_comparison_row(psm_id, result, query_df, ref_spec, ref_df, metrics):
    """Build a single comparison row dict from query/reference pair.

    *result* is a ``QueryResult`` from ``load_spectra_streaming``.
    *metrics* is the dict returned by ``compute_all_metrics``.
    """
    row: Dict[str, Any] = {
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

    return row


def _annotate_query_spectrum(result, annotation_kwargs, psm_id):
    """Phase 1 of ``process_query_result``: annotate the query spectrum."""
    try:
        return prepare_spectrum_for_plot(
            result.query_spectrum.mz,
            result.query_spectrum.intensity,
            result.query_spectrum.peptidoform,
            result.query_spectrum.precursor_charge,
            **annotation_kwargs,
        )
    except Exception as e:
        logger.warning(f"Annotation failed psm_id={psm_id}: {e}")
        return None


def _cache_reference_spectra(result, reference_store, annotation_kwargs):
    """Phase 2: annotate and cache references for this variant (once)."""
    variant_key = result.variant_key
    if variant_key in reference_store or not result.reference_spectra:
        return
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
    reference_store[variant_key] = ref_list


def _try_koina_fallback(result, reference_store, annotation_kwargs,
                        koina_model, koina_ce, koina_instrument):
    """Phase 3: if no experimental refs, try a Koina AI prediction. Returns the ref list."""
    variant_key = result.variant_key
    refs = reference_store.get(variant_key, [])
    if refs:
        return refs
    koina_ref = get_koina_prediction(
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        annotation_kwargs,
        koina_model=koina_model,
        koina_ce=koina_ce,
        koina_instrument=koina_instrument,
    )
    if koina_ref is not None:
        reference_store[variant_key] = [koina_ref]
        return reference_store[variant_key]
    return refs


def _build_metric_rows(psm_id, result, query_df, refs):
    """Phase 4: compute metrics + build summary rows for each query/ref pair."""
    rows: List[Dict[str, Any]] = []
    for ref_spec, ref_df in refs:
        try:
            metrics = compute_all_metrics(query_df, ref_df)
            rows.append(build_comparison_row(
                psm_id, result, query_df, ref_spec, ref_df, metrics,
            ))
        except Exception as e:
            logger.warning(f"Metric computation failed: {e}")
    return rows


def process_query_result(result, psm_lookup, annotation_kwargs,
                         reference_store, koina_enabled=False,
                         koina_model="ms2pip_timsTOF2024",
                         koina_ce=35, koina_instrument="TIMSTOF"):
    """Annotate query, cache refs, compute metrics. Returns (psm_id, spec, df, rows) or None.

    Orchestrates four phase helpers: ``_annotate_query_spectrum``,
    ``_cache_reference_spectra``, ``_try_koina_fallback``, ``_build_metric_rows``.
    """
    psm_id = psm_lookup.get(spectrum_lookup_key(result.query_spectrum))
    if psm_id is None:
        return None

    query_df = _annotate_query_spectrum(result, annotation_kwargs, psm_id)
    if query_df is None:
        return None

    _cache_reference_spectra(result, reference_store, annotation_kwargs)

    if koina_enabled:
        refs = _try_koina_fallback(
            result, reference_store, annotation_kwargs,
            koina_model, koina_ce, koina_instrument,
        )
    else:
        refs = reference_store.get(result.variant_key, [])

    comp_rows = _build_metric_rows(psm_id, result, query_df, refs)
    return (psm_id, result.query_spectrum, query_df, comp_rows)


def build_plot_tasks(annotated_queries, reference_store, variant_groups,
                     normalize, sanitize_fn, fmt_dirs=None, figsize=(12, 8),
                     dpi=300, html_dir="", plot_unpaired=False):
    """Build list of ``ComparisonTask`` instances for parallel/sequential generation.

    Metrics + ``metrics_text`` are pre-computed here (in the main process) so
    workers can re-use them via :pyattr:`ComparisonTask.metrics_text` instead
    of duplicating the work per task.
    """
    if fmt_dirs is None:
        fmt_dirs = {}

    export_formats = list(fmt_dirs.keys())
    if html_dir:
        export_formats = export_formats + ["html"]

    psm_id_to_variant_key = build_psm_to_variant_mapping(
        variant_groups, annotated_queries,
    )

    tasks: List[ComparisonTask] = []
    task_id = 0
    for psm_id, query_spec, query_df in sorted(annotated_queries, key=lambda x: x[0]):
        variant_key = psm_id_to_variant_key.get(psm_id)
        refs = reference_store.get(variant_key, []) if variant_key else []
        folder_name = sanitize_fn(query_spec.peptidoform)
        plot_title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"

        if refs:
            for ref_spec, ref_df in refs:
                metrics_text = format_metrics_text(compute_all_metrics(query_df, ref_df))
                query_label, ref_label = format_spectrum_labels(query_spec, ref_spec)
                safe_filename = sanitize_fn(
                    f"{get_source_label(query_spec)}_vs_{get_source_label(ref_spec)}"
                )
                html_path = (
                    str(Path(html_dir) / folder_name / f"{safe_filename}.html")
                    if html_dir else ""
                )
                task_id += 1
                tasks.append(ComparisonTask(
                    task_id=task_id,
                    task_type="mirror",
                    query_df=query_df,
                    ref_df=ref_df,
                    query_label=query_label,
                    ref_label=ref_label,
                    plot_title=plot_title,
                    folder_name=folder_name,
                    safe_filename=safe_filename,
                    png_dir=fmt_dirs.get("png", ""),
                    svg_dir=fmt_dirs.get("svg", ""),
                    pdf_dir=fmt_dirs.get("pdf", ""),
                    export_formats=export_formats,
                    figsize=figsize,
                    dpi=dpi,
                    normalize=normalize,
                    html_plot_path=html_path,
                    metrics_text=metrics_text,
                    peptide=query_spec.peptidoform,
                    charge=query_spec.precursor_charge,
                ))
        elif plot_unpaired:
            safe_filename = sanitize_fn(get_source_label(query_spec))
            query_label = f"{get_source_label(query_spec)} (RT: {query_spec.rt:.1f})"
            # Unpaired/single plots never emit HTML (cli used to enforce this).
            single_formats = [f for f in export_formats if f != "html"]
            task_id += 1
            tasks.append(ComparisonTask(
                task_id=task_id,
                task_type="single",
                query_df=query_df,
                ref_df=None,
                query_label=query_label,
                ref_label="",
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=fmt_dirs.get("png", ""),
                svg_dir=fmt_dirs.get("svg", ""),
                pdf_dir=fmt_dirs.get("pdf", ""),
                export_formats=single_formats,
                figsize=figsize,
                dpi=dpi,
                normalize=normalize,
                peptide=query_spec.peptidoform,
                charge=query_spec.precursor_charge,
            ))

    return tasks
