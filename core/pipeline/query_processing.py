"""Per-query workflow layer: annotate one query spectrum, resolve its references,
score every (query, reference) pair, and package the results either as comparison
rows or as ComparisonTask plot jobs.

`process_query_result` is the per-query step of the streaming pipeline; it runs
through the annotation phases (annotate query, cache experimental references,
optional Koina AI fallback) and the scoring via core.annotation.comparison_metrics.
`build_plot_tasks` converts already-annotated queries into one plot task per pair,
pre-computing the metrics overlay so the parallel plot workers never recompute it.

process_query_result called by core.pipeline.streaming_orchestration;
build_plot_tasks called by core.export.run_export;
build_psm_to_variant_mapping and format_spectrum_labels also used by
streamlit_app.views.explore.mirror_plots.
"""
from pathlib import Path

import numpy as np
from psm_utils import Peptidoform

from core.annotation.comparison_metrics import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    format_metrics_text,
)
from core.annotation.peak_matching import (
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.defaults import DEFAULTS
from core.plotting.plot_worker import ComparisonTask
from core.spectra.koina_prediction import predict_koina_spectrum
from core.spectra.loaded_spectrum import LoadedSpectrum, source_label, spectrum_lookup_key


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
    # Predict: call Koina with model + CE + instrument, bail on empty frame
    predicted_frame = predict_koina_spectrum(
        peptidoform_str=peptidoform,
        charge=charge,
        model_id=koina_model,
        collision_energy=koina_collision_energy,
        instrument_type=koina_instrument_type,
    )
    if predicted_frame is None or predicted_frame.empty:
        return None

    # Wrap prediction: build a synthetic LoadedSpectrum tagged as reference
    reference_spectrum = LoadedSpectrum(
        spectrum_id="koina_prediction",
        mz=predicted_frame["mz"].to_numpy(dtype=np.float64),
        intensity=predicted_frame["intensity"].to_numpy(dtype=np.float64),
        precursor_mz=0.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence=Peptidoform(peptidoform).sequence,
        rt=0.0,
        source_file=f"Koina:{koina_model}",
        is_reference=True,
    )

    # Annotate: reuse Koina's native AI-predicted ion labels directly
    reference_frame = prepare_koina_spectrum_for_plot(predicted_frame, **annotation_kwargs)
    return (reference_spectrum, reference_frame)


def _annotate_query_spectrum(result, annotation_kwargs):
    """Match the query spectrum's peaks against theoretical fragment ions; return the annotated frame.
    """
    return prepare_spectrum_for_plot(
        result.query_spectrum.mz,
        result.query_spectrum.intensity,
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        precursor_mz=result.query_spectrum.precursor_mz,
        **annotation_kwargs,
    )


def _cache_reference_spectra(result, reference_store, annotation_kwargs):
    """Annotate every reference spectrum for this peptide variant and stash them in
    reference_store (once per variant).
    """
    variant_key = result.variant_key
    # Guard: skip if this variant is already cached or has no references to annotate
    if variant_key in reference_store or not result.reference_spectra:
        return
    # Annotate each experimental reference against its own theoretical fragments
    ref_list = []
    for reference_spectrum in result.reference_spectra:
        reference_frame = prepare_spectrum_for_plot(
            reference_spectrum.mz, reference_spectrum.intensity,
            reference_spectrum.peptidoform, reference_spectrum.precursor_charge,
            precursor_mz=reference_spectrum.precursor_mz,
            **annotation_kwargs,
        )
        ref_list.append((reference_spectrum, reference_frame))
    reference_store[variant_key] = ref_list


def _try_koina_fallback(result, reference_store, annotation_kwargs,
                        koina_model, koina_collision_energy, koina_instrument_type):
    """When no experimental reference exists for this variant, request a Koina prediction; return the resolved list.
    """
    variant_key = result.variant_key
    # Prefer experimental references: only fall back to Koina when the cache is empty
    references = reference_store.get(variant_key, [])
    if references:
        return references
    # Request an AI-predicted reference for this peptidoform + charge
    koina_ref = get_koina_prediction(
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        annotation_kwargs,
        koina_model=koina_model,
        koina_collision_energy=koina_collision_energy,
        koina_instrument_type=koina_instrument_type,
    )
    # On success, cache the synthetic reference; on failure, return the empty list
    if koina_ref is not None:
        reference_store[variant_key] = [koina_ref]
        return reference_store[variant_key]
    return references


# ============================================================================
# Label/lookup helpers (also imported by streamlit_app.views.explore.mirror_plots)
# ============================================================================

def format_spectrum_labels(query_spectrum, reference_spectrum=None) -> tuple[str, str]:
    """Build the human-readable plot labels for a (query, reference) pair.
    Returns (query_label, reference_label); each label is the source file
    name with its retention time appended, e.g. "Query: run01.mzML (RT: 42.3)".
    The reference label is "" when reference_spectrum is None (a
    single, unpaired query).
    """
    query_label = f"Query: {source_label(query_spectrum)} (RT: {query_spectrum.rt:.1f})"
    # Guard clause: an unpaired query has no reference label to build
    if reference_spectrum is None:
        return query_label, ""
    reference_label = f"Ref: {source_label(reference_spectrum)} (RT: {reference_spectrum.rt:.1f})"
    return query_label, reference_label


def build_psm_to_variant_mapping(variant_groups, annotated_queries) -> dict[str, str]:
    """Map each PSM id to its peptide-variant key.
    Returns {psm_id: variant_key}. annotated_queries may arrive in two
    shapes: the Streamlit dict {psm_id: (spectrum, frame)} or the CLI
    list-of-tuples [(psm_id, spectrum, frame), ...]; both are normalised
    before matching. A PSM whose spectrum key has no annotated query is simply
    omitted from the result.
    """
    # Normalise input: handle Streamlit's dict shape and the CLI list-of-tuples shape uniformly
    if isinstance(annotated_queries, dict):
        items = [(psm_id, spec) for psm_id, (spec, _) in annotated_queries.items()]
    else:
        items = [(psm_id, spec) for psm_id, spec, _ in annotated_queries]
    # Invert to a spectrum-key lookup so each variant PSM can find its PSM id
    spec_to_psm_id = {spectrum_lookup_key(spec): psm_id for psm_id, spec in items}

    # Walk variants: rebuild each PSM's spectrum key and resolve it to a PSM id
    psm_id_to_variant_key: dict[str, str] = {}
    for variant_key, psm_list in variant_groups.items():
        for psm in psm_list:
            key = f"{psm.spectrum_id}|{psm.peptidoform.modified_sequence}"
            psm_id = spec_to_psm_id.get(key)
            if psm_id is not None:
                psm_id_to_variant_key[psm_id] = variant_key
    return psm_id_to_variant_key


# ============================================================================
# Per-query orchestration + plot-task builder
# ============================================================================

def process_query_result(result, psm_lookup, annotation_kwargs,
                         reference_store, koina_enabled=False,
                         koina_model: str = DEFAULTS.koina_model,
                         koina_collision_energy: int = DEFAULTS.koina_collision_energy,
                         koina_instrument_type: str = DEFAULTS.koina_instrument_type):
    """Annotate one query spectrum, resolve its references, and score every pair.
    Resolves references from the experimental reference_store cache, with an
    optional Koina AI fallback when no experimental reference exists. Returns
    (psm_id, query_spectrum, query_frame, comparison_rows), where
    comparison_rows holds one metrics dict per (query, reference) pair, or
    None when the query has no matching PSM in the lookup results.
    """
    # PSM lookup: skip spectra that have no matching PSM in the lookup results
    psm_id = psm_lookup.get(spectrum_lookup_key(result.query_spectrum))
    if psm_id is None:
        return None

    # Annotate query + warm reference cache for this variant (once)
    query_frame = _annotate_query_spectrum(result, annotation_kwargs)
    _cache_reference_spectra(result, reference_store, annotation_kwargs)

    # Resolve references: experimental cache first, optional Koina AI fallback
    if koina_enabled:
        references = _try_koina_fallback(
            result, reference_store, annotation_kwargs,
            koina_model, koina_collision_energy, koina_instrument_type,
        )
    else:
        references = reference_store.get(result.variant_key, [])

    # Compute metrics: one comparison row per (query, reference) pair
    query_spectrum = result.query_spectrum
    base_sequence = query_spectrum.base_sequence
    # Modification column is blank for the unmodified variant (variant key == bare sequence)
    modification = result.variant_key if result.variant_key != base_sequence else ""
    comparison_rows = []
    for reference_spectrum, reference_frame in references:
        metrics = compute_all_metrics(query_frame, reference_frame)
        # Identity columns: per-PSM metadata + query/reference labels and RTs
        row: dict = {
            "psm_id": psm_id,
            "Peptide": base_sequence,
            "Modification": modification,
            "Charge": query_spectrum.precursor_charge,
            "Query": source_label(query_spectrum),
            "Reference": source_label(reference_spectrum),
            "query_spectrum_id": query_spectrum.spectrum_id,
            "ref_spectrum_id": reference_spectrum.spectrum_id,
            "query_rt": query_spectrum.rt,
            "ref_rt": reference_spectrum.rt,
        }
        # Metric columns: short keys -> display names, NaN for missing values
        for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
            row[display_name] = metrics.get(short_key, float("nan"))
        # Count columns: peak/annotation tallies on both sides plus shared
        for count_key in (
            "query_peak_count", "reference_peak_count",
            "query_annotation_count", "reference_annotation_count", "shared",
        ):
            row[count_key] = metrics.get(count_key, 0)
        comparison_rows.append(row)
    return (psm_id, query_spectrum, query_frame, comparison_rows)


def build_plot_tasks(annotated_queries, reference_store, variant_groups,
                     sanitize_fn, format_dirs=None, figure_size=(12, 8),
                     export_resolution=300, html_dir="", plot_unpaired_queries=False):
    """Build the list of ComparisonTask plot jobs for a set of annotated queries.
    Emits one "mirror" task per (query, reference) pair; a query with no
    reference yields one "single" task only when plot_unpaired_queries is
    set. The metrics overlay text is pre-computed here so the parallel plot
    workers never recompute it. Returns [ComparisonTask, ...] ordered by psm_id.
    """
    # Setup: normalise format_dirs, derive export_formats, build PSM->variant lookup
    if format_dirs is None:
        format_dirs = {}
    export_formats = [*format_dirs] + (["html"] if html_dir else [])
    psm_id_to_variant_key = build_psm_to_variant_mapping(variant_groups, annotated_queries)
    png_dir = format_dirs.get("png", "")
    svg_dir = format_dirs.get("svg", "")
    pdf_dir = format_dirs.get("pdf", "")

    # Build tasks: iterate PSMs in stable order, mirror when refs exist else single
    tasks: list[ComparisonTask] = []
    task_id = 0
    for psm_id, query_spectrum, query_frame in sorted(annotated_queries, key=lambda x: x[0]):
        variant_key = psm_id_to_variant_key.get(psm_id)
        references = reference_store.get(variant_key, []) if variant_key else []
        folder_name = sanitize_fn(query_spectrum.peptidoform)
        plot_title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"
        query_src = source_label(query_spectrum)

        if references:
            # Mirror tasks: one per reference, with metrics + HTML output
            for reference_spectrum, reference_frame in references:
                task_id += 1
                metrics_text = format_metrics_text(compute_all_metrics(query_frame, reference_frame))
                query_label, reference_label = format_spectrum_labels(query_spectrum, reference_spectrum)
                safe_filename = sanitize_fn(f"{query_src}_vs_{source_label(reference_spectrum)}")
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
                    html_plot_path=html_path, metrics_text=metrics_text,
                    peptide=query_spectrum.peptidoform, charge=query_spectrum.precursor_charge,
                ))
        elif plot_unpaired_queries:
            # Single task: unpaired query, no metrics overlay, no HTML
            task_id += 1
            tasks.append(ComparisonTask(
                task_id=task_id, task_type="single",
                query_frame=query_frame, reference_frame=None,
                query_label=f"{query_src} (RT: {query_spectrum.rt:.1f})", reference_label="",
                plot_title=plot_title, folder_name=folder_name,
                safe_filename=sanitize_fn(query_src),
                png_dir=png_dir, svg_dir=svg_dir, pdf_dir=pdf_dir,
                export_formats=[f for f in export_formats if f != "html"],
                figure_size=figure_size, export_resolution=export_resolution,
                html_plot_path="", metrics_text="",
                peptide=query_spectrum.peptidoform, charge=query_spectrum.precursor_charge,
            ))
    return tasks
