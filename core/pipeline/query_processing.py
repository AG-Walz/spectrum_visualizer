"""Per-query workflow (annotate, match, score) plus plot-task building (called by ``pipeline.streaming_orchestration``, ``cli.main``, and ``streamlit_app.views``)."""
from pathlib import Path

from core.annotation.comparison_metrics import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    format_metrics_text,
)
from core.defaults import DEFAULTS
from core.pipeline.annotation_steps import (
    _annotate_query_spectrum,
    _cache_reference_spectra,
    _try_koina_fallback,
)
from core.plotting.plot_task_dataclasses import ComparisonTask
from core.spectra.loaded_spectrum import source_label, spectrum_lookup_key


def format_spectrum_labels(query_spectrum, reference_spectrum=None) -> tuple[str, str]:
    """Build ``"Query: <file> (RT: ...)"`` and matching reference label; reference label is ``""`` when no reference is given."""
    query_label = f"Query: {source_label(query_spectrum)} (RT: {query_spectrum.rt:.1f})"
    if reference_spectrum is None:
        return query_label, ""
    reference_label = f"Ref: {source_label(reference_spectrum)} (RT: {reference_spectrum.rt:.1f})"
    return query_label, reference_label


def build_psm_to_variant_mapping(variant_groups, annotated_queries) -> dict[str, str]:
    """Map each PSM id to its peptide-variant key; accepts the Streamlit dict shape or the CLI list-of-tuples shape."""
    # Normalise input: handle Streamlit's dict shape and the CLI list-of-tuples shape uniformly
    if isinstance(annotated_queries, dict):
        items = [(psm_id, spec) for psm_id, (spec, _) in annotated_queries.items()]
    else:
        items = [(psm_id, spec) for psm_id, spec, _ in annotated_queries]
    spec_to_psm_id = {spectrum_lookup_key(spec): psm_id for psm_id, spec in items}

    # Walk variants: match each PSM's spectrum key back to the annotated query's PSM id
    psm_id_to_variant_key: dict[str, str] = {}
    for variant_key, psm_list in variant_groups.items():
        for psm in psm_list:
            key = f"{psm.spectrum_id}|{psm.peptidoform.modified_sequence}"
            psm_id = spec_to_psm_id.get(key)
            if psm_id is not None:
                psm_id_to_variant_key[psm_id] = variant_key
    return psm_id_to_variant_key


def process_query_result(result, psm_lookup, annotation_kwargs,
                         reference_store, koina_enabled=False,
                         koina_model: str = DEFAULTS.koina_model,
                         koina_collision_energy: int = DEFAULTS.koina_collision_energy,
                         koina_instrument_type: str = DEFAULTS.koina_instrument_type):
    """Annotate one query, resolve its references (experimental or Koina), and score every pair; return ``(psm_id, spectrum, frame, rows)`` or ``None``."""
    # PSM lookup: skip spectra that have no matching PSM in the search results
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
    """Build one ``ComparisonTask`` per (query, reference) pair; metrics text is pre-computed here so plot workers don't recompute it."""
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
