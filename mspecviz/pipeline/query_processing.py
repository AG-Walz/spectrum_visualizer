"""Pipeline orchestration: query processing + plot-task building.

Used by both ``cli/main.py`` (terminal mode) and the Streamlit app in
``streamlit_app/views/`` (interactive mode). Keeps shared helpers in a single
authoritative location so both frontends stay in sync.
"""
from pathlib import Path

from mspecviz.annotation.comparison_metrics import (
    compute_all_metrics,
    format_metrics_text,
)
from mspecviz.common.defaults import DEFAULTS
from mspecviz.pipeline.annotation_steps import (
    _annotate_query_spectrum,
    _build_metric_rows,
    _cache_reference_spectra,
    _try_koina_fallback,
)
from mspecviz.plotting.plot_task_dataclasses import ComparisonTask
from mspecviz.spectra.loaded_spectrum import source_label, spectrum_lookup_key


def format_spectrum_labels(query_spectrum, reference_spectrum=None) -> tuple[str, str]:
    """Format query/reference spectrum labels with RT.

    Returns ``(query_label, reference_label)`` tuple.  *reference_label* is ``""``
    when *reference_spectrum* is ``None``.
    """
    query_label = f"Query: {source_label(query_spectrum)} (RT: {query_spectrum.rt:.1f})"
    if reference_spectrum is None:
        return query_label, ""
    reference_label = f"Ref: {source_label(reference_spectrum)} (RT: {reference_spectrum.rt:.1f})"
    return query_label, reference_label


def build_psm_to_variant_mapping(variant_groups, annotated_queries) -> dict[str, str]:
    """Build mapping from PSM ID to variant key.

    *annotated_queries* can be:

    * ``dict``  -- ``{psm_id: (spec, df)}``   (Streamlit format)
    * ``list``  -- ``[(psm_id, spec, df)]``    (CLI format)
    """
    # Normalise to iterable of (psm_id, spec); pre-build O(1) lookup by spectrum key.
    if isinstance(annotated_queries, dict):
        items = [(psm_id, spec) for psm_id, (spec, _) in annotated_queries.items()]
    else:
        items = [(psm_id, spec) for psm_id, spec, _ in annotated_queries]
    spec_to_psm_id = {spectrum_lookup_key(spec): psm_id for psm_id, spec in items}

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
    """Annotate query, cache references, compute metrics. Returns (psm_id, spec, df, rows) or None."""
    psm_id = psm_lookup.get(spectrum_lookup_key(result.query_spectrum))
    if psm_id is None:
        return None

    query_frame = _annotate_query_spectrum(result, annotation_kwargs)
    _cache_reference_spectra(result, reference_store, annotation_kwargs)

    if koina_enabled:
        references = _try_koina_fallback(
            result, reference_store, annotation_kwargs,
            koina_model, koina_collision_energy, koina_instrument_type,
        )
    else:
        references = reference_store.get(result.variant_key, [])

    comparison_rows = _build_metric_rows(psm_id, result, query_frame, references)
    return (psm_id, result.query_spectrum, query_frame, comparison_rows)


def _build_task(
    *, task_id: int, query_spectrum, query_frame,
    reference_spectrum=None, reference_frame=None,
    folder_name: str, plot_title: str, sanitize_fn,
    format_dirs: dict, export_formats: list[str], html_dir: str,
    figure_size, export_resolution: int, normalize: bool,
) -> ComparisonTask:
    """Construct a ComparisonTask; ``reference_spectrum=None`` yields a single-plot task."""
    is_mirror = reference_spectrum is not None
    if is_mirror:
        metrics_text = format_metrics_text(compute_all_metrics(query_frame, reference_frame))
        query_label, reference_label = format_spectrum_labels(query_spectrum, reference_spectrum)
        safe_filename = sanitize_fn(
            f"{source_label(query_spectrum)}_vs_{source_label(reference_spectrum)}"
        )
        html_path = (
            str(Path(html_dir) / folder_name / f"{safe_filename}.html")
            if html_dir else ""
        )
        formats = export_formats
    else:
        # Unpaired/single plots never emit HTML and skip metrics overlay.
        metrics_text = ""
        query_label = f"{source_label(query_spectrum)} (RT: {query_spectrum.rt:.1f})"
        reference_label = ""
        safe_filename = sanitize_fn(source_label(query_spectrum))
        html_path = ""
        formats = [f for f in export_formats if f != "html"]

    return ComparisonTask(
        task_id=task_id, task_type="mirror" if is_mirror else "single",
        query_frame=query_frame, reference_frame=reference_frame,
        query_label=query_label, reference_label=reference_label,
        plot_title=plot_title, folder_name=folder_name, safe_filename=safe_filename,
        png_dir=format_dirs.get("png", ""), svg_dir=format_dirs.get("svg", ""),
        pdf_dir=format_dirs.get("pdf", ""),
        export_formats=formats,
        figure_size=figure_size, export_resolution=export_resolution,
        normalize=normalize, html_plot_path=html_path, metrics_text=metrics_text,
        peptide=query_spectrum.peptidoform, charge=query_spectrum.precursor_charge,
    )


def build_plot_tasks(annotated_queries, reference_store, variant_groups,
                     normalize, sanitize_fn, format_dirs=None, figure_size=(12, 8),
                     export_resolution=300, html_dir="", plot_unpaired_queries=False):
    """Build list of ``ComparisonTask`` instances for parallel/sequential generation.

    Metrics + ``metrics_text`` are pre-computed here (in the main process) so
    workers can re-use them via :pyattr:`ComparisonTask.metrics_text` instead
    of duplicating the work per task.
    """
    if format_dirs is None:
        format_dirs = {}
    export_formats = [*format_dirs] + (["html"] if html_dir else [])
    psm_id_to_variant_key = build_psm_to_variant_mapping(variant_groups, annotated_queries)

    tasks: list[ComparisonTask] = []
    task_id = 0
    for psm_id, query_spectrum, query_frame in sorted(annotated_queries, key=lambda x: x[0]):
        variant_key = psm_id_to_variant_key.get(psm_id)
        references = reference_store.get(variant_key, []) if variant_key else []
        folder_name = sanitize_fn(query_spectrum.peptidoform)
        plot_title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"

        if references:
            for reference_spectrum, reference_frame in references:
                task_id += 1
                tasks.append(_build_task(
                    task_id=task_id, query_spectrum=query_spectrum, query_frame=query_frame,
                    reference_spectrum=reference_spectrum, reference_frame=reference_frame,
                    folder_name=folder_name, plot_title=plot_title, sanitize_fn=sanitize_fn,
                    format_dirs=format_dirs, export_formats=export_formats, html_dir=html_dir,
                    figure_size=figure_size, export_resolution=export_resolution,
                    normalize=normalize,
                ))
        elif plot_unpaired_queries:
            task_id += 1
            tasks.append(_build_task(
                task_id=task_id, query_spectrum=query_spectrum, query_frame=query_frame,
                folder_name=folder_name, plot_title=plot_title, sanitize_fn=sanitize_fn,
                format_dirs=format_dirs, export_formats=export_formats, html_dir=html_dir,
                figure_size=figure_size, export_resolution=export_resolution,
                normalize=normalize,
            ))
    return tasks
