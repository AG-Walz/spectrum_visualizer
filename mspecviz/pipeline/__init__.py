"""Shared pipeline logic for CLI and Streamlit frontends.

Re-exports the public API of the submodules so callers can keep using
``from mspecviz.pipeline import ...`` without caring which file a name
lives in.
"""
from mspecviz.pipeline.annotation_steps import (
    _annotate_query_spectrum,
    _build_metric_rows,
    _cache_reference_spectra,
    _try_koina_fallback,
    build_comparison_row,
    get_koina_prediction,
)
from mspecviz.pipeline.query_processing import (
    build_plot_tasks,
    build_psm_to_variant_mapping,
    format_spectrum_labels,
    process_query_result,
)

__all__ = [
    # annotation_steps
    "_annotate_query_spectrum",
    "_build_metric_rows",
    "_cache_reference_spectra",
    "_try_koina_fallback",
    "get_koina_prediction",
    # orchestration
    "build_comparison_row",
    "build_plot_tasks",
    "build_psm_to_variant_mapping",
    "format_spectrum_labels",
    "process_query_result",
]
