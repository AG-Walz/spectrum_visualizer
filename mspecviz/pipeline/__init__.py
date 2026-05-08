"""Shared pipeline logic for CLI and Streamlit frontends.

Re-exports the public API of the submodules so callers can use
``from mspecviz.pipeline import ...`` without caring which file a name lives in.
"""
from mspecviz.pipeline.query_processing import (
    build_plot_tasks,
    build_psm_to_variant_mapping,
    format_spectrum_labels,
    process_query_result,
)

__all__ = [
    "build_plot_tasks",
    "build_psm_to_variant_mapping",
    "format_spectrum_labels",
    "process_query_result",
]
