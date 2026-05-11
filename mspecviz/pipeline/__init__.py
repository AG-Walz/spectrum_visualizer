"""Shared pipeline logic for CLI and Streamlit frontends.

Re-exports the public API of the submodules so callers can use
``from mspecviz.pipeline import ...`` without caring which file a name lives in.
"""
from mspecviz.pipeline.kwargs_builders import (
    build_annotation_kwargs,
    build_koina_kwargs,
)
from mspecviz.pipeline.query_processing import (
    build_plot_tasks,
    build_psm_to_variant_mapping,
    format_spectrum_labels,
    process_query_result,
)
from mspecviz.pipeline.streaming_orchestration import (
    make_query_ready_callback,
    run_streaming_annotation,
)

__all__ = [
    "build_annotation_kwargs",
    "build_koina_kwargs",
    "build_plot_tasks",
    "build_psm_to_variant_mapping",
    "format_spectrum_labels",
    "make_query_ready_callback",
    "process_query_result",
    "run_streaming_annotation",
]
