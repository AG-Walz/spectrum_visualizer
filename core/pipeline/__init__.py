"""Streaming-annotation pipeline shared by cli.main and streamlit_app.views.search_load.
"""
from core.pipeline.kwargs_builders import (
    build_annotation_kwargs,
    build_koina_kwargs,
)
from core.pipeline.query_processing import (
    build_plot_tasks,
    build_psm_to_variant_mapping,
    format_spectrum_labels,
    process_query_result,
)
from core.pipeline.streaming_orchestration import (
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
