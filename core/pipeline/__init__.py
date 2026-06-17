"""Streaming-annotation pipeline shared by cli.main and streamlit_app.views.load_annotate.
"""
from core.pipeline.kwargs_builders import (
    build_annotation_kwargs,
    build_koina_kwargs,
)
from core.pipeline.load_annotate import (
    build_comparison_frame,
    build_plot_tasks,
    build_psm_to_peptidoform_mapping,
    format_spectrum_labels,
    make_query_ready_callback,
    process_query_result,
    run_streaming_annotation,
)

__all__ = [
    "build_annotation_kwargs",
    "build_koina_kwargs",
    "build_comparison_frame",
    "build_plot_tasks",
    "build_psm_to_peptidoform_mapping",
    "format_spectrum_labels",
    "make_query_ready_callback",
    "process_query_result",
    "run_streaming_annotation",
]
