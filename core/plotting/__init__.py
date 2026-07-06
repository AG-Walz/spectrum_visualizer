"""Plotting public API: mirror plots, summary plots, parallel render worker
(imported by core.export, cli.main, streamlit_app.views.explore).
"""
from core.plotting.mirror import create_spectrum_plot
from core.plotting.plot_worker import ComparisonTask, process_single_task
from core.plotting.summary import (
    MEASURE_VALUE_COLUMN,
    METRIC_COLUMNS,
    SUMMARY_RENDERERS,
    SUMMARY_SORT_MODES,
    generate_all_summary_plots,
    peptidoform_label_series,
    plot_multi_boxplot,
    plot_ridgeline,
    prepare_single_peptide_measures,
    prepare_summary_data,
)

__all__ = [
    # plot_worker — picklable work-unit + worker that renders one ComparisonTask
    "ComparisonTask",
    "process_single_task",
    # mirror — single/mirror spectrum figure builder
    "create_spectrum_plot",
    # summary — metric columns + per-renderer + end-to-end summary pipeline
    "METRIC_COLUMNS",
    "SUMMARY_RENDERERS",
    "SUMMARY_SORT_MODES",
    "MEASURE_VALUE_COLUMN",
    "generate_all_summary_plots",
    "peptidoform_label_series",
    "plot_multi_boxplot",
    "plot_ridgeline",
    "prepare_single_peptide_measures",
    "prepare_summary_data",
]
