"""Plotting public API: mirror plots, summary plots, parallel render worker 
(imported by core.export, cli.main, streamlit_app.views.explore).
"""
from core.plotting._summary_shared import METRIC_COLUMNS
from core.plotting.mirror import create_spectrum_plot
from core.plotting.plot_task_dataclasses import ComparisonTask
from core.plotting.plot_worker import process_single_task
from core.plotting.summary_dual_axis import plot_classic_dual_axis
from core.plotting.summary_pipeline import generate_all_summary_plots, prepare_summary_data
from core.plotting.summary_ridgeline import plot_ridgeline

__all__ = [
    # plot_task_dataclasses — picklable work-unit passed to the parallel render worker
    "ComparisonTask",
    # _summary_shared — metric columns the summary plots draw from
    "METRIC_COLUMNS",
    # mirror — single/mirror spectrum figure builder
    "create_spectrum_plot",
    # summary_pipeline — end-to-end summary-figure generation and its input preparation
    "generate_all_summary_plots",
    # summary_dual_axis / summary_ridgeline — individual summary-figure renderers
    "plot_classic_dual_axis",
    "plot_ridgeline",
    "prepare_summary_data",
    # plot_worker — renders one ComparisonTask to disk (runs inside a worker process)
    "process_single_task",
]
