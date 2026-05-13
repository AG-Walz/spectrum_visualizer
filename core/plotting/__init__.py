"""Plot generation: mirror plots (pyopenms-viz wrapper), summary plots, parallel workers.

Re-exports the public API of the submodules so callers can use
``from core.plotting import ...`` without caring which file a name lives in.
"""
from core.plotting._summary_shared import METRIC_COLUMNS
from core.plotting.mirror import create_spectrum_plot
from core.plotting.plot_task_dataclasses import ComparisonTask
from core.plotting.plot_worker import process_single_task
from core.plotting.summary_dual_axis import plot_classic_dual_axis
from core.plotting.summary_pipeline import generate_all_summary_plots, prepare_summary_data
from core.plotting.summary_ridgeline import plot_ridgeline

__all__ = [
    "ComparisonTask",
    "METRIC_COLUMNS",
    "create_spectrum_plot",
    "generate_all_summary_plots",
    "plot_classic_dual_axis",
    "plot_ridgeline",
    "prepare_summary_data",
    "process_single_task",
]
