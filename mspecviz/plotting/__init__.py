"""Plot generation: Plotly + Matplotlib mirror plots, summary plots, parallel workers.

Re-exports the public API of the submodules so callers can use
``from mspecviz.plotting import ...`` without caring which file a name lives in.
"""
from mspecviz.plotting._summary_shared import METRIC_COLUMNS
from mspecviz.plotting.mirror_matplotlib import (
    create_mirror_plot,
    create_single_spectrum_plot,
    figure_to_bytes,
)
from mspecviz.plotting.mirror_plotly import (
    create_plotly_mirror_plot,
    create_plotly_single_plot,
)
from mspecviz.plotting.plot_task_dataclasses import ComparisonTask
from mspecviz.plotting.plot_worker import process_single_task
from mspecviz.plotting.summary_dual_axis import plot_classic_dual_axis
from mspecviz.plotting.summary_pipeline import generate_all_summary_plots, prepare_summary_data
from mspecviz.plotting.summary_ridgeline import plot_ridgeline

__all__ = [
    "ComparisonTask",
    "METRIC_COLUMNS",
    "create_mirror_plot",
    "create_plotly_mirror_plot",
    "create_plotly_single_plot",
    "create_single_spectrum_plot",
    "figure_to_bytes",
    "generate_all_summary_plots",
    "plot_classic_dual_axis",
    "plot_ridgeline",
    "prepare_summary_data",
    "process_single_task",
]
