"""Plot generation: Plotly + Matplotlib mirror plots, summary plots, parallel workers.

Re-exports the public API of the submodules so callers can keep using
``from mspecviz.plotting import ...`` without caring which file a name
lives in.
"""
from mspecviz.plotting.mirror_matplotlib import (
    create_mirror_plot,
    create_single_spectrum_plot,
    figure_to_bytes,
)
from mspecviz.plotting.mirror_plotly import (
    create_plotly_mirror_plot,
    create_plotly_single_plot,
)
from mspecviz.plotting.summary_dual_axis import plot_classic_dual_axis
from mspecviz.plotting.summary_orchestration import (
    DEFAULT_METRIC,
    METRIC_COLUMNS,
    PLOT_TYPES,
    generate_all_summary_plots,
    prepare_summary_data,
)
from mspecviz.plotting.summary_ridgeline import plot_ridgeline
from mspecviz.plotting.worker_execution import (
    process_single_task,
    run_parallel_comparisons,
)
from mspecviz.plotting.worker_tasks import ComparisonResult, ComparisonTask

# Re-export ion colour palette from common (plotting is the public home)
from mspecviz.common.ion_utils import ION_COLORS, get_ion_color

__all__ = [
    # mirror plots
    "create_mirror_plot",
    "create_plotly_mirror_plot",
    "create_plotly_single_plot",
    "create_single_spectrum_plot",
    "figure_to_bytes",
    "ION_COLORS",
    "get_ion_color",
    # summary plots
    "DEFAULT_METRIC",
    "METRIC_COLUMNS",
    "PLOT_TYPES",
    "generate_all_summary_plots",
    "plot_classic_dual_axis",
    "plot_ridgeline",
    "prepare_summary_data",
    # parallel workers
    "ComparisonResult",
    "ComparisonTask",
    "process_single_task",
    "run_parallel_comparisons",
]
