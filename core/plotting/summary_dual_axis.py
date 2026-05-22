"""Dual-axis summary plot: peptide-count bars twinned with per-metric boxplots.
Role — a leaf-level renderer in the summary-plot layer; produces one Matplotlib
figure and holds no aggregation logic of its own. What it produces — a single
Figure showing, for every peptidoform, a horizontal count bar on the bottom
axis overlaid with a similarity-metric boxplot on a twinned top axis, rows sorted
most-prevalent first. Call relationships — plot_classic_dual_axis is dispatched
by core.plotting.summary_pipeline.generate_all_summary_plots and re-exported
through core.plotting.__init__; it consumes the aggregated/raw frames built by
summary_pipeline.prepare_summary_data and the _auto_height / _metric_label
helpers from core.plotting._summary_shared.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from core.defaults import DEFAULTS  # noqa: E402
from core.plotting._summary_shared import _auto_height, _metric_label  # noqa: E402

# Palette for the two overlaid series — warm orange for the count bars, slate
# grey for the metric boxplots — chosen so the twinned axes read distinctly.
_BAR_COLOR = "#E8945A"  # Count-bar fill
_BAR_EDGE = "#D07030"   # Count-bar outline
_BOX_COLOR = "#2C3E50"  # Boxplot lines (median, whiskers, caps, outliers)
_BOX_FILL = "#BDC3C7"   # Boxplot box fill


def plot_classic_dual_axis(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render the dual-axis summary figure: count bars twinned with metric boxplots.
    Draws one horizontal row per peptidoform — a peptide-count bar on the bottom
    x axis overlaid with a boxplot of the chosen similarity metric on a
    twinned top x axis. raw_frame supplies the per-psm metric values the
    boxplots summarise; aggregated_frame supplies the per-peptidoform counts.
    Rows are ordered most-prevalent first (descending count).
    Returns the assembled Matplotlib Figure; the caller is responsible for
    saving and closing it. figure_size defaults to a fixed 14-inch width with
    a height that scales with the peptidoform count when not supplied; an explicit
    title suppresses the metric axis label and adds a figure suptitle instead.
    """
    # Order rows so the most-counted peptidoform sits at the top of the chart.
    aggregated_frame = aggregated_frame.sort_values("count", ascending=False).reset_index(drop=True)
    variant_count = len(aggregated_frame)
    # Auto-size height to the row count when the caller gave no explicit size.
    if figure_size is None:
        figure_size = (14, _auto_height(variant_count))

    fig, ax = plt.subplots(figsize=figure_size, dpi=export_resolution)
    y_pos = np.arange(variant_count)
    labels = aggregated_frame["peptidoform"].tolist()

    # Count bars on the bottom axis (one bar per peptidoform)
    ax.barh(y_pos, aggregated_frame["count"].values, height=0.65,
            color=_BAR_COLOR, alpha=0.55, edgecolor=_BAR_EDGE, linewidth=0.5, zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.spines["right"].set_visible(False)

    # Metric boxplots on the twinned top axis (seaborn draws IQR / whiskers / outliers)
    ax2 = ax.twiny()
    sns.boxplot(
        data=raw_frame, y="peptidoform", x=metric, order=labels, ax=ax2,
        color=_BOX_FILL, width=0.35, linewidth=0.8, fliersize=3,
        boxprops=dict(facecolor=_BOX_FILL, edgecolor=_BOX_COLOR, alpha=0.85),
        medianprops=dict(color=_BOX_COLOR, linewidth=1.5),
        whiskerprops=dict(color=_BOX_COLOR, linewidth=0.8),
        capprops=dict(color=_BOX_COLOR, linewidth=0.8),
        flierprops=dict(marker=".", markerfacecolor=_BOX_COLOR,
                        markeredgecolor="none", alpha=0.4),
    )
    if not title:
        ax2.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylabel("")
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Y-axis ticks label every row with its peptidoform; the legend is built from
    # hand-made proxy handles (a Patch + a Line2D) since the two series live on
    # separate twinned axes and have no shared auto-legend.
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.legend(handles=[
        Patch(facecolor=_BAR_COLOR, alpha=0.55, edgecolor=_BAR_EDGE, label="Peptide Count"),
        Line2D([0], [0], color=_BOX_COLOR, marker="s", markersize=6,
               markerfacecolor=_BOX_FILL, linestyle="-", linewidth=0.8,
               label=_metric_label(metric)),
    ], loc="upper left", bbox_to_anchor=(0.0, -0.04), fontsize=9, framealpha=0.9, ncol=2)

    # Reserve top margin for a suptitle when one is given, else keep it tight.
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.22, top=0.92 if title else 0.96)
    return fig
