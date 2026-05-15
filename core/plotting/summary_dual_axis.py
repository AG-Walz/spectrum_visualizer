"""Classic dual-axis summary plot: horizontal count bars + metric boxplots.

Bottom axis is hand-drawn (matplotlib ``barh`` for counts); the top twin axis
delegates to ``seaborn.boxplot`` for the per-variant IQR rendering.
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

from core.common.defaults import DEFAULTS  # noqa: E402
from core.plotting._summary_shared import _auto_height, _metric_label  # noqa: E402

_BAR_COLOR = "#E8945A"
_BAR_EDGE = "#D07030"
_BOX_COLOR = "#2C3E50"
_BOX_FILL = "#BDC3C7"


def plot_classic_dual_axis(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Horizontal count bars + metric boxplots."""
    aggregated_frame = aggregated_frame.sort_values("count", ascending=False).reset_index(drop=True)
    variant_count = len(aggregated_frame)
    if figure_size is None:
        figure_size = (14, _auto_height(variant_count))

    fig, ax = plt.subplots(figsize=figure_size, dpi=export_resolution)
    y_pos = np.arange(variant_count)
    labels = aggregated_frame["peptidoform"].tolist()

    # Count bars (bottom axis).
    ax.barh(y_pos, aggregated_frame["count"].values, height=0.65,
            color=_BAR_COLOR, alpha=0.55, edgecolor=_BAR_EDGE, linewidth=0.5, zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.spines["right"].set_visible(False)

    # Metric boxplots (twin top axis) — seaborn handles IQR / whiskers / outliers.
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

    # Y-axis ticks + legend.
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.legend(handles=[
        Patch(facecolor=_BAR_COLOR, alpha=0.55, edgecolor=_BAR_EDGE, label="Peptide Count"),
        Line2D([0], [0], color=_BOX_COLOR, marker="s", markersize=6,
               markerfacecolor=_BOX_FILL, linestyle="-", linewidth=0.8,
               label=_metric_label(metric)),
    ], loc="upper left", bbox_to_anchor=(0.0, -0.04), fontsize=9, framealpha=0.9, ncol=2)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.22, top=0.92 if title else 0.96)
    return fig
