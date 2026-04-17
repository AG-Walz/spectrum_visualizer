"""Classic dual-axis summary plot: horizontal count bars + metric boxplots."""
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mspecviz.plotting._summary_shared import (
    DEFAULT_METRIC,
    SUMMARY_COLORS,
    _auto_height,
    _metric_label,
    _style_ax,
    _subsample_group,
)


def plot_classic_dual_axis(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULT_METRIC,
    figure_size: Tuple[float, float] = None,
    export_resolution: int = 300,
    title: str = None,
) -> Figure:
    """Horizontal count bars + metric boxplots."""
    aggregated_frame = aggregated_frame.sort_values("count", ascending=False).reset_index(drop=True)

    n = len(aggregated_frame)
    if figure_size is None:
        figure_size = (14, _auto_height(n))

    fig, ax = plt.subplots(figsize=figure_size, dpi=export_resolution)
    y_pos = np.arange(n)
    labels = aggregated_frame["peptidoform"].values

    ax.barh(y_pos, aggregated_frame["count"].values, height=0.65, color=SUMMARY_COLORS["bar"],
            alpha=0.55, edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5, zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax, remove_top=False)

    ax2 = ax.twiny()
    for i, peptidoform in enumerate(labels):
        values = _subsample_group(raw_frame[raw_frame["peptidoform"] == peptidoform], metric)
        if len(values) < 2:
            ax2.scatter(values, [i] * len(values), s=15, color=SUMMARY_COLORS["scatter"],
                        alpha=0.8, zorder=4)
            continue

        q1, median, q3 = np.percentile(values, [25, 50, 75])
        iqr = q3 - q1
        whisker_low = max(values.min(), q1 - 1.5 * iqr)
        whisker_high = min(values.max(), q3 + 1.5 * iqr)
        outliers = values[(values < whisker_low) | (values > whisker_high)]

        ax2.bxp(
            [dict(med=median, q1=q1, q3=q3, whislo=whisker_low, whishi=whisker_high, fliers=outliers)],
            positions=[i], vert=False, widths=0.35, patch_artist=True, showfliers=True,
            boxprops=dict(facecolor=SUMMARY_COLORS["box_fill"], edgecolor=SUMMARY_COLORS["box"],
                          linewidth=0.8, alpha=0.85),
            whiskerprops=dict(color=SUMMARY_COLORS["box"], linewidth=0.8),
            capprops=dict(color=SUMMARY_COLORS["box"], linewidth=0.8),
            medianprops=dict(color=SUMMARY_COLORS["box"], linewidth=1.5),
            flierprops=dict(marker=".", markersize=3, markerfacecolor=SUMMARY_COLORS["scatter"],
                            markeredgecolor="none", alpha=0.4),
        )

    if not title:
        ax2.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ax2.set_xlim(-0.05, 1.05)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)

    legend_elements = [
        Patch(facecolor=SUMMARY_COLORS["bar"], alpha=0.55, edgecolor=SUMMARY_COLORS["bar_edge"],
              label="Peptide Count"),
        Line2D([0], [0], color=SUMMARY_COLORS["box"], marker="s", markersize=6,
               markerfacecolor=SUMMARY_COLORS["box_fill"], linestyle="-", linewidth=0.8,
               label=f"{_metric_label(metric)}"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0.0, -0.04),
              fontsize=9, framealpha=0.9, ncol=2)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
        fig.subplots_adjust(left=0.22, top=0.92)
    else:
        fig.subplots_adjust(left=0.22, top=0.96)
    return fig
