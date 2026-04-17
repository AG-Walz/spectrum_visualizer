"""Ridgeline / Joy plot: stacked KDE densities + count sidebar."""
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from mspecviz.plotting._summary_shared import (
    DEFAULT_METRIC,
    SUMMARY_COLORS,
    _auto_height,
    _metric_label,
    _style_ax,
)


def plot_ridgeline(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULT_METRIC,
    figure_size: Tuple[float, float] = None,
    export_resolution: int = 300,
    title: str = None,
) -> Figure:
    """Stacked KDE density plots with count sidebar."""
    n = len(aggregated_frame)
    overlap = 0.65
    if figure_size is None:
        figure_size = (14, _auto_height(n, per_row=0.5))

    fig, (ax_ridge, ax_bar) = plt.subplots(
        1, 2, figsize=figure_size, dpi=export_resolution,
        gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05}, sharey=True,
    )

    labels = aggregated_frame["peptidoform"].values
    counts = aggregated_frame["count"].values
    x_grid = np.linspace(0, 1, 300)
    medians = aggregated_frame["median"].values
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0.0, vmax=1.0)
    y_offsets = np.arange(n) * overlap

    for i in range(n - 1, -1, -1):
        peptidoform = labels[i]
        values = raw_frame.loc[raw_frame["peptidoform"] == peptidoform, metric].dropna().values
        y_base = y_offsets[i]
        color = cmap(norm(medians[i]))

        if len(values) < 5:
            counts_h, edges = np.histogram(values, bins=20, range=(0, 1))
            density = counts_h / (counts_h.max() + 1e-9) * 0.4
            centers = (edges[:-1] + edges[1:]) / 2
            ax_ridge.bar(centers, density, bottom=y_base, width=edges[1] - edges[0],
                         color=color, alpha=0.7, edgecolor="white", linewidth=0.3)
            continue

        try:
            kde = gaussian_kde(values, bw_method="scott")
            density = kde(x_grid)
            density = density / density.max() * 0.55
        except np.linalg.LinAlgError:
            continue

        ax_ridge.fill_between(x_grid, y_base, y_base + density,
                              color=color, alpha=0.75, zorder=n - i)
        median = medians[i]
        kde_at_median = density[np.argmin(np.abs(x_grid - median))]
        ax_ridge.vlines(median, y_base, y_base + kde_at_median, color="#2C3E50",
                        linewidth=0.8, linestyle="--", alpha=0.6, zorder=n - i + 2)

    ax_ridge.set_yticks(y_offsets)
    ax_ridge.set_yticklabels(labels, fontsize=8)
    ax_ridge.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ax_ridge.set_xlim(-0.02, 1.02)
    ax_ridge.spines["top"].set_visible(False)
    ax_ridge.spines["right"].set_visible(False)
    ax_ridge.spines["left"].set_visible(False)
    ax_ridge.tick_params(axis="y", length=0)

    bar_y = y_offsets + overlap * 0.25
    ax_bar.barh(bar_y, counts, height=overlap * 0.7,
                color=SUMMARY_COLORS["bar"], alpha=0.7,
                edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5)
    ax_bar.set_xlabel("Count", fontsize=10, fontweight="bold")
    ax_bar.tick_params(axis="y", left=False, labelleft=False)
    _style_ax(ax_bar)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    if title:
        fig.subplots_adjust(left=0.20, top=0.92, bottom=0.14, wspace=0.05)
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    else:
        fig.subplots_adjust(left=0.20, top=0.97, bottom=0.14, wspace=0.05)
    cbar_ax = fig.add_axes([0.25, 0.02, 0.35, 0.015])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(f"Median {_metric_label(metric)}", fontsize=9)
    return fig
