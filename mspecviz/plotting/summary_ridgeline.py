"""Ridgeline / Joy plot: stacked KDE densities + count sidebar."""
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from mspecviz.plotting.summary_orchestration import (
    DEFAULT_METRIC,
    SUMMARY_COLORS,
    _auto_height,
    _metric_label,
    _style_ax,
)


def plot_ridgeline(
    raw_df: pd.DataFrame,
    agg_df: pd.DataFrame,
    metric: str = DEFAULT_METRIC,
    figsize: Tuple[float, float] = None,
    dpi: int = 300,
    title: str = None,
) -> Figure:
    """Stacked KDE density plots with count sidebar."""
    n = len(agg_df)
    overlap = 0.65
    if figsize is None:
        figsize = (14, _auto_height(n, per_row=0.5))

    fig, (ax_ridge, ax_bar) = plt.subplots(
        1, 2, figsize=figsize, dpi=dpi,
        gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05}, sharey=True,
    )

    labels = agg_df["peptidoform"].values
    counts = agg_df["count"].values
    x_grid = np.linspace(0, 1, 300)
    medians = agg_df["median"].values
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0.0, vmax=1.0)
    y_offsets = np.arange(n) * overlap

    for i in range(n - 1, -1, -1):
        pf = labels[i]
        vals = raw_df.loc[raw_df["peptidoform"] == pf, metric].dropna().values
        y_base = y_offsets[i]
        color = cmap(norm(medians[i]))

        if len(vals) < 5:
            counts_h, edges = np.histogram(vals, bins=20, range=(0, 1))
            density = counts_h / (counts_h.max() + 1e-9) * 0.4
            centers = (edges[:-1] + edges[1:]) / 2
            ax_ridge.bar(centers, density, bottom=y_base, width=edges[1] - edges[0],
                         color=color, alpha=0.7, edgecolor="white", linewidth=0.3)
            continue

        try:
            kde = gaussian_kde(vals, bw_method="scott")
            density = kde(x_grid)
            density = density / density.max() * 0.55
        except Exception:
            continue

        ax_ridge.fill_between(x_grid, y_base, y_base + density,
                              color=color, alpha=0.75, zorder=n - i)
        med = medians[i]
        kde_at_med = density[np.argmin(np.abs(x_grid - med))]
        ax_ridge.vlines(med, y_base, y_base + kde_at_med, color="#2C3E50",
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
