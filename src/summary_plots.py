"""
Summary plots for comparison metric data.

Generates per-peptidoform metric distribution visualizations:
  - Classic dual-axis (count bars + metric boxplots)
  - Ridgeline / Joy plot (KDE densities + count sidebar)

Can be called programmatically::

    from src.summary_plots import generate_all_summary_plots
    results = generate_all_summary_plots(source=df, output_dir="exports/run1")
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from src.common.text_utils import sanitize_filename
from src.peptide_search import format_inline_sequence

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METRIC_COLUMNS = [
    "R\u00b2 match",
    "R\u00b2 binned",
    "Cosine",
    "Entropy",
    "Dice",
    "Jaccard",
    "SharedRatio",
]

DEFAULT_METRIC = "R\u00b2 binned"

SUMMARY_COLORS = {
    "bar": "#E8945A",
    "bar_edge": "#D07030",
    "box": "#2C3E50",
    "box_fill": "#BDC3C7",
    "scatter": "#2C3E50",
}

PLOT_TYPES = ["dual_axis", "ridgeline"]

MAX_SCATTER_PER_GROUP = 200


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def prepare_summary_data(
    source: Union[str, Path, pd.DataFrame],
    metric: str = DEFAULT_METRIC,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and aggregate comparison data into per-peptidoform statistics.

    Returns (raw_df, agg_df).
    """
    if isinstance(source, (str, Path)):
        raw_df = pd.read_csv(source, low_memory=False)
    else:
        raw_df = source.copy()

    if metric not in raw_df.columns:
        raise ValueError(
            f"Metric '{metric}' not found. Available: "
            f"{[c for c in raw_df.columns if c in METRIC_COLUMNS]}"
        )

    # Build peptidoform label
    mod = raw_df["Modification"].fillna("").astype(str)
    raw_df["peptidoform"] = (
        raw_df.apply(
            lambda row: format_inline_sequence(row["Peptide"], mod[row.name]),
            axis=1,
        )
        + "/"
        + raw_df["Charge"].astype(str)
    )

    raw_df = raw_df.dropna(subset=[metric]).copy()

    grouped = raw_df.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg(["count", "median", "mean", "std"])
    agg.columns = ["count", "median", "mean", "std"]
    agg["q25"] = grouped[metric].quantile(0.25)
    agg["q75"] = grouped[metric].quantile(0.75)
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])

    charge_mode = grouped["Charge"].agg(lambda x: x.mode().iloc[0] if len(x) > 0 else 0)
    agg["charge"] = charge_mode

    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_df, agg


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    return max(minimum, per_row * n)


def _style_ax(ax, remove_top=True, remove_right=True):
    if remove_top:
        ax.spines["top"].set_visible(False)
    if remove_right:
        ax.spines["right"].set_visible(False)


def _subsample_group(group_df, metric, max_n=MAX_SCATTER_PER_GROUP, seed=42):
    vals = group_df[metric].values
    if len(vals) <= max_n:
        return vals
    rng = np.random.RandomState(seed)
    return rng.choice(vals, size=max_n, replace=False)


def _metric_label(metric: str) -> str:
    labels = {
        "R\u00b2 match": "R\u00b2 (match)",
        "R\u00b2 binned": "R\u00b2 (binned)",
        "Cosine": "Cosine similarity",
        "Entropy": "Spectral entropy similarity",
        "Dice": "Dice coefficient",
        "Jaccard": "Jaccard index",
        "SharedRatio": "Shared peak ratio",
    }
    return labels.get(metric, metric)


# ---------------------------------------------------------------------------
# Plot 1: Classic Dual-Axis
# ---------------------------------------------------------------------------

def plot_classic_dual_axis(
    raw_df: pd.DataFrame,
    agg_df: pd.DataFrame,
    metric: str = DEFAULT_METRIC,
    figsize: Tuple[float, float] = None,
    dpi: int = 300,
    title: str = None,
) -> Figure:
    """Horizontal count bars + metric boxplots."""
    agg_df = agg_df.sort_values("count", ascending=False).reset_index(drop=True)

    n = len(agg_df)
    if figsize is None:
        figsize = (14, _auto_height(n))

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    y_pos = np.arange(n)
    labels = agg_df["peptidoform"].values

    ax.barh(y_pos, agg_df["count"].values, height=0.65, color=SUMMARY_COLORS["bar"],
            alpha=0.55, edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5, zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax, remove_top=False)

    ax2 = ax.twiny()
    for i, pf in enumerate(labels):
        vals = _subsample_group(raw_df[raw_df["peptidoform"] == pf], metric)
        if len(vals) < 2:
            ax2.scatter(vals, [i] * len(vals), s=15, color=SUMMARY_COLORS["scatter"],
                        alpha=0.8, zorder=4)
            continue

        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        wlo = max(vals.min(), q1 - 1.5 * iqr)
        whi = min(vals.max(), q3 + 1.5 * iqr)
        outliers = vals[(vals < wlo) | (vals > whi)]

        ax2.bxp(
            [dict(med=med, q1=q1, q3=q3, whislo=wlo, whishi=whi, fliers=outliers)],
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


# ---------------------------------------------------------------------------
# Plot 2: Ridgeline
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Plot 3: Multi-Metric Heatmap
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

_PLOT_REGISTRY = {
    "dual_axis": plot_classic_dual_axis,
    "ridgeline": plot_ridgeline,
}


def generate_all_summary_plots(
    source: Union[str, Path, pd.DataFrame],
    output_dir: Union[str, Path],
    metric: str = DEFAULT_METRIC,
    plot_types: Optional[List[str]] = None,
    formats: Optional[List[str]] = None,
    dpi: int = 300,
    title_prefix: str = "",
) -> Dict[str, Dict[str, str]]:
    """Generate all (or selected) summary plots and save them.

    Returns ``{plot_type: {format: file_path}}``.
    """
    if plot_types is None:
        plot_types = list(PLOT_TYPES)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    raw_df, agg_df = prepare_summary_data(source, metric=metric)
    logger.info("Summary data: %d peptidoforms, %d rows, metric=%s",
                len(agg_df), len(raw_df), metric)

    results: Dict[str, Dict[str, str]] = {}

    for pt in plot_types:
        func = _PLOT_REGISTRY.get(pt)
        if func is None:
            logger.warning("Unknown plot type: %s (skipping)", pt)
            continue

        plot_title = (f"{title_prefix} -- {_metric_label(metric)}".strip(" -")
                      if title_prefix else None)

        try:
            fig = func(raw_df, agg_df, metric=metric, dpi=dpi, title=plot_title)
        except Exception as exc:
            logger.error("Failed to create %s plot: %s", pt, exc)
            continue

        paths: Dict[str, str] = {}
        filename = f"ridgeline_r2_binned" if pt == "ridgeline" and metric == DEFAULT_METRIC else pt
        for fmt in formats:
            out_path = summary_dir / f"{filename}.{fmt}"
            try:
                fig.savefig(str(out_path), format=fmt, dpi=dpi,
                            bbox_inches="tight", facecolor="white")
                paths[fmt] = str(out_path)
            except Exception as exc:
                logger.error("Failed to save %s.%s: %s", pt, fmt, exc)

        plt.close(fig)
        if paths:
            results[pt] = paths

    return results
