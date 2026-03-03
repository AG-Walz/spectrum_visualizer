"""
Summary plots for MSpecViz comparison_summary.csv data.

Generates per-peptidoform metric distribution visualisations with peptide
count overlays. Three plot styles are available: classic dual-axis,
ridgeline, and heatmap.

Additionally provides per-peptide detail plots that show R² (binned) and
R² (match) for each query × reference pair within a single peptidoform.

Can be run standalone::

    python -m src.summary_plots path/to/comparison_summary.csv
    python -m src.summary_plots path/to/csv -m "Cosine" --plots dual_axis heatmap
    python -m src.summary_plots path/to/csv --peptide-detail

Or called programmatically::

    from src.summary_plots import generate_all_summary_plots
    results = generate_all_summary_plots(source=df, output_dir="exports/run1")

    from src.summary_plots import generate_peptide_detail_plots
    results = generate_peptide_detail_plots(source=df, output_dir="exports/run1")
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METRIC_COLUMNS = [
    "R² match",
    "R² binned",
    "Cosine",
    "Entropy",
    "Dice",
    "Jaccard",
    "SharedRatio",
]

DEFAULT_METRIC = "R² binned"

SUMMARY_COLORS = {
    "bar": "#E8945A",        # warm orange for count bars
    "bar_edge": "#D07030",
    "box": "#2C3E50",        # dark slate for boxplots
    "box_fill": "#BDC3C7",   # light gray fill
    "scatter": "#2C3E50",    # dark for scatter dots
    "heatmap_cmap": "RdYlGn",
}

PLOT_TYPES = [
    "dual_axis",
    "ridgeline",
    "heatmap",
]

MAX_SCATTER_PER_GROUP = 200  # subsample scatter to avoid slow renders

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def prepare_summary_data(
    source: Union[str, Path, pd.DataFrame],
    metric: str = DEFAULT_METRIC,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and aggregate comparison data into per-peptidoform statistics.

    Parameters
    ----------
    source : str, Path, or DataFrame
        Path to ``comparison_summary.csv`` or an already-loaded DataFrame.
    metric : str
        Which metric column to focus on (must exist in data).

    Returns
    -------
    (raw_df, agg_df)
        *raw_df* has an added ``peptidoform`` column.
        *agg_df* has columns: ``peptidoform, count, median, mean, std,
        q25, q75, ci_lo, ci_hi, charge`` sorted by count descending.

    Raises
    ------
    ValueError
        If *metric* is not found in the data columns.
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
        raw_df["Peptide"]
        + mod.apply(lambda x: f"[{x}]" if x else "")
        + "/"
        + raw_df["Charge"].astype(str)
    )

    # Drop rows where metric is NaN
    raw_df = raw_df.dropna(subset=[metric]).copy()

    # Aggregate
    grouped = raw_df.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg(["count", "median", "mean", "std"])
    agg.columns = ["count", "median", "mean", "std"]
    agg["q25"] = grouped[metric].quantile(0.25)
    agg["q75"] = grouped[metric].quantile(0.75)
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])

    # Predominant charge per peptidoform
    charge_mode = grouped["Charge"].agg(lambda x: x.mode().iloc[0] if len(x) > 0 else 0)
    agg["charge"] = charge_mode

    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_df, agg


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    """Compute figure height from number of peptidoforms."""
    return max(minimum, per_row * n)


def _style_ax(ax, remove_top=True, remove_right=True):
    """Apply house style: remove spines, set white background."""
    if remove_top:
        ax.spines["top"].set_visible(False)
    if remove_right:
        ax.spines["right"].set_visible(False)


def _subsample_group(group_df, metric, max_n=MAX_SCATTER_PER_GROUP, seed=42):
    """Subsample metric values to at most *max_n* for scatter plots."""
    vals = group_df[metric].values
    if len(vals) <= max_n:
        return vals
    rng = np.random.RandomState(seed)
    return rng.choice(vals, size=max_n, replace=False)


def _metric_label(metric: str) -> str:
    """Human-readable axis label for a metric."""
    labels = {
        "R² match": "R² (match)",
        "R² binned": "R² (binned)",
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
    """Horizontal count bars (bottom) + metric boxplots (top) — reference style."""

    # Highest count at top → sort descending for this plot
    agg_df = agg_df.sort_values("count", ascending=False).reset_index(drop=True)

    n = len(agg_df)
    if figsize is None:
        figsize = (14, _auto_height(n))

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    y_pos = np.arange(n)
    labels = agg_df["peptidoform"].values

    # --- Count bars (bottom axis) ---
    ax.barh(y_pos, agg_df["count"].values, height=0.65, color=SUMMARY_COLORS["bar"],
            alpha=0.55, edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5,
            label="Peptide count", zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax, remove_top=False)

    # --- Metric boxplots (top axis) ---
    ax2 = ax.twiny()
    for i, pf in enumerate(labels):
        vals = _subsample_group(
            raw_df[raw_df["peptidoform"] == pf], metric
        )
        if len(vals) < 2:
            ax2.scatter(vals, [i] * len(vals), s=15, color=SUMMARY_COLORS["scatter"],
                        alpha=0.8, zorder=4)
            continue

        # Manual boxplot stats
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        wlo = max(vals.min(), q1 - 1.5 * iqr)
        whi = min(vals.max(), q3 + 1.5 * iqr)
        outliers = vals[(vals < wlo) | (vals > whi)]

        ax2.bxp(
            [dict(med=med, q1=q1, q3=q3, whislo=wlo, whishi=whi, fliers=outliers)],
            positions=[i], vert=False, widths=0.35, patch_artist=True,
            showfliers=True,
            boxprops=dict(facecolor=SUMMARY_COLORS["box_fill"], edgecolor=SUMMARY_COLORS["box"],
                          linewidth=0.8, alpha=0.85),
            whiskerprops=dict(color=SUMMARY_COLORS["box"], linewidth=0.8),
            capprops=dict(color=SUMMARY_COLORS["box"], linewidth=0.8),
            medianprops=dict(color=SUMMARY_COLORS["box"], linewidth=1.5),
            flierprops=dict(marker=".", markersize=3, markerfacecolor=SUMMARY_COLORS["scatter"],
                            markeredgecolor="none", alpha=0.4),
        )

    # Remove redundant top-axis label when title already contains the metric
    if not title:
        ax2.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ax2.set_xlim(-0.05, 1.05)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Re-set y-ticks AFTER bxp (twiny shares y-axis and bxp resets labels)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)

    # Legend — place outside the plot area to avoid overlap
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
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
        fig.subplots_adjust(left=0.22, top=0.96)
    else:
        fig.subplots_adjust(left=0.22, top=0.96)
    return fig


# ---------------------------------------------------------------------------
# Plot 2: Ridgeline / Joy Plot (with count sidebar barplot)
# ---------------------------------------------------------------------------

def plot_ridgeline(
    raw_df: pd.DataFrame,
    agg_df: pd.DataFrame,
    metric: str = DEFAULT_METRIC,
    figsize: Tuple[float, float] = None,
    dpi: int = 300,
    title: str = None,
) -> Figure:
    """Stacked overlapping KDE density plots with count barplot sidebar.

    Sorted by count (highest to lowest). Count barplot on the right.
    """
    from scipy.stats import gaussian_kde

    n = len(agg_df)
    overlap = 0.65
    if figsize is None:
        figsize = (14, _auto_height(n, per_row=0.5))

    fig, (ax_ridge, ax_bar) = plt.subplots(
        1, 2, figsize=figsize, dpi=dpi,
        gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05},
        sharey=True,
    )

    labels = agg_df["peptidoform"].values
    counts = agg_df["count"].values
    x_grid = np.linspace(0, 1, 300)

    # Colormap based on median metric value
    medians = agg_df["median"].values
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0.0, vmax=1.0)

    y_offsets = np.arange(n) * overlap

    for i in range(n - 1, -1, -1):  # draw bottom to top
        pf = labels[i]
        vals = raw_df.loc[raw_df["peptidoform"] == pf, metric].dropna().values
        y_base = y_offsets[i]
        color = cmap(norm(medians[i]))

        if len(vals) < 5:
            # Histogram fallback
            counts_h, edges = np.histogram(vals, bins=20, range=(0, 1))
            density = counts_h / (counts_h.max() + 1e-9) * 0.4
            centers = (edges[:-1] + edges[1:]) / 2
            ax_ridge.bar(centers, density, bottom=y_base, width=edges[1] - edges[0],
                         color=color, alpha=0.7, edgecolor="white", linewidth=0.3)
            continue

        try:
            kde = gaussian_kde(vals, bw_method="scott")
            density = kde(x_grid)
            density = density / density.max() * 0.55  # scale to row height
        except Exception:
            continue

        ax_ridge.fill_between(x_grid, y_base, y_base + density,
                              color=color, alpha=0.75, zorder=n - i)

        # Median line
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

    # --- Count barplot sidebar ---
    # Shift bars up to align with the visual centre of each KDE ridge
    # (density fills upward from y_base, so the visual midpoint is above it)
    bar_y = y_offsets + overlap * 0.25
    ax_bar.barh(bar_y, counts, height=overlap * 0.7,
                color=SUMMARY_COLORS["bar"], alpha=0.7,
                edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5)
    ax_bar.set_xlabel("Count", fontsize=10, fontweight="bold")
    ax_bar.tick_params(axis="y", left=False, labelleft=False)
    _style_ax(ax_bar)

    # Colorbar — manually positioned well below the x-axis labels
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.subplots_adjust(left=0.20, top=0.97, bottom=0.14, wspace=0.05)
    cbar_ax = fig.add_axes([0.25, 0.02, 0.35, 0.015])  # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(f"Median {_metric_label(metric)}", fontsize=9)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)
    return fig


# ---------------------------------------------------------------------------
# Plot 3: Multi-Metric Heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(
    raw_df: pd.DataFrame,
    agg_df: pd.DataFrame,
    metric: str = DEFAULT_METRIC,  # ignored — shows all metrics
    metrics: List[str] = None,
    figsize: Tuple[float, float] = None,
    dpi: int = 300,
    title: str = None,
) -> Figure:
    """All metrics as heatmap columns with count sidebar."""

    # Highest count at top → sort descending for this plot
    agg_df = agg_df.sort_values("count", ascending=False).reset_index(drop=True)

    if metrics is None:
        metrics = [m for m in METRIC_COLUMNS if m in raw_df.columns]

    n = len(agg_df)
    n_metrics = len(metrics)
    if figsize is None:
        figsize = (4 + n_metrics * 1.2, _auto_height(n, per_row=0.42))

    fig, (ax_heat, ax_bar) = plt.subplots(
        1, 2, figsize=figsize, dpi=dpi,
        gridspec_kw={"width_ratios": [n_metrics, 1.5], "wspace": 0.08},
        sharey=True,
    )

    labels = agg_df["peptidoform"].values
    y_pos = np.arange(n)

    # Build median matrix
    grouped = raw_df.groupby("peptidoform", sort=False)
    matrix = np.zeros((n, n_metrics))
    for j, m in enumerate(metrics):
        medians_series = grouped[m].median()
        for i, pf in enumerate(labels):
            matrix[i, j] = medians_series.get(pf, np.nan)

    # Heatmap
    im = ax_heat.imshow(matrix, aspect="auto", cmap=SUMMARY_COLORS["heatmap_cmap"],
                        vmin=0, vmax=1, interpolation="nearest")

    # Cell annotations
    for i in range(n):
        for j in range(n_metrics):
            val = matrix[i, j]
            if not np.isnan(val):
                text_color = "white" if val < 0.3 or val > 0.85 else "black"
                ax_heat.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=6.5, color=text_color, fontweight="bold")

    ax_heat.set_xticks(range(n_metrics))
    ax_heat.set_xticklabels([_metric_label(m) for m in metrics], fontsize=8,
                            rotation=35, ha="right")
    ax_heat.set_yticks(y_pos)
    ax_heat.set_yticklabels(labels, fontsize=8)
    ax_heat.tick_params(axis="both", length=0)

    # Count sidebar — shared y-axis ensures alignment
    ax_bar.barh(y_pos, agg_df["count"].values, height=0.7,
                color=SUMMARY_COLORS["bar"], alpha=0.7,
                edgecolor=SUMMARY_COLORS["bar_edge"], linewidth=0.5)
    ax_bar.set_xlabel("Count", fontsize=10, fontweight="bold")
    ax_bar.tick_params(axis="y", left=False, labelleft=False)
    _style_ax(ax_bar)

    # Colorbar — manually positioned well below the x-axis labels
    fig.subplots_adjust(left=0.22, bottom=0.18, top=0.96, wspace=0.08)
    cbar_ax = fig.add_axes([0.25, 0.04, 0.35, 0.015])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Median metric value", fontsize=9)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    return fig


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

# Registry mapping plot type names to functions
_PLOT_REGISTRY = {
    "dual_axis": plot_classic_dual_axis,
    "ridgeline": plot_ridgeline,
    "heatmap": plot_heatmap,
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

    Parameters
    ----------
    source : str, Path, or DataFrame
        CSV path or pre-loaded DataFrame.
    output_dir : str or Path
        Parent export directory; a ``summary/`` subdirectory is created.
    metric : str
        Metric to visualise (default: "R² binned").
    plot_types : list of str, optional
        Which plot styles to generate. Default: all from :data:`PLOT_TYPES`.
    formats : list of str, optional
        Export formats, e.g. ``["png", "svg"]``. Default: ``["png"]``.
    dpi : int
        Resolution (default 300).
    title_prefix : str
        Optional prefix for figure titles.

    Returns
    -------
    dict
        ``{plot_type: {format: file_path}}``.
    """
    if plot_types is None:
        plot_types = list(PLOT_TYPES)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    raw_df, agg_df = prepare_summary_data(source, metric=metric)
    logger.info(
        "Summary data: %d peptidoforms, %d rows, metric=%s",
        len(agg_df), len(raw_df), metric,
    )

    results: Dict[str, Dict[str, str]] = {}

    for pt in plot_types:
        func = _PLOT_REGISTRY.get(pt)
        if func is None:
            logger.warning("Unknown plot type: %s (skipping)", pt)
            continue

        plot_title = f"{title_prefix} — {_metric_label(metric)}".strip(" —") if title_prefix else None

        try:
            fig = func(raw_df, agg_df, metric=metric, dpi=dpi, title=plot_title)
        except Exception as exc:
            logger.error("Failed to create %s plot: %s", pt, exc)
            continue

        paths: Dict[str, str] = {}
        # Ridgeline for default metric gets metric-specific filename
        filename = f"ridgeline_r2_binned" if pt == "ridgeline" and metric == "R² binned" else pt
        for fmt in formats:
            out_path = summary_dir / f"{filename}.{fmt}"
            try:
                fig.savefig(str(out_path), format=fmt, dpi=dpi,
                            bbox_inches="tight", facecolor="white")
                paths[fmt] = str(out_path)
                logger.info("Saved %s", out_path)
            except Exception as exc:
                logger.error("Failed to save %s.%s: %s", pt, fmt, exc)

        plt.close(fig)
        if paths:
            results[pt] = paths

    # --- Ridgeline for R² match (in addition to the selected metric) ---
    if "ridgeline" in plot_types and "R² match" in raw_df.columns and metric != "R² match":
        try:
            _, agg_match = prepare_summary_data(source, metric="R² match")
            match_title = f"{title_prefix} — {_metric_label('R² match')}".strip(" —") if title_prefix else None
            fig = plot_ridgeline(raw_df, agg_match, metric="R² match", dpi=dpi,
                                 title=match_title)
            paths_match: Dict[str, str] = {}
            for fmt in formats:
                out_path = summary_dir / f"ridgeline_r2_match.{fmt}"
                try:
                    fig.savefig(str(out_path), format=fmt, dpi=dpi,
                                bbox_inches="tight", facecolor="white")
                    paths_match[fmt] = str(out_path)
                    logger.info("Saved %s", out_path)
                except Exception as exc:
                    logger.error("Failed to save ridgeline_r2_match.%s: %s", fmt, exc)
            plt.close(fig)
            if paths_match:
                results["ridgeline_r2_match"] = paths_match
        except Exception as exc:
            logger.error("Failed to create ridgeline R² match plot: %s", exc)

    return results


# ===========================================================================
# Per-Peptide Detail Plots
# ===========================================================================

PEPTIDE_DETAIL_TYPES = ["dumbbell"]

METRIC_COLORS = {
    "R² binned": SUMMARY_COLORS["bar"],       # warm orange
    "R² match":  SUMMARY_COLORS["box"],        # dark slate
}

# ---------------------------------------------------------------------------
# Per-peptide data preparation
# ---------------------------------------------------------------------------

def _shorten_labels(
    names: List[str],
    max_len: int = 35,
) -> List[str]:
    """Two-tier label shortening: strip common tokens, then truncate+hash."""
    from src.spectrum_metrics import shorten_filenames
    from src.export import shorten_stem

    shortened = shorten_filenames(names)
    return [shorten_stem(s, max_len=max_len) for s in shortened]


def prepare_peptide_detail_data(
    source: Union[str, Path, pd.DataFrame],
    peptidoform: str = None,
) -> pd.DataFrame:
    """Load raw CSV, add peptidoform column, optionally filter to one peptidoform.

    Aggregates per (Query, Reference) pair to median R² binned and R² match.

    Returns
    -------
    DataFrame with columns: Query, Reference, R² binned, R² match, peptidoform
    """
    if isinstance(source, (str, Path)):
        raw_df = pd.read_csv(source, low_memory=False)
    else:
        raw_df = source.copy()

    # Build peptidoform label if missing
    if "peptidoform" not in raw_df.columns:
        mod = raw_df["Modification"].fillna("").astype(str)
        raw_df["peptidoform"] = (
            raw_df["Peptide"]
            + mod.apply(lambda x: f"[{x}]" if x else "")
            + "/"
            + raw_df["Charge"].astype(str)
        )

    if peptidoform is not None:
        raw_df = raw_df[raw_df["peptidoform"] == peptidoform]

    if raw_df.empty:
        return pd.DataFrame(columns=["Query", "Reference", "R² binned", "R² match", "peptidoform"])

    metrics_to_agg = {}
    if "R² binned" in raw_df.columns:
        metrics_to_agg["R² binned"] = ("R² binned", "median")
    if "R² match" in raw_df.columns:
        metrics_to_agg["R² match"] = ("R² match", "median")

    agg = (
        raw_df.groupby(["peptidoform", "Query", "Reference"], sort=False)
        .agg(**metrics_to_agg)
        .reset_index()
    )
    return agg


# ---------------------------------------------------------------------------
# Peptide detail plot: Dumbbell / Connected Dot Plot
# ---------------------------------------------------------------------------

def plot_dumbbell(
    detail_df: pd.DataFrame,
    peptidoform: str,
    figsize: Tuple[float, float] = None,
    dpi: int = 300,
    title: str = None,
) -> Figure:
    """Faceted by reference: two dots per query (R² binned + R² match)."""
    queries = sorted(detail_df["Query"].unique())
    refs = sorted(detail_df["Reference"].unique())
    n_q, n_r = len(queries), len(refs)

    short_q = _shorten_labels(queries, max_len=40)
    short_r = _shorten_labels(refs, max_len=40)

    lookup = detail_df.set_index(["Query", "Reference"])

    if figsize is None:
        w_per_panel = 5 if n_r <= 2 else 4.5
        w = max(8, w_per_panel * n_r + 2)
        # Compact rows: 0.25 per query, small padding for title/legend
        title_room = 1.0 if n_r > 1 else 0.6
        h = max(3.5, n_q * 0.25 + 1.5 + title_room)
        figsize = (min(w, 24), h)

    fig, axes = plt.subplots(
        1, n_r, sharey=True, figsize=figsize, dpi=dpi,
        gridspec_kw={"wspace": 0.12},
        squeeze=False,
    )
    axes = axes.ravel()

    col_binned = METRIC_COLORS["R² binned"]
    col_match = METRIC_COLORS["R² match"]

    for j, ref in enumerate(refs):
        ax = axes[j]
        for i, q in enumerate(queries):
            r2b, r2m = np.nan, np.nan
            if (q, ref) in lookup.index:
                row = lookup.loc[(q, ref)]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                r2b = row.get("R² binned", np.nan)
                r2m = row.get("R² match", np.nan)

            if np.isnan(r2b) and np.isnan(r2m):
                continue

            # Connecting line
            if not np.isnan(r2b) and not np.isnan(r2m):
                ax.plot([r2b, r2m], [i, i], color=SUMMARY_COLORS["box_fill"],
                        linewidth=1.5, zorder=1)

            s = 40 if n_q <= 20 else 25
            if not np.isnan(r2b):
                ax.scatter(r2b, i, color=col_binned, s=s, zorder=3, marker="o",
                           edgecolors="white", linewidths=0.3)
            if not np.isnan(r2m):
                ax.scatter(r2m, i, color=col_match, s=s, zorder=3, marker="D",
                           edgecolors="white", linewidths=0.3)

        ax.set_xlim(-0.05, 1.05)
        ax.axvline(x=0.8, color="#E0E0E0", linewidth=0.5, linestyle="--", zorder=0)
        ax.set_xlabel("R² value", fontsize=10)
        panel_title = f"vs {short_r[j]}" if n_r > 1 else short_r[j]
        ax.set_title(panel_title, fontsize=8, pad=8)
        _style_ax(ax)
        ax.grid(axis="x", alpha=0.15, linewidth=0.5)

    axes[0].set_yticks(range(n_q))
    axes[0].set_yticklabels(short_q, fontsize=7 if n_q > 20 else 8)

    # Legend inside the first axes at bottom-left (data clusters right)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=col_binned,
               markersize=8, label="R² (binned)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=col_match,
               markersize=8, label="R² (match)"),
    ]
    axes[0].legend(handles=legend_elements, loc="lower left", ncol=1,
                   fontsize=8, framealpha=0.9)

    fig.subplots_adjust(top=0.92, bottom=0.10)
    fig.suptitle(title or peptidoform, fontsize=13, fontweight="bold", y=1.03)
    return fig


# ---------------------------------------------------------------------------
# Per-peptide orchestrator
# ---------------------------------------------------------------------------

_PEPTIDE_PLOT_REGISTRY = {
    "dumbbell": plot_dumbbell,
}


def generate_peptide_detail_plots(
    source: Union[str, Path, pd.DataFrame],
    output_dir: Union[str, Path],
    peptidoforms: Optional[List[str]] = None,
    formats: Optional[List[str]] = None,
    dpi: int = 300,
    title_prefix: str = "",
) -> Dict[str, Dict[str, str]]:
    """Generate per-peptide dumbbell plots alongside mirror plots.

    Saves into ``{output_dir}/png/{Peptide}/`` and/or
    ``{output_dir}/svg/{Peptide}/`` — the same directories that already
    contain mirror plots.

    Parameters
    ----------
    source : str, Path, or DataFrame
        CSV path or pre-loaded DataFrame.
    output_dir : str or Path
        Root export directory (containing ``png/`` and/or ``svg/``).
    peptidoforms : list of str, optional
        Which peptidoforms to plot. Default: all found in data.
    formats : list of str, optional
        Export formats. Default: ``["png"]``.
    dpi : int
        Resolution (default 300).
    title_prefix : str
        Optional prefix for figure titles.

    Returns
    -------
    dict
        ``{peptidoform: {format: file_path}}``.
    """
    from src.export import sanitize_filename, shorten_stem

    if formats is None:
        formats = ["png"]

    # Load and build peptidoform column once
    full_df = prepare_peptide_detail_data(source)

    if peptidoforms is None:
        peptidoforms = sorted(full_df["peptidoform"].unique())

    logger.info(
        "Generating peptide dumbbell plots: %d peptidoforms",
        len(peptidoforms),
    )

    results: Dict[str, Dict[str, str]] = {}
    base = Path(output_dir)

    for pf in peptidoforms:
        pf_data = full_df[full_df["peptidoform"] == pf]
        if pf_data.empty:
            logger.warning("No data for peptidoform '%s' — skipping", pf)
            continue

        # Derive folder name from peptide (without charge) to match mirror plot dirs
        # Mirror plots use Peptide name as folder, e.g. svg/KEEPGAFVI/
        peptide_name = pf.rsplit("/", 1)[0]  # strip /charge
        folder_name = sanitize_filename(peptide_name)

        # Filename: PEPTIDE_charge_dumbbell.ext
        safe_pf = sanitize_filename(pf)
        filename = f"{shorten_stem(safe_pf, max_len=60)}_dumbbell"

        pf_title = f"{title_prefix} — {pf}".strip(" —") if title_prefix else pf

        try:
            fig = plot_dumbbell(pf_data, pf, dpi=dpi, title=pf_title)
        except Exception as exc:
            logger.error("Failed dumbbell for %s: %s", pf, exc)
            continue

        paths: Dict[str, str] = {}
        for fmt in formats:
            fmt_dir = base / fmt / folder_name
            fmt_dir.mkdir(parents=True, exist_ok=True)
            out_path = fmt_dir / f"{filename}.{fmt}"
            try:
                fig.savefig(str(out_path), format=fmt, dpi=dpi,
                            bbox_inches="tight", facecolor="white")
                paths[fmt] = str(out_path)
            except Exception as exc:
                logger.error("Failed to save %s/%s.%s: %s", pf, filename, fmt, exc)

        plt.close(fig)
        if paths:
            results[pf] = paths

    n_total = sum(len(v) for v in results.values())
    logger.info("Peptide dumbbell plots: %d files generated", n_total)
    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        prog="summary-plots",
        description="Generate summary plots from MSpecViz comparison_summary.csv",
    )
    parser.add_argument("csv", help="Path to comparison_summary.csv")
    parser.add_argument(
        "-o", "--output-dir", default=None,
        help="Output directory (default: same directory as CSV)",
    )
    parser.add_argument(
        "-m", "--metric", default=DEFAULT_METRIC,
        help=f"Metric to visualise (default: '{DEFAULT_METRIC}')",
    )
    parser.add_argument(
        "--plots", nargs="+", default=None, choices=PLOT_TYPES,
        help="Which plot types to generate (default: all)",
    )
    parser.add_argument(
        "--formats", nargs="+", default=None, choices=["png", "svg", "pdf"],
        help="Export formats (default: svg)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300, help="DPI (default: 300)",
    )
    parser.add_argument(
        "--title", default="", help="Title prefix for plots",
    )
    parser.add_argument(
        "--list-metrics", action="store_true",
        help="List available metrics in the CSV and exit",
    )
    parser.add_argument(
        "--peptide-detail", action="store_true",
        help="Generate per-peptide dumbbell plots into png/svg dirs",
    )
    parser.add_argument(
        "--peptidoforms", nargs="+", default=None,
        help="Specific peptidoforms to plot (default: all)",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"Error: {csv_path} not found", file=sys.stderr)
        sys.exit(1)

    if args.list_metrics:
        df = pd.read_csv(csv_path, nrows=0)
        available = [c for c in df.columns if c in METRIC_COLUMNS]
        print("Available metrics:")
        for m in available:
            print(f"  {m}")
        sys.exit(0)

    output_dir = args.output_dir or str(csv_path.parent)

    # Single format default for standalone CLI (user can override with --formats)
    formats = args.formats if args.formats else ["svg"]

    if not args.peptide_detail:
        print(f"Loading {csv_path.name} ...")
        results = generate_all_summary_plots(
            source=str(csv_path),
            output_dir=output_dir,
            metric=args.metric,
            plot_types=args.plots,
            formats=formats,
            dpi=args.dpi,
            title_prefix=args.title,
        )

        print(f"\nGenerated {sum(len(v) for v in results.values())} files:")
        for plot_type, paths in results.items():
            for fmt, path in paths.items():
                print(f"  {plot_type}.{fmt}: {path}")
    else:
        print(f"Loading {csv_path.name} for peptide detail plots ...")
        results = generate_peptide_detail_plots(
            source=str(csv_path),
            output_dir=output_dir,
            peptidoforms=args.peptidoforms,
            formats=formats,
            dpi=args.dpi,
            title_prefix=args.title,
        )

        n_total = sum(len(v) for v in results.values())
        print(f"\nGenerated {n_total} files for {len(results)} peptidoforms:")
        for pf, paths in sorted(results.items()):
            for fmt, path in paths.items():
                print(f"  {path}")
