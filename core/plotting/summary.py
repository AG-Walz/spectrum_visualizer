"""Summary-plot layer: aggregate the comparison table and render the ridgeline
and multi-boxplot summaries.

`prepare_summary_data` turns a comparison table into a per-peptidoform aggregate
frame (count, median, mean, std, 95 % confidence interval, dominant charge).
`generate_all_summary_plots` runs that aggregation once and writes the requested
figure types to <output_dir>/summary/ in each requested file format, dispatching
to `plot_ridgeline` (stacked KDE density per peptidoform with count sidebar) and
`plot_multi_boxplot` (one horizontal boxplot per peptidoform with all data points).

generate_all_summary_plots is called by core.export.run_export (and indirectly
cli.main); prepare_summary_data is also called directly by
streamlit_app.views.explore.explore_page for the in-app preview.
"""
import logging
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from core.defaults import DEFAULTS  # noqa: E402

logger = logging.getLogger(__name__)


# ============================================================================
# Shared metric metadata + figure-height helper
# ============================================================================

# Ordered comparison-metric columns expected in the summary DataFrame — consumed
# by generate_all_summary_plots to select metric columns and re-exported via
# core.plotting.__init__. Order defines the layout order of summary plots.
METRIC_COLUMNS = [
    "R² match", "R² binned", "Cosine", "Entropy",
    "Dice", "Jaccard", "SharedRatio",
]

# Maps each raw METRIC_COLUMNS name to a verbose axis-label string used when
# rendering summary plots.
_METRIC_LABELS = {
    "R² match": "R² (match)",
    "R² binned": "R² (binned)",
    "Cosine": "Cosine similarity",
    "Entropy": "Spectral entropy similarity",
    "Dice": "Dice coefficient",
    "Jaccard": "Jaccard index",
    "SharedRatio": "Shared peak ratio",
}


# Peptidoform ordering modes for the summary plots: mode key -> (aggregate column,
# ascending). Shared single source of truth for prepare_summary_data and the two front
# ends. "ascending" reads top→bottom in the rendered chart (row 0 sits at the top).
SUMMARY_SORT_MODES = {
    "alpha_asc": ("peptidoform", True),
    "alpha_desc": ("peptidoform", False),
    "count_asc": ("count", True),
    "count_desc": ("count", False),
    "median_asc": ("median", True),
    "median_desc": ("median", False),
}


def _metric_label(metric: str) -> str:
    """Return the human-readable axis label for metric.
    Falls back to the raw metric name itself when it is not a known key in
    _METRIC_LABELS, so an unrecognised metric still renders sensibly.
    """
    return _METRIC_LABELS.get(metric, metric)


def _metric_slug(metric: str) -> str:
    """Filesystem-safe lowercase slug for a metric name ('R² binned' -> 'r2_binned')."""
    return re.sub(r"[^a-z0-9]+", "_", metric.lower().replace("²", "2")).strip("_")


def peptidoform_label_series(frame: pd.DataFrame) -> pd.Series:
    """Return the per-row peptidoform label ("<modified-or-base-sequence>/<charge>").
    Single source of truth for the label prepare_summary_data groups on and the Explore
    peptide pickers list; uses the Modification column when present, else the base Peptide.
    """
    modification = frame["Modification"].fillna("").astype(str)
    inline = modification.mask(
        modification.isin(["", "(unmodified)"]), frame["Peptide"].astype(str),
    )
    return inline + "/" + frame["Charge"].astype(str)


def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    """Compute a figure height (in inches) that scales with the number of rows.
    Allocates per_row inches per peptidoform row but never returns less than
    minimum, keeping small plots legible while letting tall ones grow.
    """
    return max(minimum, per_row * n)


# ============================================================================
# Shared median-coloured two-panel scaffolding (ridgeline + multi-boxplot)
# ============================================================================
# Both renderers draw the same frame: a wide main panel (one row per peptidoform,
# coloured by median) beside a narrow peptide-count sidebar, under a median
# colorbar. _summary_panels builds the frame and pulls the shared per-peptidoform
# arrays; _finalize_summary applies the shared axis dressing, sidebar, colorbar,
# and title once a renderer has drawn its own row content. Only the middle drawing
# loop (KDE ridge vs boxplot + points) differs between the two.

# Count-sidebar bar colours — a warm orange that stays distinct from the RdYlGn
# median colormap.
_COUNT_BAR_COLOR = "#E8945A"
_COUNT_BAR_EDGE = "#D07030"
# Median colour scale: high metric = green, low = red, fixed to the [0, 1] range
# every comparison metric is normalised to.
_MEDIAN_CMAP = plt.cm.RdYlGn
_MEDIAN_NORM = plt.Normalize(vmin=0.0, vmax=1.0)
# Faint horizontal rule drawn once per peptidoform row (behind the box/ridge) so the
# eye can trace a peptide's y-label across to its own row.
_ROW_GUIDE_COLOR = "#CCCCCC"


def _finite_values(raw_frame: pd.DataFrame, peptidoform: str, metric: str) -> np.ndarray:
    """Return one peptidoform's metric values as a float array with NaN/inf dropped."""
    values = raw_frame.loc[raw_frame["peptidoform"] == peptidoform, metric].to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _summary_panels(aggregated_frame, figure_size, export_resolution, per_row):
    """Build the shared two-panel figure and pull the per-peptidoform arrays.
    Returns (fig, main_ax, bar_ax, peptidoforms, medians, counts): a wide main axis (5)
    beside a narrow count sidebar (1) sharing the y-axis. figure_size defaults to 14 inches
    wide by per_row inches per peptidoform. Arrays stay in prepare_summary_data's ascending-count order.
    """
    if figure_size is None:
        figure_size = (14, _auto_height(len(aggregated_frame), per_row=per_row))
    fig, (main_ax, bar_ax) = plt.subplots(
        1, 2, figsize=figure_size, sharey=True,
        gridspec_kw={"width_ratios": [5, 1], "wspace": 0.03},
    )
    fig.set_dpi(export_resolution)
    return (
        fig, main_ax, bar_ax,
        aggregated_frame["peptidoform"].tolist(),
        aggregated_frame["median"].to_numpy(dtype=float),
        aggregated_frame["count"].to_numpy(),
    )


def _finalize_summary(fig, main_ax, bar_ax, peptidoforms, counts, metric, xlim, ylim, title):
    """Dress the figure once a renderer has drawn its rows: main-axis limits/labels, the
    peptide-count sidebar, the median colorbar, and an optional figure title.
    """
    # === Main axis ===
    main_ax.set_xlim(*xlim)
    main_ax.set_ylim(*ylim)
    main_ax.set_yticks(range(len(peptidoforms)))
    main_ax.set_yticklabels(peptidoforms[::-1])  # y = 0 is the last (highest-count) peptidoform
    main_ax.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    main_ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        main_ax.spines[side].set_visible(False)

    # === Count sidebar (bar i aligns with main-axis row i, top row = highest count) ===
    bar_ax.barh(range(len(counts)), counts[::-1], height=0.6,
                color=_COUNT_BAR_COLOR, alpha=0.8, edgecolor=_COUNT_BAR_EDGE, linewidth=0.5)
    bar_ax.set_xlabel("Peptide Count", fontsize=10, fontweight="bold")
    bar_ax.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "right", "left"):
        bar_ax.spines[side].set_visible(False)

    # === Median colorbar ===
    fig.subplots_adjust(bottom=0.16)
    scalar_mappable = plt.cm.ScalarMappable(cmap=_MEDIAN_CMAP, norm=_MEDIAN_NORM)
    scalar_mappable.set_array([])
    cbar_ax = fig.add_axes([0.25, 0.04, 0.35, 0.018])
    fig.colorbar(scalar_mappable, cax=cbar_ax, orientation="horizontal").set_label(
        f"Median {_metric_label(metric)}", fontsize=9,
    )

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)


# ============================================================================
# Ridgeline renderer: stacked KDE density per peptidoform with count sidebar
# ============================================================================

# Evaluation grid for every KDE: 256 evenly-spaced points across the [0, 1]
# metric range (all comparison metrics are normalised to that range).
_RIDGE_GRID = np.linspace(0.0, 1.0, 256)
_RIDGE_HEIGHT = 2.0      # peak height of a ridge in baseline units; > 1 makes neighbours overlap
_RIDGE_DEGENERATE_SIGMA = 0.02  # bump width for peptidoforms with one (or near-constant) value, where a KDE is undefined


def plot_ridgeline(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render a ridgeline plot: one KDE ridge per peptidoform, coloured by its median, with
    a count sidebar and a median colorbar. raw_frame supplies the per-psm values each KDE
    ridge is estimated from; aggregated_frame holds one row per peptidoform — ordered by
    descending count — and controls ridge order, colour, and the sidebar. Returns the
    Matplotlib Figure (DPI export_resolution) for the caller to save or close.
    """
    # Shared two-panel scaffold + per-peptidoform arrays (ridges need more row room).
    fig, ridge_ax, bar_ax, peptidoforms, medians, counts = _summary_panels(
        aggregated_frame, figure_size, export_resolution, per_row=0.5)
    peptidoform_count = len(peptidoforms)

    # === Row guide lines ===
    # A faint rule at each ridge baseline so every ridge ties back to its y-label.
    for baseline in range(peptidoform_count):
        ridge_ax.axhline(baseline, color=_ROW_GUIDE_COLOR, linewidth=0.6, zorder=0)

    # === Ridges ===
    # Drawn bottom (highest count) to top (lowest count): aggregated row i sits on
    # baseline y = peptidoform_count - 1 - i, and a lower ridge overdraws the one above.
    for i in range(peptidoform_count):
        baseline = peptidoform_count - 1 - i
        # Gather this peptidoform's finite metric values; skip the ridge entirely
        # when none survive (the KDE has nothing to estimate from).
        values = _finite_values(raw_frame, peptidoforms[i], metric)
        if values.size == 0:
            continue

        # KDE over the [0, 1] grid; fall back to a narrow bump when a single value
        # (or zero spread) makes the KDE covariance singular.
        if values.size < 2 or np.ptp(values) < 1e-9:
            density = np.exp(-0.5 * ((_RIDGE_GRID - values.mean()) / _RIDGE_DEGENERATE_SIGMA) ** 2)
        else:
            try:
                density = gaussian_kde(values)(_RIDGE_GRID)
            except np.linalg.LinAlgError:
                density = np.exp(-0.5 * ((_RIDGE_GRID - values.mean()) / _RIDGE_DEGENERATE_SIGMA) ** 2)
        # A flat (all-zero) density carries no shape — skip rather than divide
        # by zero when normalising below.
        peak = float(density.max())
        if peak <= 0.0:
            continue

        # Normalise each density to a unit peak then scale to _RIDGE_HEIGHT, so
        # every ridge has the same maximum height regardless of its spread.
        ridge = baseline + density / peak * _RIDGE_HEIGHT
        ridge_ax.fill_between(_RIDGE_GRID, baseline, ridge, color=_MEDIAN_CMAP(_MEDIAN_NORM(medians[i])),
                              alpha=0.85, zorder=i)
        ridge_ax.plot(_RIDGE_GRID, ridge, color="white", linewidth=0.8, zorder=i)

    # Ridges overlap upward by _RIDGE_HEIGHT, so the top row needs extra headroom.
    _finalize_summary(fig, ridge_ax, bar_ax, peptidoforms, counts, metric,
                      xlim=(0.0, 1.0), ylim=(-0.4, peptidoform_count - 1 + _RIDGE_HEIGHT + 0.4),
                      title=title)
    return fig


# ============================================================================
# Multi-boxplot renderer: one horizontal boxplot per peptidoform with all points
# ============================================================================

# Edge/median/point colours shared by every box — neutral greys that keep the
# RdYlGn median fill (and the orange count sidebar) the only carriers of meaning.
_BOX_EDGE_COLOR = "#444444"     # Box outline, whiskers, caps
_BOX_MEDIAN_COLOR = "#222222"   # Median line inside each box
_BOX_POINT_COLOR = "#333333"    # Overlaid per-psm data dots
_BOX_POINT_JITTER = 0.18        # Vertical jitter (in row units) spread for the dots


def plot_multi_boxplot(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render one horizontal boxplot per peptidoform, the box filled by its median, with every
    per-psm value overlaid as a jittered dot, plus a count sidebar and median colorbar. Shares
    the ridgeline layout: raw_frame supplies the per-psm values, aggregated_frame holds one row
    per peptidoform ordered by descending count. Returns the Matplotlib Figure for the caller.
    """
    # Shared two-panel scaffold + per-peptidoform arrays (boxes pack tighter than ridges).
    fig, box_ax, bar_ax, peptidoforms, medians, counts = _summary_panels(
        aggregated_frame, figure_size, export_resolution, per_row=0.4)
    peptidoform_count = len(peptidoforms)
    # Seeded jitter so repeated renders of the same data are pixel-stable.
    rng = np.random.default_rng(0)

    # === Row guide lines ===
    # A faint rule through each row's centre so every box ties back to its y-label.
    for position in range(peptidoform_count):
        box_ax.axhline(position, color=_ROW_GUIDE_COLOR, linewidth=0.6, zorder=0)

    # === Boxplots ===
    # Position 0 sits at the bottom; draw the highest-count peptidoform there so the
    # ordering matches the ridgeline (most prevalent at the bottom of the chart).
    for position in range(peptidoform_count):
        index = peptidoform_count - 1 - position
        values = _finite_values(raw_frame, peptidoforms[index], metric)
        if values.size == 0:
            continue

        box_ax.boxplot(
            values, positions=[position], orientation="horizontal", widths=0.6,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor=_MEDIAN_CMAP(_MEDIAN_NORM(medians[index])),
                          edgecolor=_BOX_EDGE_COLOR, alpha=0.85, linewidth=0.8),
            medianprops=dict(color=_BOX_MEDIAN_COLOR, linewidth=1.4),
            whiskerprops=dict(color=_BOX_EDGE_COLOR, linewidth=0.8),
            capprops=dict(color=_BOX_EDGE_COLOR, linewidth=0.8),
        )
        # Overlay every value as a jittered dot so no measurement is hidden.
        jitter = rng.uniform(-_BOX_POINT_JITTER, _BOX_POINT_JITTER, size=values.size)
        box_ax.scatter(values, position + jitter, s=8, color=_BOX_POINT_COLOR,
                       alpha=0.5, linewidths=0, zorder=3)

    # Slight x padding so boxes touching 0 or 1 aren't clipped at the spine.
    _finalize_summary(fig, box_ax, bar_ax, peptidoforms, counts, metric,
                      xlim=(-0.02, 1.02), ylim=(-0.6, peptidoform_count - 0.4),
                      title=title)
    return fig


# ============================================================================
# Single-peptide reshape: measures on the y-axis for ONE peptidoform
# ============================================================================
# The renderers above lay one peptidoform per y-row for a single measure. For the
# single-peptide view (#105) we keep those exact renderers but swap what a "row" means:
# reshape one peptidoform's rows so each MEASURE becomes a row and its per-spectrum values
# fill a shared value column. Feed the result to plot_ridgeline / plot_multi_boxplot with
# metric=MEASURE_VALUE_COLUMN — no new plotting code.

# The value column the reshaped frame carries; also the x-axis label the renderers show.
MEASURE_VALUE_COLUMN = "Metric value"


def prepare_single_peptide_measures(
    frame: pd.DataFrame, metrics: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reshape one peptidoform's comparison rows into the (raw, aggregated) pair the population
    renderers expect, but with each measure in metrics as a "peptidoform" row (values under
    MEASURE_VALUE_COLUMN). Measures keep the given order. Raises ValueError when nothing is finite.
    """
    records = []
    for metric in metrics:
        if metric not in frame.columns:
            continue
        values = frame[metric].to_numpy(dtype=float)
        records.extend({"peptidoform": _metric_label(metric), MEASURE_VALUE_COLUMN: value}
                       for value in values[np.isfinite(values)])
    if not records:
        raise ValueError("no finite values for the selected measures")

    raw_frame = pd.DataFrame(records)
    aggregated_frame = (
        raw_frame.groupby("peptidoform", sort=False)[MEASURE_VALUE_COLUMN]
        .agg(["count", "median"]).reset_index()
    )
    return raw_frame, aggregated_frame


# ============================================================================
# Summary-plot orchestrator: aggregate + dispatch to renderers
# ============================================================================

def prepare_summary_data(
    source: str | Path | pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    sort_by: str = DEFAULTS.summary_sort,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the comparison table and aggregate per-peptidoform stats for metric.
    source is either a CSV path or an already-built comparison DataFrame (the latter is
    copied, so callers' frames are never mutated). sort_by picks the peptidoform ordering
    from SUMMARY_SORT_MODES (default: ascending count, the historical layout). Returns
    (raw_frame, aggregated_frame): the per-row table augmented with a peptidoform column,
    plus one row per peptidoform (count, median, mean, std, quartiles, 95% CI, dominant
    charge). Raises ValueError when metric is absent.
    """
    # Load: accept either an on-disk CSV path or an already-built DataFrame
    if isinstance(source, (str, Path)):
        raw_frame = pd.read_csv(source, low_memory=False)
    else:
        raw_frame = source.copy()

    # Validate: fail fast if the requested metric column is missing
    if metric not in raw_frame.columns:
        raise ValueError(
            f"Metric '{metric}' not found. Available: "
            f"{[c for c in raw_frame.columns if c in METRIC_COLUMNS]}"
        )

    # Build peptidoform label: modification-or-base sequence joined with charge
    raw_frame["peptidoform"] = peptidoform_label_series(raw_frame)

    # Drop rows with a missing value in metric before aggregating
    raw_frame = raw_frame.dropna(subset=[metric]).copy()

    # Aggregate: per-peptidoform count/quartile stats + 95% CI + dominant charge
    grouped = raw_frame.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg([
        "count", "median", "mean", "std",
        ("q25", lambda x: x.quantile(0.25)),
        ("q75", lambda x: x.quantile(0.75)),
    ])
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["charge"] = grouped["Charge"].agg(lambda x: x.mode().iloc[0])

    # Finalise: order peptidoforms per the requested mode. Renderers place aggregated
    # row 0 at the top of the chart, so "ascending" reads top→bottom (the default,
    # count_asc, keeps the historical least→most-prevalent layout).
    sort_column, ascending = SUMMARY_SORT_MODES.get(sort_by, SUMMARY_SORT_MODES[DEFAULTS.summary_sort])
    agg = agg.reset_index().sort_values(sort_column, ascending=ascending).reset_index(drop=True)
    return raw_frame, agg


# Maps each summary plot-type key to its renderer function. Public so the
# Streamlit explore view dispatches through the same single source of truth.
SUMMARY_RENDERERS = {
    "ridgeline": plot_ridgeline,
    "multi_boxplot": plot_multi_boxplot,
}


def generate_all_summary_plots(
    source: str | Path | pd.DataFrame,
    output_dir: str | Path,
    metrics: list[str] | None = None,
    plot_types: list[str] | None = None,
    formats: list[str] | None = None,
    export_resolution: int = 300,
    sort_by: str = DEFAULTS.summary_sort,
) -> dict[str, dict[str, str]]:
    """Render every requested plot type for every requested metric into <output_dir>/summary/.
    For each metric, aggregates once via prepare_summary_data (peptidoforms ordered per sort_by),
    then renders each recognised type (ridgeline, multi_boxplot) and saves it as
    <plot_type>_<metric_slug>.<format>. Unknown types, metrics absent from the table, renderer
    ValueErrors, and per-file OSErrors are logged and skipped, never raised.
    Returns {plot_type: {"<metric_slug>.<format>": path}}.
    """
    # Resolve defaults: configured metrics + plot types, PNG-only output
    if metrics is None:
        metrics = [DEFAULTS.summary_metric]
    if plot_types is None:
        plot_types = list(DEFAULTS.summary_plot_types)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Render & save: one figure per (metric, plot type), written under every format
    results: dict[str, dict[str, str]] = {}
    for metric in metrics:
        # Aggregate once per metric; a metric missing from the table is skipped, not fatal.
        try:
            raw_frame, aggregated_frame = prepare_summary_data(source, metric=metric, sort_by=sort_by)
        except ValueError as error:
            logger.error("Skipping metric %s: %s", metric, error)
            continue
        logger.info(f"Summary data: {len(aggregated_frame)} peptidoforms, {len(raw_frame)} rows, metric={metric}")

        for plot_type in plot_types:
            renderer = SUMMARY_RENDERERS.get(plot_type)
            if renderer is None:
                logger.warning("Unknown plot type: %s (skipping)", plot_type)
                continue
            try:
                fig = renderer(raw_frame, aggregated_frame, metric=metric,
                               export_resolution=export_resolution, title=None)
            except ValueError as error:
                logger.error("Failed to create %s plot for %s: %s", plot_type, metric, error)
                continue

            basename = f"{plot_type}_{_metric_slug(metric)}"
            for output_format in formats:
                output_path = summary_dir / f"{basename}.{output_format}"
                try:
                    fig.savefig(str(output_path), format=output_format, dpi=export_resolution,
                                bbox_inches="tight", facecolor="white")
                    results.setdefault(plot_type, {})[f"{_metric_slug(metric)}.{output_format}"] = str(output_path)
                except OSError as error:
                    logger.error("Failed to save %s.%s: %s", basename, output_format, error)
            plt.close(fig)
    return results
