"""Summary-plot layer: aggregate the comparison table and render the dual-axis
and ridgeline summaries.

`prepare_summary_data` turns a comparison table into a per-peptidoform aggregate
frame (count, median, mean, std, 95 % confidence interval, dominant charge).
`generate_all_summary_plots` runs that aggregation once and writes the requested
figure types to <output_dir>/summary/ in each requested file format, dispatching
to `plot_classic_dual_axis` (peptide-count bars twinned with metric boxplots)
and `plot_ridgeline` (stacked KDE density per peptidoform with count sidebar).

generate_all_summary_plots is called by core.export.run_export (and indirectly
cli.main); prepare_summary_data is also called directly by
streamlit_app.views.explore.explore_page for the in-app preview.
"""
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
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


def _metric_label(metric: str) -> str:
    """Return the human-readable axis label for metric.
    Falls back to the raw metric name itself when it is not a known key in
    _METRIC_LABELS, so an unrecognised metric still renders sensibly.
    """
    return _METRIC_LABELS.get(metric, metric)


def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    """Compute a figure height (in inches) that scales with the number of rows.
    Allocates per_row inches per peptidoform row but never returns less than
    minimum, keeping small plots legible while letting tall ones grow.
    """
    return max(minimum, per_row * n)


# ============================================================================
# Dual-axis renderer: peptide-count bars twinned with metric boxplots
# ============================================================================

# Palette for the two overlaid series — warm orange for the count bars, slate
# grey for the metric boxplots — chosen so the twinned axes read distinctly.
_DUAL_BAR_COLOR = "#E8945A"  # Count-bar fill
_DUAL_BAR_EDGE = "#D07030"   # Count-bar outline
_DUAL_BOX_COLOR = "#2C3E50"  # Boxplot lines (median, whiskers, caps, outliers)
_DUAL_BOX_FILL = "#BDC3C7"   # Boxplot box fill


def plot_classic_dual_axis(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render the dual-axis summary figure: count bars twinned with metric boxplots.
    Draws one horizontal row per peptidoform — a count bar on the bottom x axis overlaid
    with a boxplot of the chosen metric on a twinned top axis, ordered most-prevalent first.
    raw_frame supplies the per-psm metric values the boxplots summarise; aggregated_frame
    the per-peptidoform counts. Returns the Matplotlib Figure for the caller to save and close.
    """
    # Order rows so the most-counted peptidoform sits at the top of the chart.
    aggregated_frame = aggregated_frame.sort_values("count", ascending=False).reset_index(drop=True)
    peptidoform_count = len(aggregated_frame)
    # Auto-size height to the row count when the caller gave no explicit size.
    if figure_size is None:
        figure_size = (14, _auto_height(peptidoform_count))

    fig, ax = plt.subplots(figsize=figure_size, dpi=export_resolution)
    y_pos = np.arange(peptidoform_count)
    labels = aggregated_frame["peptidoform"].tolist()

    # Count bars on the bottom axis (one bar per peptidoform)
    ax.barh(y_pos, aggregated_frame["count"].values, height=0.65,
            color=_DUAL_BAR_COLOR, alpha=0.55, edgecolor=_DUAL_BAR_EDGE, linewidth=0.5, zorder=1)
    ax.set_xlabel("Peptide count", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.spines["right"].set_visible(False)

    # Metric boxplots on the twinned top axis (seaborn draws IQR / whiskers / outliers)
    ax2 = ax.twiny()
    sns.boxplot(
        data=raw_frame, y="peptidoform", x=metric, order=labels, ax=ax2,
        color=_DUAL_BOX_FILL, width=0.35, linewidth=0.8, fliersize=3,
        boxprops=dict(facecolor=_DUAL_BOX_FILL, edgecolor=_DUAL_BOX_COLOR, alpha=0.85),
        medianprops=dict(color=_DUAL_BOX_COLOR, linewidth=1.5),
        whiskerprops=dict(color=_DUAL_BOX_COLOR, linewidth=0.8),
        capprops=dict(color=_DUAL_BOX_COLOR, linewidth=0.8),
        flierprops=dict(marker=".", markerfacecolor=_DUAL_BOX_COLOR,
                        markeredgecolor="none", alpha=0.4),
    )
    # Without an explicit title, label the metric axis here; an explicit title suppresses
    # this label and is added as a figure suptitle below instead.
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
        Patch(facecolor=_DUAL_BAR_COLOR, alpha=0.55, edgecolor=_DUAL_BAR_EDGE, label="Peptide Count"),
        Line2D([0], [0], color=_DUAL_BOX_COLOR, marker="s", markersize=6,
               markerfacecolor=_DUAL_BOX_FILL, linestyle="-", linewidth=0.8,
               label=_metric_label(metric)),
    ], loc="upper left", bbox_to_anchor=(0.0, -0.04), fontsize=9, framealpha=0.9, ncol=2)

    # Reserve top margin for a suptitle when one is given, else keep it tight.
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.22, top=0.92 if title else 0.96)
    return fig


# ============================================================================
# Ridgeline renderer: stacked KDE density per peptidoform with count sidebar
# ============================================================================

# Fill and edge colours for the count sidebar bars — a warm orange that stays
# distinct from the RdYlGn ridge colormap.
_RIDGE_BAR_COLOR = "#E8945A"
_RIDGE_BAR_EDGE = "#D07030"
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
    peptidoform_count = len(aggregated_frame)
    # Auto-size the figure height so every ridge gets vertical room; width is
    # fixed at 14 inches.
    if figure_size is None:
        figure_size = (14, _auto_height(peptidoform_count, per_row=0.5))

    # Colour each ridge by its median value (RdYlGn so high = green, low = red).
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0.0, vmax=1.0)

    # Pull the per-peptidoform columns into plain arrays — already ordered by
    # descending count by prepare_summary_data.
    peptidoforms = aggregated_frame["peptidoform"].tolist()
    medians = aggregated_frame["median"].to_numpy(dtype=float)
    counts = aggregated_frame["count"].to_numpy()

    # Two panels sharing the y-axis: ridges on the left, a count sidebar on the right.
    fig, (ridge_ax, bar_ax) = plt.subplots(
        1, 2, figsize=figure_size, sharey=True,
        gridspec_kw={"width_ratios": [5, 1], "wspace": 0.03},
    )
    fig.set_dpi(export_resolution)

    # === Ridges ===
    # Drawn bottom (highest count) to top (lowest count): aggregated row i sits on
    # baseline y = peptidoform_count - 1 - i, and a lower ridge overdraws the one above.
    for i in range(peptidoform_count):
        baseline = peptidoform_count - 1 - i
        # Gather this peptidoform's metric values from the raw rows (one row per
        # measurement, one column per metric), dropping NaN/inf so the KDE only
        # sees usable numbers.
        values = raw_frame.loc[raw_frame["peptidoform"] == peptidoforms[i], metric].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        # No finite values — skip drawing this ridge entirely.
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
        ridge_ax.fill_between(_RIDGE_GRID, baseline, ridge, color=cmap(norm(medians[i])),
                              alpha=0.85, zorder=i)
        ridge_ax.plot(_RIDGE_GRID, ridge, color="white", linewidth=0.8, zorder=i)

    ridge_ax.set_xlim(0.0, 1.0)
    ridge_ax.set_ylim(-0.4, peptidoform_count - 1 + _RIDGE_HEIGHT + 0.4)
    ridge_ax.set_yticks(range(peptidoform_count))
    ridge_ax.set_yticklabels(peptidoforms[::-1])  # y = 0 is the last (highest-count) peptidoform
    ridge_ax.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ridge_ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ridge_ax.spines[side].set_visible(False)

    # === Count sidebar ===
    # Shares the ridge y-axis, so bar i lines up with the baseline of ridge i.
    bar_ax.barh(range(peptidoform_count), counts[::-1], height=0.6,
                color=_RIDGE_BAR_COLOR, alpha=0.8, edgecolor=_RIDGE_BAR_EDGE, linewidth=0.5)
    bar_ax.set_xlabel("Count", fontsize=10, fontweight="bold")
    bar_ax.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "right", "left"):
        bar_ax.spines[side].set_visible(False)

    # === Median colorbar ===
    fig.subplots_adjust(bottom=0.16)
    scalar_mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    scalar_mappable.set_array([])
    cbar_ax = fig.add_axes([0.25, 0.04, 0.35, 0.018])
    fig.colorbar(scalar_mappable, cax=cbar_ax, orientation="horizontal").set_label(
        f"Median {_metric_label(metric)}", fontsize=9,
    )

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.99)
    return fig


# ============================================================================
# Summary-plot orchestrator: aggregate + dispatch to renderers
# ============================================================================

def prepare_summary_data(
    source: str | Path | pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the comparison table and aggregate per-peptidoform stats for metric.
    source is either a CSV path or an already-built comparison DataFrame (the latter is
    copied, so callers' frames are never mutated). Returns (raw_frame, aggregated_frame):
    the per-row table augmented with a peptidoform column, plus one row per peptidoform
    (count, median, mean, std, quartiles, 95% CI, dominant charge). Raises ValueError when metric is absent.
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
    modification = raw_frame["Modification"].fillna("").astype(str)
    inline = modification.mask(
        modification.isin(["", "(unmodified)"]), raw_frame["Peptide"].astype(str),
    )
    raw_frame["peptidoform"] = inline + "/" + raw_frame["Charge"].astype(str)

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

    # Finalise: sort ascending by count so plots draw least→most prevalent
    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_frame, agg


def generate_all_summary_plots(
    source: str | Path | pd.DataFrame,
    output_dir: str | Path,
    metric: str = DEFAULTS.summary_metric,
    plot_types: list[str] | None = None,
    formats: list[str] | None = None,
    export_resolution: int = 300,
) -> dict[str, dict[str, str]]:
    """Render every requested summary plot type into <output_dir>/summary/.
    Aggregates once via prepare_summary_data, then renders each recognised type
    (dual_axis, ridgeline) once and saves it under every requested format. Unknown types,
    renderer ValueErrors, and per-file OSErrors are logged and skipped, never raised.
    Returns {plot_type: {format: path}} covering only the files written successfully.
    """
    # Resolve defaults: fall back to configured plot types + PNG-only output
    if plot_types is None:
        plot_types = list(DEFAULTS.summary_plot_types)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Prep data: aggregate once, reuse across every requested plot type
    raw_frame, aggregated_frame = prepare_summary_data(source, metric=metric)
    logger.info(f"Summary data: {len(aggregated_frame)} peptidoforms, {len(raw_frame)} rows, metric={metric}")

    # Render & save: dispatch each plot type directly, write per-format files
    results: dict[str, dict[str, str]] = {}
    for plot_type in plot_types:
        if plot_type == "dual_axis":
            renderer = plot_classic_dual_axis
        elif plot_type == "ridgeline":
            renderer = plot_ridgeline
        else:
            logger.warning("Unknown plot type: %s (skipping)", plot_type)
            continue

        try:
            fig = renderer(raw_frame, aggregated_frame, metric=metric,
                           export_resolution=export_resolution, title=None)
        except ValueError as error:
            logger.error("Failed to create %s plot: %s", plot_type, error)
            continue

        # The ridgeline for the default metric gets a special basename; everything
        # else is named after its plot type.
        basename = ("ridgeline_r2_binned"
                    if plot_type == "ridgeline" and metric == DEFAULTS.summary_metric
                    else plot_type)
        paths: dict[str, str] = {}
        for output_format in formats:
            output_path = summary_dir / f"{basename}.{output_format}"
            try:
                fig.savefig(str(output_path), format=output_format, dpi=export_resolution,
                            bbox_inches="tight", facecolor="white")
                paths[output_format] = str(output_path)
            except OSError as error:
                logger.error("Failed to save %s.%s: %s", basename, output_format, error)
        plt.close(fig)
        if paths:
            results[plot_type] = paths
    return results
