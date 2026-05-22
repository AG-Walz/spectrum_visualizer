"""Ridgeline summary: stacked KDE density per peptidoform with count sidebar (dispatched by plotting.summary_pipeline).
Role — a leaf-level renderer in the summary-plot layer; one of two figure
renderers (alongside summary_dual_axis) that summary_pipeline can
dispatch to. What it produces — a single Matplotlib Figure showing, for one
chosen comparison metric, a stacked ridgeline: each peptidoform contributes a
KDE-shaped ridge coloured by its median, with a count sidebar and a shared
median colorbar. Call relationships — plot_ridgeline is imported by
core.plotting.__init__ and dispatched by summary_pipeline.create_summary_plot
on the pre-aggregated frames from prepare_summary_data (called by
core.export.export_pipeline and cli.main); it consumes _auto_height
and _metric_label from core.plotting._summary_shared and DEFAULTS
from core.defaults.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

from core.defaults import DEFAULTS  # noqa: E402
from core.plotting._summary_shared import _auto_height, _metric_label  # noqa: E402

# Fill and edge colours for the count sidebar bars — a warm orange that stays
# distinct from the RdYlGn ridge colormap.
_BAR_COLOR = "#E8945A"
_BAR_EDGE = "#D07030"
# Evaluation grid for every KDE: 256 evenly-spaced points across the [0, 1]
# metric range (all comparison metrics are normalised to that range).
_GRID = np.linspace(0.0, 1.0, 256)
_RIDGE_HEIGHT = 2.0      # peak height of a ridge in baseline units; > 1 makes neighbours overlap
_DEGENERATE_SIGMA = 0.02  # bump width for peptidoforms with one (or near-constant) value, where a KDE is undefined


def plot_ridgeline(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render a ridgeline plot: one KDE ridge per peptidoform, coloured by median, with a count sidebar and median colorbar.
    raw_frame holds the per-psm rows (one row per measurement, with a
    peptidoform column and one column per metric) and supplies the values
    that each KDE ridge is estimated from; aggregated_frame holds one row
    per peptidoform — ordered by descending count — with peptidoform,
    median, and count columns and controls ridge order, colour, and the
    sidebar. metric selects which column of raw_frame to plot;
    figure_size defaults to 14 inches wide with a height that scales with
    the peptidoform count. Returns the assembled Matplotlib Figure (DPI set
    to export_resolution); the caller is responsible for saving or closing it.
    """
    variant_count = len(aggregated_frame)
    # Auto-size the figure height so every ridge gets vertical room; width is
    # fixed at 14 inches.
    if figure_size is None:
        figure_size = (14, _auto_height(variant_count, per_row=0.5))

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
    # baseline y = variant_count - 1 - i, and a lower ridge overdraws the one above.
    for i in range(variant_count):
        baseline = variant_count - 1 - i
        # Gather this peptidoform's metric values from the raw rows, dropping
        # NaN/inf so the KDE only sees usable numbers.
        values = raw_frame.loc[raw_frame["peptidoform"] == peptidoforms[i], metric].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        # No finite values — skip drawing this ridge entirely.
        if values.size == 0:
            continue

        # KDE over the [0, 1] grid; fall back to a narrow bump when a single value
        # (or zero spread) makes the KDE covariance singular.
        if values.size < 2 or np.ptp(values) < 1e-9:
            density = np.exp(-0.5 * ((_GRID - values.mean()) / _DEGENERATE_SIGMA) ** 2)
        else:
            try:
                density = gaussian_kde(values)(_GRID)
            except np.linalg.LinAlgError:
                density = np.exp(-0.5 * ((_GRID - values.mean()) / _DEGENERATE_SIGMA) ** 2)
        # A flat (all-zero) density carries no shape — skip rather than divide
        # by zero when normalising below.
        peak = float(density.max())
        if peak <= 0.0:
            continue

        # Normalise each density to a unit peak then scale to _RIDGE_HEIGHT, so
        # every ridge has the same maximum height regardless of its spread.
        ridge = baseline + density / peak * _RIDGE_HEIGHT
        ridge_ax.fill_between(_GRID, baseline, ridge, color=cmap(norm(medians[i])),
                              alpha=0.85, zorder=i)
        ridge_ax.plot(_GRID, ridge, color="white", linewidth=0.8, zorder=i)

    ridge_ax.set_xlim(0.0, 1.0)
    ridge_ax.set_ylim(-0.4, variant_count - 1 + _RIDGE_HEIGHT + 0.4)
    ridge_ax.set_yticks(range(variant_count))
    ridge_ax.set_yticklabels(peptidoforms[::-1])  # y = 0 is the last (highest-count) peptidoform
    ridge_ax.set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")
    ridge_ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ridge_ax.spines[side].set_visible(False)

    # === Count sidebar ===
    # Shares the ridge y-axis, so bar i lines up with the baseline of ridge i.
    bar_ax.barh(range(variant_count), counts[::-1], height=0.6,
                color=_BAR_COLOR, alpha=0.8, edgecolor=_BAR_EDGE, linewidth=0.5)
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