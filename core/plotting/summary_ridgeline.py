"""Ridgeline summary: stacked KDE density per peptidoform with count sidebar (dispatched by ``plotting.summary_pipeline``)."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import pandas as pd  # noqa: E402
from joypy import joyplot  # noqa: E402

from core.defaults import DEFAULTS  # noqa: E402
from core.plotting._summary_shared import _auto_height, _metric_label  # noqa: E402

_BAR_COLOR = "#E8945A"
_BAR_EDGE = "#D07030"


def plot_ridgeline(
    raw_frame: pd.DataFrame,
    aggregated_frame: pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
    figure_size: tuple[float, float] | None = None,
    export_resolution: int = 300,
    title: str | None = None,
) -> Figure:
    """Render a ridgeline plot: one KDE ridge per peptidoform, coloured by median, with a count sidebar and median colorbar."""
    variant_count = len(aggregated_frame)
    if figure_size is None:
        figure_size = (14, _auto_height(variant_count, per_row=0.5))

    # Colour each ridge by its median value (RdYlGn so high = green, low = red).
    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=0.0, vmax=1.0)
    colors = [cmap(norm(m)) for m in aggregated_frame["median"]]

    # Stacked KDE: joypy draws one ridge per peptidoform in aggregated_frame order
    fig, _ = joyplot(
        data=raw_frame, by="peptidoform", column=metric,
        labels=aggregated_frame["peptidoform"].tolist(),
        color=colors, overlap=0.65, x_range=(0, 1),
        figsize=figure_size, linewidth=0.8, fill=True,
    )
    fig.set_dpi(export_resolution)
    fig.axes[-1].set_xlabel(_metric_label(metric), fontsize=11, fontweight="bold")

    # Count sidebar: shrink the ridges, add a horizontal-bar axis in the freed strip
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.78, pos.height])
    bar_ax = fig.add_axes([0.78, 0.14, 0.14, 0.78])
    bar_ax.barh(range(variant_count), aggregated_frame["count"][::-1].values,
                color=_BAR_COLOR, alpha=0.7, edgecolor=_BAR_EDGE, linewidth=0.5)
    bar_ax.set_xlabel("Count", fontsize=10, fontweight="bold")
    bar_ax.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "right", "left"):
        bar_ax.spines[side].set_visible(False)

    # Median colorbar at the bottom of the figure
    scalar_mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    scalar_mappable.set_array([])
    cbar_ax = fig.add_axes([0.25, 0.02, 0.35, 0.015])
    fig.colorbar(scalar_mappable, cax=cbar_ax, orientation="horizontal").set_label(
        f"Median {_metric_label(metric)}", fontsize=9,
    )

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    return fig
