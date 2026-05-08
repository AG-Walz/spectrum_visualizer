"""Static Matplotlib mirror + single-spectrum plots (for PNG/SVG/PDF export).

Shares data-prep with the Plotly renderer via ``_mirror_shared``. Uses the
``Agg`` backend so this module is safe to import in worker processes.
"""
from io import BytesIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker  # noqa: E402 (kept after matplotlib.use)
from matplotlib.figure import Figure as MplFigure  # noqa: E402
import pandas as pd  # noqa: E402

from mspecviz.plotting._mirror_shared import (
    ANN_OFFSET,
    MZ_WINDOW_PADDING,
    _get_annotated_mask,
    _prepare_spectrum_data,
)


def _add_matplotlib_spectrum(ax, df: pd.DataFrame, y_multiplier: float, ann_offset: float):
    """Add matplotlib vlines and per-peak annotation text for one spectrum."""
    if df.empty:
        return

    for color, group in df.groupby("_color", sort=False):
        ax.vlines(
            group["mz"].values, 0, group["intensity"].values,
            colors=color, linewidth=1.5, alpha=0.8,
        )

    sign = 1 if y_multiplier > 0 else -1
    va_label = "bottom" if y_multiplier > 0 else "top"
    for _, peak in df[_get_annotated_mask(df)].iterrows():
        y_pos = peak["intensity"] + sign * ann_offset
        ax.text(peak["mz"], y_pos, peak["annotation"],
                fontsize=7, color=peak["_color"],
                ha="center", va=va_label, rotation=0)


def _set_mirror_axes(ax, upper_spectrum: pd.DataFrame, lower_spectrum: pd.DataFrame,
                    normalize: bool) -> None:
    """Apply x/y limits and tick formatting for a mirror plot."""
    all_mz = []
    if not upper_spectrum.empty:
        all_mz.extend(upper_spectrum["mz"].tolist())
    if not lower_spectrum.empty:
        all_mz.extend(lower_spectrum["mz"].tolist())
    if all_mz:
        ax.set_xlim(min(all_mz) - MZ_WINDOW_PADDING, max(all_mz) + MZ_WINDOW_PADDING)

    if normalize:
        ax.set_ylim(-1.2, 1.2)
        ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
        ax.set_yticklabels(["1.0", "0.5", "0", "0.5", "1.0"])
    else:
        upper_max = upper_spectrum["intensity"].max() if not upper_spectrum.empty else 1
        lower_max = lower_spectrum["intensity"].max() if not lower_spectrum.empty else 1
        y_max = max(upper_max, lower_max) * 1.2
        ax.set_ylim(-y_max, y_max)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):.0f}")
        )


def _set_single_axes(ax, spectrum: pd.DataFrame, normalize: bool) -> None:
    """Apply x/y limits + tick layout for a single (non-mirrored) plot."""
    if not spectrum.empty:
        all_mz = spectrum["mz"].tolist()
        ax.set_xlim(min(all_mz) - MZ_WINDOW_PADDING, max(all_mz) + MZ_WINDOW_PADDING)
    if normalize:
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    else:
        y_max = spectrum["intensity"].max() * 1.15 if not spectrum.empty else 1
        ax.set_ylim(0, y_max)


def _draw_corner_label(ax, text: str, *, top: bool) -> None:
    """Draw a small bold corner label (top-left when ``top`` else bottom-left)."""
    y, va = (0.99, "top") if top else (0.01, "bottom")
    ax.text(0.02, y, text, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va=va, ha="left", color="#333",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="none", alpha=0.8))


def _draw_metrics_box(ax, metrics_text: str) -> None:
    """Draw the upper-right metrics monospace text box."""
    ax.text(0.98, 0.82, metrics_text, transform=ax.transAxes,
            fontsize=10, fontfamily="monospace", va="top", ha="right", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#CCC", alpha=0.85))


def _style_axes(ax, title: str, normalize: bool) -> None:
    """Apply the shared title/labels/spines styling for spectrum plots."""
    ax.set_xlabel("m/z", fontsize=12, fontweight="bold")
    ax.set_ylabel("Relative Intensity" if normalize else "Intensity",
                  fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _create_spectrum_figure(
    spectra: list[tuple[pd.DataFrame, str, float]],
    *, title: str, normalize: bool, figure_size: tuple[int, int],
    export_resolution: int, mirror: bool, metrics_text: str | None = None,
) -> MplFigure:
    """Build a Matplotlib figure; ``spectra`` is ``[(df, label, y_multiplier)]``."""
    fig = MplFigure(figsize=figure_size, dpi=export_resolution)
    ax = fig.add_subplot(111)

    # Annotation y-offset: constant when normalised, else fraction of max intensity.
    if normalize:
        ann_offset = ANN_OFFSET
    else:
        max_int = max(
            (df["intensity"].max() for df, _, _ in spectra if not df.empty),
            default=1,
        )
        ann_offset = max_int * ANN_OFFSET

    for df, label, y_mult in spectra:
        if not df.empty:
            prep = _prepare_spectrum_data(df, normalize=normalize)
            if y_mult < 0:
                prep["intensity"] = prep["intensity"] * -1
            _add_matplotlib_spectrum(ax, prep, y_multiplier=y_mult, ann_offset=ann_offset)
        if label:
            _draw_corner_label(ax, label, top=(y_mult > 0))

    if mirror:
        _set_mirror_axes(ax, spectra[0][0], spectra[1][0], normalize)
        ax.axhline(y=0, color="black", linewidth=1, alpha=0.5)
    else:
        _set_single_axes(ax, spectra[0][0], normalize)

    if metrics_text:
        _draw_metrics_box(ax, metrics_text)
    _style_axes(ax, title, normalize)
    fig.tight_layout()
    return fig


def create_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    figure_size: tuple[int, int] = (12, 8),
    export_resolution: int = 300,
    metrics_text: str | None = None,
    normalize: bool = True,
) -> MplFigure:
    """Create a publication-ready static matplotlib mirror plot."""
    return _create_spectrum_figure(
        [(upper_spectrum, upper_label, 1.0), (lower_spectrum, lower_label, -1.0)],
        title=title, normalize=normalize,
        figure_size=figure_size, export_resolution=export_resolution,
        mirror=True, metrics_text=metrics_text,
    )


def create_single_spectrum_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    figure_size: tuple[int, int] = (12, 8),
    export_resolution: int = 300,
    normalize: bool = True,
) -> MplFigure:
    """Create a static matplotlib single annotated spectrum plot."""
    return _create_spectrum_figure(
        [(spectrum, label, 1.0)],
        title=title, normalize=normalize,
        figure_size=figure_size, export_resolution=export_resolution,
        mirror=False,
    )


def figure_to_bytes(fig, format: str = "png") -> bytes:
    """Convert a matplotlib figure to bytes."""
    buffer = BytesIO()
    fig.savefig(buffer, format=format, bbox_inches="tight", facecolor="white")
    buffer.seek(0)
    return buffer.read()
