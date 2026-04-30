"""Static Matplotlib mirror + single-spectrum plots (for PNG/SVG/PDF export).

Shares data-prep with the Plotly renderer via ``_mirror_shared``. Uses the
``Agg`` backend so this module is safe to import in worker processes.
"""
from io import BytesIO
from typing import Optional, Tuple

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


def _add_matplotlib_spectrum(ax, df: pd.DataFrame, y_multiplier: float,
                             show_annotations: bool = True,
                             annotation_fontsize: int = 7,
                             ann_offset: float = ANN_OFFSET):
    """Add matplotlib vlines and annotation text for one spectrum.

    Args:
        ax: matplotlib Axes to draw on.
        df: prepared DataFrame (must already have '_color' column).
        y_multiplier: +1 for upper/single, -1 for lower (mirror).
        show_annotations: whether to add text annotations on peaks.
        annotation_fontsize: font size for annotations.
        ann_offset: vertical offset for annotation text.
    """
    if df.empty:
        return

    for color, group in df.groupby("_color", sort=False):
        ax.vlines(
            group["mz"].values, 0, group["intensity"].values,
            colors=color, linewidth=1.5, alpha=0.8,
        )

    if show_annotations:
        sign = 1 if y_multiplier > 0 else -1
        va_label = "bottom" if y_multiplier > 0 else "top"
        ann_mask = _get_annotated_mask(df)
        for _, peak in df[ann_mask].iterrows():
            y_pos = peak["intensity"] + sign * ann_offset
            ax.text(peak["mz"], y_pos, peak["annotation"],
                    fontsize=annotation_fontsize, color=peak["_color"],
                    ha="center", va=va_label, rotation=0)


def create_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    figure_size: Tuple[int, int] = (12, 8),
    export_resolution: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7,
    metrics_text: Optional[str] = None,
    normalize: bool = True,
):
    """Create a publication-ready static matplotlib mirror plot."""
    fig = MplFigure(figsize=figure_size, dpi=export_resolution)
    ax = fig.add_subplot(111)

    if normalize:
        ann_offset = ANN_OFFSET
    else:
        upper_max = upper_spectrum["intensity"].max() if not upper_spectrum.empty else 1
        lower_max = lower_spectrum["intensity"].max() if not lower_spectrum.empty else 1
        ann_offset = max(upper_max, lower_max) * ANN_OFFSET

    if not upper_spectrum.empty:
        upper_frame = _prepare_spectrum_data(upper_spectrum, normalize=normalize)
        _add_matplotlib_spectrum(ax, upper_frame, y_multiplier=1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

    if not lower_spectrum.empty:
        lower_frame = _prepare_spectrum_data(lower_spectrum, normalize=normalize)
        lower_frame["intensity"] = lower_frame["intensity"] * -1
        _add_matplotlib_spectrum(ax, lower_frame, y_multiplier=-1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

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

    ax.axhline(y=0, color="black", linewidth=1, alpha=0.5)

    ax.text(0.02, 0.99, upper_label, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha="left", color="#333",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8))
    ax.text(0.02, 0.01, lower_label, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="left", color="#333",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8))

    if metrics_text:
        ax.text(0.98, 0.82, metrics_text, transform=ax.transAxes,
                fontsize=10, fontfamily="monospace", va="top", ha="right", color="#333",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#CCC", alpha=0.85))

    ax.set_xlabel("m/z", fontsize=12, fontweight="bold")
    ax.set_ylabel("Relative Intensity" if normalize else "Intensity",
                  fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    return fig


def create_single_spectrum_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    figure_size: Tuple[int, int] = (12, 8),
    export_resolution: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7,
    normalize: bool = True,
):
    """Create a static matplotlib single annotated spectrum plot."""
    fig = MplFigure(figsize=figure_size, dpi=export_resolution)
    ax = fig.add_subplot(111)

    if not spectrum.empty:
        df = _prepare_spectrum_data(spectrum, normalize=normalize)
        if normalize:
            ann_offset = ANN_OFFSET
        else:
            ann_offset = df["intensity"].max() * ANN_OFFSET if len(df) > 0 else 1
        _add_matplotlib_spectrum(ax, df, y_multiplier=1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

    if not spectrum.empty:
        all_mz = spectrum["mz"].tolist()
        ax.set_xlim(min(all_mz) - MZ_WINDOW_PADDING, max(all_mz) + MZ_WINDOW_PADDING)

    if normalize:
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    else:
        y_max = spectrum["intensity"].max() * 1.15 if not spectrum.empty else 1
        ax.set_ylim(0, y_max)

    if label:
        ax.text(0.02, 0.99, label, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top", ha="left", color="#333",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8))

    ax.set_xlabel("m/z", fontsize=12, fontweight="bold")
    ax.set_ylabel("Relative Intensity" if normalize else "Intensity",
                  fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    return fig


def figure_to_bytes(fig, format: str = "png") -> bytes:
    """Convert a matplotlib figure to bytes."""
    buffer = BytesIO()
    fig.savefig(buffer, format=format, bbox_inches="tight", facecolor="white")
    buffer.seek(0)
    return buffer.read()
