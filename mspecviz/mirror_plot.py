"""
Mirror plot generation for MS/MS spectrum comparison.

Provides both static Matplotlib mirror plots (for PNG/SVG/PDF export)
and interactive Plotly mirror plots (for Streamlit display and HTML export).
"""
import logging
from io import BytesIO
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure as MplFigure
import pandas as pd
import plotly.graph_objects as go

from mspecviz.common.ion_utils import ION_COLORS, get_ion_color  # noqa: F401

logger = logging.getLogger(__name__)


def _normalize_intensities(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize df['intensity'] to [0, 1] in-place and return df."""
    max_int = df["intensity"].max()
    if max_int > 0:
        df["intensity"] = df["intensity"] / max_int
    return df


def _ensure_annotation_col(df: pd.DataFrame) -> pd.DataFrame:
    """Return df with an 'annotation' column, adding 'unannotated' if absent."""
    if "annotation" not in df.columns:
        df = df.copy()
        df["annotation"] = "unannotated"
    return df


def _prepare_spectrum_data(df: pd.DataFrame, normalize: bool = True,
                           color_map: Dict[str, str] = None) -> pd.DataFrame:
    """Prepare spectrum DataFrame for plotting (shared by Plotly and Matplotlib)."""
    df = _ensure_annotation_col(df.copy())
    if normalize:
        df = _normalize_intensities(df)
    if color_map is None:
        color_map = ION_COLORS
    df["_color"] = df["annotation"].map(lambda a: get_ion_color(a, color_map))
    return df


def _get_annotated_mask(group: pd.DataFrame) -> pd.Series:
    """Return boolean mask for peaks that have a real annotation (not '' or 'unannotated')."""
    return group["annotation"].ne("unannotated") & group["annotation"].ne("")


def _add_plotly_spectrum_traces(fig, df: pd.DataFrame, y_multiplier: float,
                                label: str, show_annotations: bool = True):
    """Add Plotly scatter traces and annotation labels for one spectrum.

    Args:
        fig: plotly Figure to add traces to.
        df: prepared DataFrame (must already have '_color' column).
        y_multiplier: +1 for upper/single spectrum, -1 for lower (mirror).
        label: hover label prefix (e.g. "Query", "Reference", or "").
        show_annotations: whether to add text annotations on peaks.
    """
    hover_prefix = f"{label}<br>" if label else ""

    for color, group in df.groupby("_color", sort=False):
        ann_mask = _get_annotated_mask(group)

        for _, peak in group.iterrows():
            fig.add_trace(go.Scatter(
                x=[peak["mz"], peak["mz"]],
                y=[0, y_multiplier * peak["intensity"]],
                mode="lines",
                line=dict(color=color, width=1.5),
                hovertext=f"{hover_prefix}m/z: {peak['mz']:.4f}<br>"
                          f"Int: {peak['intensity']:.4f}<br>{peak['annotation']}",
                hoverinfo="text",
                showlegend=False,
            ))

        if show_annotations:
            yshift = 10 if y_multiplier > 0 else -10
            for _, peak in group[ann_mask].iterrows():
                fig.add_annotation(
                    x=peak["mz"],
                    y=y_multiplier * peak["intensity"],
                    text=peak["annotation"],
                    showarrow=False,
                    font=dict(size=8, color=color),
                    yshift=yshift,
                )


def _add_matplotlib_spectrum(ax, df: pd.DataFrame, y_multiplier: float,
                             show_annotations: bool = True,
                             annotation_fontsize: int = 7,
                             ann_offset: float = 0.03):
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


# ---------------------------------------------------------------------------
# Plotly interactive plots
# ---------------------------------------------------------------------------

def create_plotly_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    color_map: Dict[str, str] = None,
    metrics_text: Optional[str] = None,
    normalize: bool = True,
    height: int = 600,
):
    """
    Create an interactive Plotly mirror plot.

    Args:
        upper_spectrum: DataFrame with columns [mz, intensity, annotation]
        lower_spectrum: DataFrame with columns [mz, intensity, annotation]
        title: Plot title
        upper_label: Label for upper spectrum
        lower_label: Label for lower spectrum
        color_map: Custom ion color mapping
        metrics_text: Metrics string for annotation overlay
        normalize: Normalize intensities to [0, 1]
        height: Plot height in pixels

    Returns:
        plotly.graph_objs.Figure
    """
    if color_map is None:
        color_map = ION_COLORS

    fig = go.Figure()

    # --- Upper spectrum (query, positive) ---
    if not upper_spectrum.empty:
        upper_df = _prepare_spectrum_data(upper_spectrum, normalize=normalize, color_map=color_map)
        _add_plotly_spectrum_traces(fig, upper_df, y_multiplier=1.0,
                                   label=upper_label, show_annotations=True)

    # --- Lower spectrum (reference, negative) ---
    if not lower_spectrum.empty:
        lower_df = _prepare_spectrum_data(lower_spectrum, normalize=normalize, color_map=color_map)
        _add_plotly_spectrum_traces(fig, lower_df, y_multiplier=-1.0,
                                   label=lower_label, show_annotations=True)

    # Center line
    all_mz = []
    if not upper_spectrum.empty:
        all_mz.extend(upper_spectrum["mz"].tolist())
    if not lower_spectrum.empty:
        all_mz.extend(lower_spectrum["mz"].tolist())

    if all_mz:
        mz_min = min(all_mz) - 50
        mz_max = max(all_mz) + 50
        fig.add_shape(type="line", x0=mz_min, x1=mz_max, y0=0, y1=0,
                      line=dict(color="black", width=1, dash="solid"))

    # Layout
    y_range = [-1.25, 1.25] if normalize else None
    fig.update_layout(
        title=title,
        xaxis_title="m/z",
        yaxis_title="Relative Intensity" if normalize else "Intensity",
        height=height,
        hovermode="closest",
        template="plotly_white",
        showlegend=False,
    )
    if y_range:
        fig.update_yaxes(range=y_range,
                         tickvals=[-1.0, -0.5, 0, 0.5, 1.0],
                         ticktext=["1.0", "0.5", "0", "0.5", "1.0"])

    # Labels
    fig.add_annotation(x=0.02, y=0.98, xref="paper", yref="paper",
                       text=f"<b>{upper_label}</b>", showarrow=False,
                       font=dict(size=11, color="#333"),
                       bgcolor="rgba(255,255,255,0.8)", xanchor="left", yanchor="top")
    fig.add_annotation(x=0.02, y=0.02, xref="paper", yref="paper",
                       text=f"<b>{lower_label}</b>", showarrow=False,
                       font=dict(size=11, color="#333"),
                       bgcolor="rgba(255,255,255,0.8)", xanchor="left", yanchor="bottom")

    # Metrics overlay
    if metrics_text:
        fig.add_annotation(x=0.98, y=0.98, xref="paper", yref="paper",
                           text=metrics_text.replace("\n", "<br>"),
                           showarrow=False,
                           font=dict(size=10, family="monospace", color="#333"),
                           bgcolor="rgba(255,255,255,0.85)",
                           bordercolor="#CCC", borderwidth=1,
                           xanchor="right", yanchor="top",
                           align="left")

    return fig


def create_plotly_single_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    color_map: Dict[str, str] = None,
    normalize: bool = True,
    height: int = 500,
):
    """Create an interactive Plotly single-spectrum plot."""
    if color_map is None:
        color_map = ION_COLORS

    fig = go.Figure()

    if not spectrum.empty:
        df = _prepare_spectrum_data(spectrum, normalize=normalize, color_map=color_map)
        _add_plotly_spectrum_traces(fig, df, y_multiplier=1.0,
                                   label=label, show_annotations=True)

    fig.update_layout(
        title=title,
        xaxis_title="m/z",
        yaxis_title="Relative Intensity" if normalize else "Intensity",
        height=height,
        hovermode="closest",
        template="plotly_white",
        showlegend=False,
    )
    if normalize:
        fig.update_yaxes(range=[0, 1.15])

    if label:
        fig.add_annotation(x=0.02, y=0.98, xref="paper", yref="paper",
                           text=f"<b>{label}</b>", showarrow=False,
                           font=dict(size=11, color="#333"),
                           bgcolor="rgba(255,255,255,0.8)",
                           xanchor="left", yanchor="top")

    return fig


# ---------------------------------------------------------------------------
# Matplotlib static plots (for PNG/SVG/PDF export)
# ---------------------------------------------------------------------------

def create_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    color_map: Dict[str, str] = None,
    figsize: Tuple[int, int] = (12, 8),
    dpi: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7,
    metrics_text: Optional[str] = None,
    normalize: bool = True,
):
    """Create a publication-ready static matplotlib mirror plot."""
    if color_map is None:
        color_map = ION_COLORS

    fig = MplFigure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    if normalize:
        ann_offset = 0.03
    else:
        upper_max = upper_spectrum["intensity"].max() if not upper_spectrum.empty else 1
        lower_max = lower_spectrum["intensity"].max() if not lower_spectrum.empty else 1
        ann_offset = max(upper_max, lower_max) * 0.03

    if not upper_spectrum.empty:
        upper_df = _prepare_spectrum_data(upper_spectrum, normalize=normalize, color_map=color_map)
        _add_matplotlib_spectrum(ax, upper_df, y_multiplier=1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

    if not lower_spectrum.empty:
        lower_df = _prepare_spectrum_data(lower_spectrum, normalize=normalize, color_map=color_map)
        lower_df["intensity"] = lower_df["intensity"] * -1
        _add_matplotlib_spectrum(ax, lower_df, y_multiplier=-1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

    all_mz = []
    if not upper_spectrum.empty:
        all_mz.extend(upper_spectrum["mz"].tolist())
    if not lower_spectrum.empty:
        all_mz.extend(lower_spectrum["mz"].tolist())

    if all_mz:
        ax.set_xlim(min(all_mz) - 50, max(all_mz) + 50)

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
    color_map: Dict[str, str] = None,
    figsize: Tuple[int, int] = (12, 8),
    dpi: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7,
    normalize: bool = True,
):
    """Create a static matplotlib single annotated spectrum plot."""
    if color_map is None:
        color_map = ION_COLORS

    fig = MplFigure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    if not spectrum.empty:
        df = _prepare_spectrum_data(spectrum, normalize=normalize, color_map=color_map)
        if normalize:
            ann_offset = 0.03
        else:
            ann_offset = df["intensity"].max() * 0.03 if len(df) > 0 else 1
        _add_matplotlib_spectrum(ax, df, y_multiplier=1.0,
                                show_annotations=show_annotations,
                                annotation_fontsize=annotation_fontsize,
                                ann_offset=ann_offset)

    if not spectrum.empty:
        all_mz = spectrum["mz"].tolist()
        ax.set_xlim(min(all_mz) - 50, max(all_mz) + 50)

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
    buf = BytesIO()
    fig.savefig(buf, format=format, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.read()
