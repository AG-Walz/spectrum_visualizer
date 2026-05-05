"""Interactive Plotly mirror + single-spectrum plots.

Used for Streamlit display and HTML export. Shares data-prep with the
Matplotlib backend via ``_mirror_shared``.
"""
from typing import Optional

import pandas as pd
import plotly.graph_objects as go

from mspecviz.plotting._mirror_shared import (
    MZ_WINDOW_PADDING,
    _get_annotated_mask,
    _prepare_spectrum_data,
)


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


def _add_mirror_center_line(fig, upper_spectrum: pd.DataFrame,
                            lower_spectrum: pd.DataFrame) -> None:
    """Draw the y=0 horizontal divider across the m/z range of both spectra."""
    all_mz = []
    if not upper_spectrum.empty:
        all_mz.extend(upper_spectrum["mz"].tolist())
    if not lower_spectrum.empty:
        all_mz.extend(lower_spectrum["mz"].tolist())
    if not all_mz:
        return
    fig.add_shape(
        type="line",
        x0=min(all_mz) - MZ_WINDOW_PADDING,
        x1=max(all_mz) + MZ_WINDOW_PADDING,
        y0=0, y1=0,
        line=dict(color="black", width=1, dash="solid"),
    )


def _apply_spectrum_layout(fig, title: str, normalize: bool, height: int) -> None:
    """Apply the shared axis/title/template layout used by both Plotly plots."""
    fig.update_layout(
        title=title,
        xaxis_title="m/z",
        yaxis_title="Relative Intensity" if normalize else "Intensity",
        height=height,
        hovermode="closest",
        template="plotly_white",
        showlegend=False,
    )


def _add_corner_label(fig, text: str, *, top: bool) -> None:
    """Add a bold corner label (top-left or bottom-left) in paper coordinates."""
    y, yanchor = (0.98, "top") if top else (0.02, "bottom")
    fig.add_annotation(
        x=0.02, y=y, xref="paper", yref="paper",
        text=f"<b>{text}</b>", showarrow=False,
        font=dict(size=11, color="#333"),
        bgcolor="rgba(255,255,255,0.8)",
        xanchor="left", yanchor=yanchor,
    )


def _add_metrics_overlay(fig, metrics_text: str) -> None:
    """Add the upper-right metrics text box."""
    fig.add_annotation(
        x=0.98, y=0.98, xref="paper", yref="paper",
        text=metrics_text.replace("\n", "<br>"),
        showarrow=False,
        font=dict(size=10, family="monospace", color="#333"),
        bgcolor="rgba(255,255,255,0.85)",
        bordercolor="#CCC", borderwidth=1,
        xanchor="right", yanchor="top", align="left",
    )


def create_plotly_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    metrics_text: Optional[str] = None,
    normalize: bool = True,
    height: int = 600,
):
    """Render an interactive Plotly mirror plot (query on top, reference below)."""
    fig = go.Figure()

    if not upper_spectrum.empty:
        upper_frame = _prepare_spectrum_data(upper_spectrum, normalize=normalize)
        _add_plotly_spectrum_traces(fig, upper_frame, y_multiplier=1.0,
                                    label=upper_label, show_annotations=True)
    if not lower_spectrum.empty:
        lower_frame = _prepare_spectrum_data(lower_spectrum, normalize=normalize)
        _add_plotly_spectrum_traces(fig, lower_frame, y_multiplier=-1.0,
                                    label=lower_label, show_annotations=True)

    _add_mirror_center_line(fig, upper_spectrum, lower_spectrum)
    _apply_spectrum_layout(fig, title, normalize, height)
    if normalize:
        fig.update_yaxes(range=[-1.25, 1.25],
                         tickvals=[-1.0, -0.5, 0, 0.5, 1.0],
                         ticktext=["1.0", "0.5", "0", "0.5", "1.0"])

    _add_corner_label(fig, upper_label, top=True)
    _add_corner_label(fig, lower_label, top=False)
    if metrics_text:
        _add_metrics_overlay(fig, metrics_text)
    return fig


def create_plotly_single_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    normalize: bool = True,
    height: int = 500,
):
    """Create an interactive Plotly single-spectrum plot."""
    fig = go.Figure()
    if not spectrum.empty:
        df = _prepare_spectrum_data(spectrum, normalize=normalize)
        _add_plotly_spectrum_traces(fig, df, y_multiplier=1.0,
                                    label=label, show_annotations=True)
    _apply_spectrum_layout(fig, title, normalize, height)
    if normalize:
        fig.update_yaxes(range=[0, 1.15])
    if label:
        _add_corner_label(fig, label, top=True)
    return fig
