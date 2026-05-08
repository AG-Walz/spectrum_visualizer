"""Interactive Plotly mirror + single-spectrum plots.

Used for Streamlit display and HTML export. Shares data-prep with the
Matplotlib backend via ``_mirror_shared``.
"""
import pandas as pd
import plotly.graph_objects as go

from mspecviz.plotting._mirror_shared import (
    MZ_WINDOW_PADDING,
    _get_annotated_mask,
    _prepare_spectrum_data,
)


def _add_plotly_spectrum_traces(fig, df: pd.DataFrame, y_multiplier: float, label: str):
    """Add Plotly stem traces and per-peak annotation labels for one spectrum.

    ``y_multiplier`` is +1 for upper/single, -1 for lower (mirror).
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


def _create_spectrum_figure(
    spectra: list[tuple[pd.DataFrame, str, float]],
    *, title: str, normalize: bool, height: int,
    mirror: bool, metrics_text: str | None = None,
) -> go.Figure:
    """Build a Plotly figure; ``spectra`` is ``[(df, label, y_multiplier)]``."""
    fig = go.Figure()
    for df, label, y_mult in spectra:
        if not df.empty:
            _add_plotly_spectrum_traces(
                fig, _prepare_spectrum_data(df, normalize=normalize),
                y_multiplier=y_mult, label=label,
            )
        if label:
            _add_corner_label(fig, label, top=(y_mult > 0))
    if mirror:
        _add_mirror_center_line(fig, spectra[0][0], spectra[1][0])
    _apply_spectrum_layout(fig, title, normalize, height)
    if normalize and mirror:
        fig.update_yaxes(range=[-1.25, 1.25],
                         tickvals=[-1.0, -0.5, 0, 0.5, 1.0],
                         ticktext=["1.0", "0.5", "0", "0.5", "1.0"])
    elif normalize:
        fig.update_yaxes(range=[0, 1.15])
    if metrics_text:
        _add_metrics_overlay(fig, metrics_text)
    return fig


def create_plotly_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    metrics_text: str | None = None,
    normalize: bool = True,
    height: int = 600,
) -> go.Figure:
    """Render an interactive Plotly mirror plot (query on top, reference below)."""
    return _create_spectrum_figure(
        [(upper_spectrum, upper_label, 1.0), (lower_spectrum, lower_label, -1.0)],
        title=title, normalize=normalize, height=height,
        mirror=True, metrics_text=metrics_text,
    )


def create_plotly_single_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    normalize: bool = True,
    height: int = 500,
) -> go.Figure:
    """Create an interactive Plotly single-spectrum plot."""
    return _create_spectrum_figure(
        [(spectrum, label, 1.0)],
        title=title, normalize=normalize, height=height,
        mirror=False,
    )
