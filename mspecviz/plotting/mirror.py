"""Mirror + single-spectrum plotting via pyopenms-viz.

Replaces the previous hand-rolled matplotlib/plotly backends with a single
``pyopenms_viz.plot_spectrum`` call. This wrapper layers the MSpecViz-specific
decorations on top of the library's output: the Material ion-colour palette
(pyopenms-viz' built-in covers only a/b/c/x/y/z; we add internal / precursor /
immonium / unannotated), corner Query/Reference labels, and the metrics-overlay
text box. Intensities are always rendered relative (max-normalised) — pyopenms-viz
forces this for mirror plots and we apply the same convention to single plots
for visual consistency.
"""
from typing import Literal, Optional

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure as MplFigure  # noqa: E402
import pandas as pd  # noqa: E402  (kept after matplotlib.use)
import plotly.graph_objects as go  # noqa: E402
import pyopenms_viz  # noqa: E402

from mspecviz.common.defaults import ION_COLORS  # noqa: E402

Backend = Literal["ms_matplotlib", "ms_plotly"]


def _inject_color_column(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of *df* with ``annotation``, ``ion_type``, and ``_color`` populated.

    ``_color`` overrides pyopenms-viz' built-in palette via the ``peak_color``
    column kwarg so internal / precursor / immonium / unannotated peaks render
    in the MSpecViz Material colours.
    """
    df = df.copy()
    for col in ("annotation", "ion_type"):
        if col not in df.columns:
            df[col] = "unannotated"
    df["_color"] = df["ion_type"].map(ION_COLORS).fillna(ION_COLORS["unannotated"])
    return df


def _add_decorations_plotly(
    fig: go.Figure, *,
    query_label: str,
    reference_label: str,
    metrics_text: str,
) -> None:
    """Layer MSpecViz corner labels + metrics-overlay onto a Plotly figure."""
    if query_label:
        fig.add_annotation(
            x=0.02, y=0.98, xref="paper", yref="paper",
            text=f"<b>{query_label}</b>", showarrow=False,
            font=dict(size=11, color="#333"),
            bgcolor="rgba(255,255,255,0.8)",
            xanchor="left", yanchor="top",
        )
    if reference_label:
        fig.add_annotation(
            x=0.02, y=0.02, xref="paper", yref="paper",
            text=f"<b>{reference_label}</b>", showarrow=False,
            font=dict(size=11, color="#333"),
            bgcolor="rgba(255,255,255,0.8)",
            xanchor="left", yanchor="bottom",
        )
    if metrics_text:
        fig.add_annotation(
            x=0.98, y=0.98, xref="paper", yref="paper",
            text=metrics_text.replace("\n", "<br>"),
            showarrow=False,
            font=dict(size=10, family="monospace", color="#333"),
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="#CCC", borderwidth=1,
            xanchor="right", yanchor="top", align="left",
        )


def _add_decorations_matplotlib(
    ax, *,
    query_label: str,
    reference_label: str,
    metrics_text: str,
) -> None:
    """Layer MSpecViz corner labels + metrics-overlay onto a Matplotlib Axes."""
    label_kwargs = dict(
        transform=ax.transAxes, fontsize=9, fontweight="bold",
        ha="left", color="#333",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="none", alpha=0.8),
    )
    if query_label:
        ax.text(0.02, 0.99, query_label, va="top", **label_kwargs)
    if reference_label:
        ax.text(0.02, 0.01, reference_label, va="bottom", **label_kwargs)
    if metrics_text:
        ax.text(
            0.98, 0.98, metrics_text,
            transform=ax.transAxes, fontsize=10, fontfamily="monospace",
            va="top", ha="right", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#CCC", alpha=0.85),
        )


def create_spectrum_plot(
    query_frame: pd.DataFrame,
    reference_frame: Optional[pd.DataFrame] = None,
    *,
    backend: Backend = "ms_matplotlib",
    title: str = "",
    query_label: str = "",
    reference_label: str = "",
    metrics_text: str = "",
    figure_size: tuple[int, int] = (12, 8),
):
    """Render a mirror plot when *reference_frame* is supplied, otherwise a single plot.

    Returns ``plotly.graph_objects.Figure`` for ``backend="ms_plotly"`` and
    ``matplotlib.figure.Figure`` for ``backend="ms_matplotlib"``. *figure_size*
    is in inches at 100 dpi (pyopenms-viz takes pixels — the wrapper converts).
    """
    is_mirror = reference_frame is not None

    # pyopenms-viz raises on empty frames; return a blank figure of the right type
    # so callers can still treat the result uniformly.
    if query_frame.empty and (not is_mirror or reference_frame.empty):
        if backend == "ms_matplotlib":
            fig = MplFigure(figsize=figure_size)
            fig.add_subplot(111).set_title(title)
            return fig
        return go.Figure(layout=dict(title=title))

    query = _inject_color_column(query_frame)
    reference = _inject_color_column(reference_frame) if is_mirror else None

    result = pyopenms_viz.plot_spectrum(
        query, x="mz", y="intensity", backend=backend,
        reference_spectrum=reference,
        mirror_spectrum=is_mirror,
        ion_annotation="annotation",
        peak_color="_color",
        annotate_top_n_peaks="all",
        annotate_mz=False,
        relative_intensity=True,
        title=title,
        width=figure_size[0] * 100,
        height=figure_size[1] * 100,
        show_plot=False,
    )

    effective_reference_label = reference_label if is_mirror else ""
    if backend == "ms_matplotlib":
        _add_decorations_matplotlib(
            result,
            query_label=query_label,
            reference_label=effective_reference_label,
            metrics_text=metrics_text,
        )
        return result.figure

    _add_decorations_plotly(
        result,
        query_label=query_label,
        reference_label=effective_reference_label,
        metrics_text=metrics_text,
    )
    return result
