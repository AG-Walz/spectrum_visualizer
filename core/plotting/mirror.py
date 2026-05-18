"""Mirror + single-spectrum plotting (matplotlib for export, Plotly for the Streamlit Explore view); called by ``plotting.plot_worker`` and ``streamlit_app.views.explore.mirror_plots``."""
from typing import Literal, Optional

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure as MplFigure  # noqa: E402
import pandas as pd  # noqa: E402
import plotly.graph_objects as go  # noqa: E402
import pyopenms_viz  # noqa: F401, E402  # pyright: ignore[reportUnusedImport]  # registers the .plot accessor

from core.common.defaults import ION_COLORS  # noqa: E402

# Implementation note: the Plotly path delegates stem rendering to
# pyopenms-viz' ``df.plot(kind="spectrum")``, then ``_recolor_plotly_traces``
# overrides its palette (per-row colour columns fight its grouping logic).
# The matplotlib path stays hand-rolled because pyopenms-viz' matplotlib
# backend has a colour-scrambling bug in
# ``_matplotlib/core.py`` (``groupby(self.by, sort=True)``) — peaks of
# interleaved ion types get the wrong colour, which is "always" in real
# spectra.

Backend = Literal["ms_matplotlib", "ms_plotly"]

# Vertical offset for matplotlib peak labels (fraction of normalised intensity).
_LABEL_OFFSET = 0.03
# Extra m/z padding on the x-axis (Da) so outer peaks aren't clipped.
_MZ_PADDING_DA = 50
# pyopenms-viz forces intensities to [0, 100] when mirror_spectrum=True or
# relative_intensity=True. We always use that path, so labels go in the same scale.
_PLOTLY_INTENSITY_SCALE = 100


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    """Copy ``df`` with per-peak colour + label columns added and intensity rescaled to ``[0, 1]``."""
    df = df.copy()
    for col in ("annotation", "ion_type"):
        if col not in df.columns:
            df[col] = "unannotated"
    df["_color"] = df["ion_type"].map(ION_COLORS).fillna(ION_COLORS["unannotated"])
    df["_label"] = df["annotation"].where(df["annotation"] != "unannotated", "")
    if not df.empty:
        max_intensity = df["intensity"].max()
        if max_intensity and max_intensity > 0:
            df["intensity"] = df["intensity"] / max_intensity
    return df


def _draw_matplotlib(ax, df: pd.DataFrame, y_multiplier: float) -> None:
    """Draw the peak stems and per-peak text labels for one spectrum onto ``ax`` (use ``y_multiplier=-1`` to flip a reference spectrum downwards)."""
    if df.empty:
        return
    # One vlines call per ion type so each type gets a single, correctly-coloured legend entry
    for ion_type, group in df.groupby("ion_type", sort=False):
        ax.vlines(
            group["mz"].values,
            0, group["intensity"].values * y_multiplier,
            colors=ION_COLORS.get(ion_type, ION_COLORS["unannotated"]),
            linewidth=1.5, alpha=0.85,
            label=ion_type if ion_type != "unannotated" else None,
        )
    sign = 1 if y_multiplier > 0 else -1
    va = "bottom" if y_multiplier > 0 else "top"
    for _, peak in df[df["_label"] != ""].iterrows():
        ax.text(
            peak["mz"],
            peak["intensity"] * y_multiplier + sign * _LABEL_OFFSET,
            peak["_label"],
            fontsize=7, color=peak["_color"],
            ha="center", va=va,
        )


def _recolor_plotly_traces(fig: go.Figure) -> None:
    """Recolour every pyopenms-viz trace using the MSpecViz ``ION_COLORS`` palette and keep only one legend entry per ion type."""
    legend_seen: set[str] = set()
    for trace in fig.data:
        ion_type = trace.name
        if ion_type in ION_COLORS:
            trace.line.color = ION_COLORS[ion_type]
        if ion_type == "unannotated" or ion_type in legend_seen:
            trace.showlegend = False
        else:
            trace.showlegend = True
            trace.legendgroup = ion_type
            legend_seen.add(ion_type)


def _add_peak_labels_plotly(fig: go.Figure, df: pd.DataFrame,
                            y_multiplier: float) -> None:
    """Place per-peak text annotations next to each stem, coloured to match the peak's ion type."""
    if df.empty:
        return
    yshift = 10 if y_multiplier > 0 else -10
    for _, peak in df[df["_label"] != ""].iterrows():
        fig.add_annotation(
            x=peak["mz"],
            y=peak["intensity"] * y_multiplier * _PLOTLY_INTENSITY_SCALE,
            text=peak["_label"], showarrow=False,
            font=dict(size=8, color=peak["_color"]),
            yshift=yshift,
        )


def _set_axes_matplotlib(ax, query_df: pd.DataFrame,
                         reference_df: Optional[pd.DataFrame],
                         is_mirror: bool) -> None:
    """Set x/y limits, ticks, bold m/z + intensity labels, and the mirror divider line."""
    # X-axis: pad min/max m/z across both panels so outer peaks aren't clipped
    mz_values = []
    if not query_df.empty:
        mz_values.extend(query_df["mz"].tolist())
    if reference_df is not None and not reference_df.empty:
        mz_values.extend(reference_df["mz"].tolist())
    if mz_values:
        ax.set_xlim(min(mz_values) - _MZ_PADDING_DA, max(mz_values) + _MZ_PADDING_DA)

    # Y-axis: symmetric ticks + divider for mirror, single-panel scale otherwise
    if is_mirror:
        ax.set_ylim(-1.2, 1.2)
        ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
        ax.set_yticklabels(["1.0", "0.5", "0", "0.5", "1.0"])
        ax.axhline(y=0, color="black", linewidth=1, alpha=0.5)
    else:
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    # Labels + spines: bold m/z & intensity, drop top/right spines for a clean frame
    ax.set_xlabel("m/z", fontsize=12, fontweight="bold")
    ax.set_ylabel("Relative Intensity", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _add_legend_matplotlib(ax) -> None:
    """Add a single deduplicated ion-type legend in the bottom-right (top-right is reserved for the metrics overlay box)."""
    handles, labels = ax.get_legend_handles_labels()
    if not labels:
        return
    seen = set()
    unique = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
    ax.legend([h for h, _ in unique], [l for _, l in unique],
              loc="lower right", fontsize=9, frameon=True,
              title="Ion type", title_fontsize=9)


def _add_decorations_plotly(fig: go.Figure, *, query_label: str,
                            reference_label: str, metrics_text: str) -> None:
    """Stamp the query/reference corner labels and the metrics overlay box onto a Plotly figure."""
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


def _add_decorations_matplotlib(ax, *, query_label: str,
                                reference_label: str, metrics_text: str) -> None:
    """Stamp the query/reference corner labels and the metrics overlay box onto a matplotlib Axes."""
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
    """Render a mirror plot when ``reference_frame`` is given, otherwise a single spectrum; returns a matplotlib or Plotly figure depending on ``backend``, with intensities max-normalised for visual consistency."""
    # Prep: normalise both panels, suppress reference label for single-spectrum
    is_mirror = reference_frame is not None
    query = _prep(query_frame)
    reference = _prep(reference_frame) if is_mirror else None
    effective_reference_label = reference_label if is_mirror else ""

    # Matplotlib path: draw stems, axes, legend, corner overlays, then tighten
    if backend == "ms_matplotlib":
        fig = MplFigure(figsize=figure_size, dpi=100)
        ax = fig.add_subplot(111)
        _draw_matplotlib(ax, query, y_multiplier=1.0)
        if reference is not None:
            _draw_matplotlib(ax, reference, y_multiplier=-1.0)
        _set_axes_matplotlib(ax, query, reference, is_mirror)
        ax.set_title(title, fontsize=14, fontweight="bold")
        _add_legend_matplotlib(ax)
        _add_decorations_matplotlib(
            ax, query_label=query_label,
            reference_label=effective_reference_label,
            metrics_text=metrics_text,
        )
        fig.tight_layout()
        return fig

    # Plotly path: delegate stem rendering to pyopenms-viz, then recolour + decorate
    fig: go.Figure
    if query.empty and (reference is None or reference.empty):
        # pyopenms-viz raises on empty data; return a blank figure instead.
        fig = go.Figure()
    else:
        # pandas .plot stub doesn't know the pyopenms-viz "spectrum" kind.
        fig = query.plot(
            kind="spectrum", backend="ms_plotly",  # type: ignore[arg-type]
            x="mz", y="intensity",
            ion_annotation="ion_type",
            mirror_spectrum=is_mirror,
            reference_spectrum=reference if is_mirror else None,
            relative_intensity=True,
            annotate_top_n_peaks=0,    # we add labels ourselves below
            grid=False,
            show_plot=False,
            title=title,
        )
        _recolor_plotly_traces(fig)
        _add_peak_labels_plotly(fig, query, y_multiplier=1.0)
        if reference is not None:
            _add_peak_labels_plotly(fig, reference, y_multiplier=-1.0)

    # Layout: legend + dimensions + axis ranges (overrides pyopenms-viz' own padding)
    fig.update_layout(
        title=title, xaxis_title="m/z", yaxis_title="Relative Intensity",
        width=figure_size[0] * 100, height=figure_size[1] * 100,
        template="plotly_white", hovermode="closest",
        # Bottom-right; top-right is reserved for the metrics overlay box.
        legend=dict(
            title="Ion type",
            yanchor="bottom", y=0.02,
            xanchor="right", x=0.98,
            bgcolor="rgba(255,255,255,0.85)", bordercolor="#CCC", borderwidth=1,
        ),
    )
    if is_mirror:
        fig.update_yaxes(
            range=[-_PLOTLY_INTENSITY_SCALE * 1.25, _PLOTLY_INTENSITY_SCALE * 1.25],
            tickvals=[-100, -50, 0, 50, 100],
            ticktext=["1.0", "0.5", "0", "0.5", "1.0"],
        )
    else:
        fig.update_yaxes(range=[0, _PLOTLY_INTENSITY_SCALE * 1.15])
    mz_values = []
    if not query.empty:
        mz_values.extend(query["mz"].tolist())
    if is_mirror and reference is not None and not reference.empty:
        mz_values.extend(reference["mz"].tolist())
    if mz_values:
        fig.update_xaxes(range=[min(mz_values) - _MZ_PADDING_DA,
                                max(mz_values) + _MZ_PADDING_DA])
    _add_decorations_plotly(
        fig, query_label=query_label,
        reference_label=effective_reference_label,
        metrics_text=metrics_text,
    )
    return fig
