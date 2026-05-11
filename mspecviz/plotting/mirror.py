"""Mirror + single-spectrum plotting (matplotlib + Plotly).

Custom hand-rolled renderer because pyopenms-viz' ``plot_spectrum`` cannot
honour per-peak ``peak_color`` correctly: it groups peaks by a category
column with ``sort=True`` and consumes one colour per group from a generator
iterating the row-aligned colour column — so the first peak's colour lands
on the first alphabetical group regardless of which ion type that group is.
Result: scrambled b/y/internal/etc. colouring whenever ion types are
interleaved in m/z order (i.e. always, in real spectra).

We draw the peaks ourselves with vlines / line scatters, grouped by
``ion_type`` so the legend stays one entry per type with the canonical
``ION_COLORS`` palette. Labels are coloured to match the peak.
"""
from typing import Literal, Optional

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure as MplFigure  # noqa: E402
import pandas as pd  # noqa: E402
import plotly.graph_objects as go  # noqa: E402

from mspecviz.common.defaults import ION_COLORS  # noqa: E402

Backend = Literal["ms_matplotlib", "ms_plotly"]

# Vertical offset for peak labels, as a fraction of the normalised intensity range.
_LABEL_OFFSET = 0.03
# Extra m/z padding on the plot x-axis (Da) so outer peaks aren't clipped.
_MZ_PADDING_DA = 50


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``_color`` + ``_label`` columns and normalise intensity to [0, 1]."""
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
    """Draw stems + per-peak labels for one spectrum onto ``ax``."""
    if df.empty:
        return
    # One vlines call per ion type → one legend entry per type, correct colour.
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


def _draw_plotly(fig: go.Figure, df: pd.DataFrame, y_multiplier: float,
                 panel_name: str) -> None:
    """Draw stems (one trace per ion type) + per-peak labels onto ``fig``."""
    if df.empty:
        return
    for ion_type, group in df.groupby("ion_type", sort=False):
        color = ION_COLORS.get(ion_type, ION_COLORS["unannotated"])
        x_data, y_data, hover_data = [], [], []
        for _, peak in group.iterrows():
            x_data.extend([peak["mz"], peak["mz"], None])
            y_data.extend([0, peak["intensity"] * y_multiplier, None])
            hover_data.extend([
                "",
                f"m/z: {peak['mz']:.4f}<br>Intensity: {peak['intensity']:.3f}"
                f"<br>{peak['annotation'] if peak['annotation'] != 'unannotated' else ''}",
                "",
            ])
        fig.add_trace(go.Scatter(
            x=x_data, y=y_data, mode="lines",
            line=dict(color=color, width=1.5),
            name=ion_type,
            legendgroup=f"{panel_name}_{ion_type}",
            hovertext=hover_data, hoverinfo="text",
            showlegend=(ion_type != "unannotated"),
        ))
    yshift = 10 if y_multiplier > 0 else -10
    for _, peak in df[df["_label"] != ""].iterrows():
        fig.add_annotation(
            x=peak["mz"], y=peak["intensity"] * y_multiplier,
            text=peak["_label"], showarrow=False,
            font=dict(size=8, color=peak["_color"]),
            yshift=yshift,
        )


def _set_axes_matplotlib(ax, query_df: pd.DataFrame,
                         reference_df: Optional[pd.DataFrame],
                         is_mirror: bool) -> None:
    """Apply x/y limits + axis labels + mirror divider line."""
    mz_values = []
    if not query_df.empty:
        mz_values.extend(query_df["mz"].tolist())
    if reference_df is not None and not reference_df.empty:
        mz_values.extend(reference_df["mz"].tolist())
    if mz_values:
        ax.set_xlim(min(mz_values) - _MZ_PADDING_DA, max(mz_values) + _MZ_PADDING_DA)

    if is_mirror:
        ax.set_ylim(-1.2, 1.2)
        ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
        ax.set_yticklabels(["1.0", "0.5", "0", "0.5", "1.0"])
        ax.axhline(y=0, color="black", linewidth=1, alpha=0.5)
    else:
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    ax.set_xlabel("m/z", fontsize=12, fontweight="bold")
    ax.set_ylabel("Relative Intensity", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _add_legend_matplotlib(ax) -> None:
    """Dedupe legend entries so each ion type appears once. Bottom-right
    (top-right is reserved for the metrics overlay box)."""
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


def _add_decorations_matplotlib(ax, *, query_label: str,
                                reference_label: str, metrics_text: str) -> None:
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
    """Render a mirror plot when *reference_frame* is supplied, otherwise single.

    Returns ``matplotlib.figure.Figure`` for ``ms_matplotlib`` and
    ``plotly.graph_objects.Figure`` for ``ms_plotly``. Intensities are always
    rendered relative (max-normalised) for visual consistency between the
    two panels in a mirror plot.
    """
    is_mirror = reference_frame is not None
    query = _prep(query_frame)
    reference = _prep(reference_frame) if is_mirror else None
    effective_reference_label = reference_label if is_mirror else ""

    if backend == "ms_matplotlib":
        fig = MplFigure(figsize=figure_size, dpi=100)
        ax = fig.add_subplot(111)
        _draw_matplotlib(ax, query, y_multiplier=1.0)
        if is_mirror:
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

    fig = go.Figure()
    _draw_plotly(fig, query, y_multiplier=1.0, panel_name="query")
    if is_mirror:
        _draw_plotly(fig, reference, y_multiplier=-1.0, panel_name="reference")
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
            range=[-1.25, 1.25],
            tickvals=[-1.0, -0.5, 0, 0.5, 1.0],
            ticktext=["1.0", "0.5", "0", "0.5", "1.0"],
        )
        fig.add_shape(
            type="line", xref="paper", x0=0, x1=1, y0=0, y1=0,
            line=dict(color="black", width=1),
        )
    else:
        fig.update_yaxes(range=[0, 1.15])
    _add_decorations_plotly(
        fig, query_label=query_label,
        reference_label=effective_reference_label,
        metrics_text=metrics_text,
    )
    return fig
