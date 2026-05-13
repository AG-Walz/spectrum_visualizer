"""
openms_insight Table / SequenceView / LinePlot / Heatmap fragment.

The four components share a single ``StateManager`` so row selections in the
Table drive the filters for the other three. Cache files are stored under the
configured export directory and regenerated whenever the underlying analysis
result changes (``cache_version`` bump).
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import Heatmap, LinePlot, SequenceView, StateManager, Table

from streamlit_app.utils.state import get_param, get_state_key
from core.common.defaults import DEFAULTS


def _create_component(cls, cache_id, cache_path, needs_regeneration, **kwargs):
    """Create an openms_insight component, using cache when available."""
    if needs_regeneration:
        return cls(cache_id=cache_id, cache_path=cache_path, **kwargs)
    return cls(cache_id=cache_id, cache_path=cache_path)


def _needs_cache_regen():
    """Compare this fragment's last-rendered cache version against the current one.

    Returns ``(needs_regen, mark_rendered)``. Call ``mark_rendered()`` after a
    successful render so the next rerun can short-circuit. Centralising this
    here keeps the ``cache_version`` / ``cache_rendered_version`` session-state
    key pair in one place.
    """
    cache_version = st.session_state.get(get_state_key("cache_version"), 0)
    reg_key = get_state_key("cache_rendered_version")
    needs = st.session_state.get(reg_key, -1) != cache_version

    def mark_rendered():
        st.session_state[reg_key] = cache_version

    return needs, mark_rendered


def _resolve_cache_path(params: dict) -> str:
    """Return the explorer's component-cache directory, creating it if needed."""
    export_dir = get_param(params, "export_directory", "")
    cache_path = export_dir if export_dir else str(Path(st.session_state.workspace) / "exports")
    os.makedirs(cache_path, exist_ok=True)
    return cache_path


def _get_state_manager() -> StateManager:
    """Return the cross-component StateManager, creating it on first access."""
    if "mspecviz_sm" not in st.session_state:
        st.session_state["mspecviz_sm"] = StateManager(session_key="mspecviz_state")
    return st.session_state["mspecviz_sm"]


def _show_selected_psm_banner(psm_frame: pl.LazyFrame | None, state_manager: StateManager) -> None:
    """Show an info banner naming the currently selected PSM, if any."""
    if psm_frame is None:
        return
    selected_id = state_manager.get_state_for_vue().get("selected_psm")
    if selected_id is None:
        return
    row = psm_frame.filter(pl.col("id") == selected_id).collect()
    if row.height == 0:
        return
    label = row["display_label"][0]
    source = row["source_file"][0]
    st.info(f"Viewing: **{label}** (PSM {selected_id}) from {source}")


def _render_psm_table(psm_frame, cache_path, needs_regeneration, state_manager) -> None:
    """Render the PSM Table; selecting a row updates ``selected_psm``."""
    st.subheader("PSM Table")
    st.caption("Click a row to update the spectrum and sequence views below.")
    columns = [
        {"field": "display_label", "title": "Peptide"},
        {"field": "charge", "title": "Z"},
        {"field": "score", "title": "Score", "formatter": "money", "formatterParams": {"precision": 4}},
        {"field": "source_file", "title": "Source File"},
    ]
    table = _create_component(
        Table, "psm_table", cache_path, needs_regeneration,
        data=psm_frame, interactivity={"selected_psm": "id"},
        column_definitions=columns, regenerate_cache=True, index_field="id",
    )
    table(key="psm_table", state_manager=state_manager)


def _render_sequence_view(
    sequence_frame, peaks_frame, params, cache_path, needs_regeneration, state_manager,
) -> None:
    """Render the SequenceView coverage map filtered by ``selected_psm``."""
    st.subheader("Sequence Coverage")
    if needs_regeneration:
        peaks_for_sv = peaks_frame.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
    else:
        peaks_for_sv = None
    sequence_view = _create_component(
        SequenceView, "sequence_view", cache_path, needs_regeneration,
        sequence_data=sequence_frame, peaks_data=peaks_for_sv,
        filters={"selected_psm": "psm_id"},
        annotation_config={
            "ion_types": get_param(params, "ion_types", DEFAULTS.ion_types),
            "tolerance": get_param(params, "mz_tolerance", DEFAULTS.mz_tolerance),
            "tolerance_ppm": False,
        },
        height=350,
    )
    sequence_view(key="sequence_view", state_manager=state_manager)


def _render_line_plot(peaks_frame, cache_path, needs_regeneration, state_manager) -> None:
    """Render the fragment-spectrum LinePlot, linked to the SequenceView."""
    st.subheader("Fragment Spectrum")
    lineplot = _create_component(
        LinePlot, "spectrum_plot", cache_path, needs_regeneration,
        data=peaks_frame, filters={"selected_psm": "psm_id"},
        filter_defaults={"selected_psm": 0}, regenerate_cache=True,
        x_column="x", y_column="y", highlight_column="highlight",
        annotation_column="annotation", x_label="m/z", y_label="Intensity",
    )
    lineplot(key="spectrum_plot", state_manager=state_manager, sequence_view_key="sequence_view")


def _render_heatmap(heatmap_frame, cache_path, needs_regeneration, state_manager) -> None:
    """Render the RT x m/z Heatmap (log-scaled intensity)."""
    st.subheader("RT x m/z Heatmap")
    heatmap = _create_component(
        Heatmap, "rt_heatmap", cache_path, needs_regeneration,
        data=heatmap_frame, x_column="x", y_column="y", intensity_column="intensity",
        x_label="Retention Time (s)", y_label="m/z", log_scale=True, regenerate_cache=True,
    )
    heatmap(key="rt_heatmap", state_manager=state_manager)


@st.fragment
def render_explorer(params: dict) -> None:
    """openms_insight Table + LinePlot + SequenceView + Heatmap."""
    # Setup: cache directory, regen flag, and the shared StateManager
    cache_path = _resolve_cache_path(params)
    needs_regeneration, mark_rendered = _needs_cache_regen()
    state_manager = _get_state_manager()

    # Pull frames: four LazyFrames staged by Search & Load into session state
    psm_frame = st.session_state.get(get_state_key("psm_frame"))
    peaks_frame = st.session_state.get(get_state_key("peaks_frame"))
    sequence_frame = st.session_state.get(get_state_key("sequence_frame"))
    heatmap_frame = st.session_state.get(get_state_key("heatmap_frame"))

    # Render components: banner + table + sequence view + line plot + heatmap
    _show_selected_psm_banner(psm_frame, state_manager)
    _render_psm_table(psm_frame, cache_path, needs_regeneration, state_manager)
    _render_sequence_view(
        sequence_frame, peaks_frame, params, cache_path, needs_regeneration, state_manager,
    )
    _render_line_plot(peaks_frame, cache_path, needs_regeneration, state_manager)
    _render_heatmap(heatmap_frame, cache_path, needs_regeneration, state_manager)

    # Mark rendered: skip cache regeneration on subsequent reruns
    if needs_regeneration:
        mark_rendered()