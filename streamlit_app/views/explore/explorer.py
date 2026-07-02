"""Interactive Explorer sub-view of the Explore tab — Table + SequenceView + LinePlot + Heatmap.
Rendered by explore_page.py; one shared StateManager cross-links the four components.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import Heatmap, LinePlot, SequenceView, StateManager, Table

from streamlit_app.utils.state import get_param, get_state_key
from core.defaults import DEFAULTS


def _create_component(cls, cache_id, cache_path, needs_regeneration, **kwargs):
    """Instantiate an openms_insight component, skipping kwargs when the disk cache is fresh.
    """
    if needs_regeneration:
        return cls(cache_id=cache_id, cache_path=cache_path, **kwargs)
    return cls(cache_id=cache_id, cache_path=cache_path)


def _needs_cache_regen():
    """Return (needs_regen, mark_rendered) by comparing rendered vs current cache_version.
    Call mark_rendered() after a successful render so the next rerun short-circuits.
    """
    cache_version = st.session_state.get(get_state_key("cache_version"), 0)
    rendered_version_state_key = get_state_key("cache_rendered_version")
    needs_regeneration = st.session_state.get(rendered_version_state_key, -1) != cache_version

    def mark_rendered():
        st.session_state[rendered_version_state_key] = cache_version

    return needs_regeneration, mark_rendered


def _resolve_cache_path(params: dict) -> str:
    """Return the openms_insight component-cache directory, creating it if missing.
    """
    export_dir = get_param(params, "export_directory", "")
    cache_path = export_dir if export_dir else str(Path(st.session_state.workspace) / "exports")
    os.makedirs(cache_path, exist_ok=True)
    return cache_path


def _get_state_manager() -> StateManager:
    """Lazily build the cross-component StateManager and stash it in session state.
    """
    if "explorer_state_manager" not in st.session_state:
        st.session_state["explorer_state_manager"] = StateManager(session_key="mspecviz_state")
    return st.session_state["explorer_state_manager"]


def _show_selected_psm_banner(psm_frame: pl.LazyFrame | None, state_manager: StateManager) -> None:
    """Display an info banner with the label and source file of the active PSM, if one is picked.
    """
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
    """Render the openms_insight PSM Table; row clicks publish selected_psm to the StateManager.
    """
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
    """Render the SequenceView ion-coverage map for the active PSM.
    """
    st.subheader("Sequence Coverage")
    # Rename peaks columns: SequenceView expects peak_id/mass/intensity
    if needs_regeneration:
        peaks_for_sequence_view = peaks_frame.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
    else:
        peaks_for_sequence_view = None
    sequence_view = _create_component(
        SequenceView, "sequence_view", cache_path, needs_regeneration,
        sequence_data=sequence_frame, peaks_data=peaks_for_sequence_view,
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
    """Render the fragment-spectrum LinePlot, color-linked to the SequenceView.
    """
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
    """Render the retention-time x m/z Heatmap with log-scaled intensity.
    """
    st.subheader("RT x m/z Heatmap")
    heatmap = _create_component(
        Heatmap, "rt_heatmap", cache_path, needs_regeneration,
        data=heatmap_frame, x_column="x", y_column="y", intensity_column="intensity",
        x_label="Retention Time (s)", y_label="m/z", log_scale=True, regenerate_cache=True,
    )
    heatmap(key="rt_heatmap", state_manager=state_manager)


@st.fragment
def render_explorer(params: dict) -> None:
    """Render the full Interactive Explorer — banner + Table + SequenceView + LinePlot + Heatmap.
    """
    # Prepare: cache directory, regen flag, and the shared StateManager
    cache_path = _resolve_cache_path(params)
    needs_regeneration, mark_rendered = _needs_cache_regen()
    state_manager = _get_state_manager()

    # Pull frames: four LazyFrames staged by Load & Annotate into session state
    psm_frame = st.session_state.get(get_state_key("psm_frame"))
    peaks_frame = st.session_state.get(get_state_key("peaks_frame"))
    sequence_frame = st.session_state.get(get_state_key("sequence_frame"))
    heatmap_frame = st.session_state.get(get_state_key("heatmap_frame"))

    # Render: banner + table (always), then the per-PSM spectrum views, then the heatmap
    _show_selected_psm_banner(psm_frame, state_manager)
    _render_psm_table(psm_frame, cache_path, needs_regeneration, state_manager)

    # The SequenceView + LinePlot are per-PSM; a PSM with no annotated peaks (e.g. a reference
    # whose experimental spectrum was never loaded) would render as empty plots. Skip them and
    # their headers for such a PSM, showing a note instead — the RT x m/z Heatmap is global, so
    # it stays (#94). Always render on a regeneration pass so the disk cache is built first.
    selected_psm = state_manager.get_state_for_vue().get("selected_psm")
    selected_has_peaks = peaks_frame is not None and selected_psm is not None and (
        peaks_frame.filter(pl.col("psm_id") == selected_psm).head(1).collect().height > 0
    )
    if needs_regeneration or selected_psm is None or selected_has_peaks:
        _render_sequence_view(
            sequence_frame, peaks_frame, params, cache_path, needs_regeneration, state_manager,
        )
        _render_line_plot(peaks_frame, cache_path, needs_regeneration, state_manager)
    else:
        st.info("This spectrum has no annotated fragment peaks — only the RT x m/z heatmap below is available.")

    _render_heatmap(heatmap_frame, cache_path, needs_regeneration, state_manager)

    # Mark rendered: subsequent reruns reuse the parquet caches instead of rebuilding
    if needs_regeneration:
        mark_rendered()