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

import streamlit as st
import polars as pl
from openms_insight import Heatmap, LinePlot, SequenceView, StateManager, Table

from streamlit_app.utils.state import get_param, get_state_key
from mspecviz.common.defaults import DEFAULTS


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


@st.fragment
def openms_insight_explorer_fragment(params: dict):
    """openms_insight Table + LinePlot + SequenceView + Heatmap."""
    export_dir = get_param(params, "export_directory", "")
    CACHE_PATH = export_dir if export_dir else str(Path(st.session_state.workspace) / "exports")
    os.makedirs(CACHE_PATH, exist_ok=True)

    needs_regeneration, mark_rendered = _needs_cache_regen()

    if "mspecviz_sm" not in st.session_state:
        st.session_state["mspecviz_sm"] = StateManager(session_key="mspecviz_state")
    state_manager = st.session_state["mspecviz_sm"]

    psm_lf = st.session_state.get(get_state_key("psm_lf"))
    peaks_lf = st.session_state.get(get_state_key("peaks_lf"))
    seq_lf = st.session_state.get(get_state_key("seq_lf"))
    heatmap_lf = st.session_state.get(get_state_key("heatmap_lf"))

    # Show currently selected PSM info
    current_state = state_manager.get_state_for_vue()
    selected_id = current_state.get("selected_psm")
    if selected_id is not None and psm_lf is not None:
        try:
            row = psm_lf.filter(pl.col("id") == selected_id).collect()
            if row.height > 0:
                label = row["display_label"][0]
                source = row["source_file"][0]
                st.info(f"Viewing: **{label}** (PSM {selected_id}) from {source}")
        except Exception:
            pass

    # PSM Table
    st.subheader("PSM Table")
    st.caption("Click a row to update the spectrum and sequence views below.")
    psm_columns = [
        {"field": "display_label", "title": "Peptide"},
        {"field": "charge", "title": "Z"},
        {"field": "score", "title": "Score", "formatter": "money", "formatterParams": {"precision": 4}},
        {"field": "source_file", "title": "Source File"},
    ]
    table = _create_component(Table, "psm_table", CACHE_PATH, needs_regeneration,
                              data=psm_lf, interactivity={"selected_psm": "id"},
                              column_definitions=psm_columns, regenerate_cache=True, index_field="id")
    table(key="psm_table", state_manager=state_manager)

    # SequenceView
    st.subheader("Sequence Coverage")
    if needs_regeneration:
        peaks_for_sv = peaks_lf.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
    else:
        peaks_for_sv = None
    sequence_view = _create_component(
        SequenceView, "seq_view", CACHE_PATH, needs_regeneration,
        sequence_data=seq_lf, peaks_data=peaks_for_sv,
        filters={"selected_psm": "psm_id"},
        annotation_config={
            "ion_types": get_param(params, "ion_types", DEFAULTS.ion_types),
            "tolerance": get_param(params, "mz_tolerance", DEFAULTS.mz_tolerance),
            "tolerance_ppm": False,
        },
        height=350,
    )
    sequence_view(key="seq_view", state_manager=state_manager)

    # LinePlot linked to SequenceView
    st.subheader("Fragment Spectrum")
    lineplot = _create_component(
        LinePlot, "spectrum_plot", CACHE_PATH, needs_regeneration,
        data=peaks_lf, filters={"selected_psm": "psm_id"},
        filter_defaults={"selected_psm": 0}, regenerate_cache=True,
        x_column="x", y_column="y", highlight_column="highlight",
        annotation_column="annotation", x_label="m/z", y_label="Intensity",
    )
    lineplot(key="spectrum_plot", state_manager=state_manager, sequence_view_key="seq_view")

    # Heatmap
    st.subheader("RT x m/z Heatmap")
    heatmap_comp = _create_component(
        Heatmap, "rt_heatmap", CACHE_PATH, needs_regeneration,
        data=heatmap_lf, x_column="x", y_column="y", intensity_column="intensity",
        x_label="Retention Time (s)", y_label="m/z", log_scale=True, regenerate_cache=True,
    )
    heatmap_comp(key="rt_heatmap", state_manager=state_manager)

    if needs_regeneration:
        mark_rendered()
