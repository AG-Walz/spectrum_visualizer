"""Export action layer — called by export_page.py to wire the Export button to core.export.run_export.
Collects session state into the shared dataclasses and routes to the CLI-shared orchestrator.
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import streamlit as st

from streamlit_app.utils.state import get_param, get_state_key
from streamlit_app.utils.page_state import build_run_settings
from streamlit_app.utils.progress_tracker import progress_status
from streamlit_app.views.export.gallery import _render_plot_gallery
from core.defaults import DEFAULTS
from core.export import (
    ExportCallbacks,
    ExportSelection,
    PlotSettings,
    run_export,
)

logger = logging.getLogger(__name__)


def _collect_session_frames() -> dict:
    """Pull the four Polars LazyFrames (psm/peaks/sequence/heatmap) out of session state.
    """
    return {
        "psm":      st.session_state.get(get_state_key("psm_frame")),
        "peaks":    st.session_state.get(get_state_key("peaks_frame")),
        "sequence": st.session_state.get(get_state_key("sequence_frame")),
        "heatmap":  st.session_state.get(get_state_key("heatmap_frame")),
    }


def _build_plot_settings(params: dict) -> PlotSettings:
    """Resolve PlotSettings (formats, figure size, workers) from session params + defaults.
    """
    return PlotSettings(
        export_formats=["png", "svg"],
        figure_size=DEFAULTS.figure_size,
        export_resolution=DEFAULTS.export_resolution,
        plot_unpaired_queries=get_param(
            params, "plot_unpaired_queries", DEFAULTS.plot_unpaired_queries,
        ),
        max_peak_labels=get_param(params, "max_peak_labels", DEFAULTS.max_peak_labels),
        workers=get_param(params, "workers", DEFAULTS.workers),
    )


def _show_export_summary(export_result, tracker) -> None:
    """Stamp mirror/summary plot counts into the status block, fill the bar, and mark it complete.
    """
    if export_result.mirror_plot_count:
        st.caption(f"  Generated {export_result.mirror_plot_count} mirror plot file(s)")
    if export_result.summary_plot_count:
        st.caption(f"  Generated {export_result.summary_plot_count} summary plot(s)")
    tracker.finalize(
        f"Exported {len(export_result.saved_files)} files to {export_result.export_root}"
    )


def _show_saved_files_expander(export_result) -> None:
    """Render a collapsed expander listing every file written by the export.
    """
    if not export_result.saved_files:
        return
    with st.expander(f"Exported files ({len(export_result.saved_files)})", expanded=False):
        for file_path in export_result.saved_files:
            st.caption(f"  {file_path}")


def _run_export(params: dict, export_dir: str, export_parquet: bool, export_csv: bool,
                export_mirror: bool, export_summary: bool, export_html: bool) -> Path | None:
    """Drive core.export.run_export from the Export button and render the resulting gallery.
    """
    # Guard: refuse to export when Load & Annotate hasn't populated session state
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    if not annotated_queries:
        st.warning("No annotated spectra in session -- re-run **Load & Annotate** first.")
        return None

    # Build config: format selection + plot settings shared with the CLI orchestrator
    selection = ExportSelection(
        parquet=export_parquet, csv=export_csv,
        mirror_plots=export_mirror, summary_plots=export_summary, html=export_html,
    )
    plot_settings = _build_plot_settings(params)

    # Dispatch: run shared export pipeline inside a live status block with a plot progress bar
    with progress_status("Exporting...") as tracker:
        export_result = run_export(
            export_directory=export_dir,
            project_name=DEFAULTS.project_name,
            run_settings=build_run_settings(),
            frames=_collect_session_frames(),
            annotated_queries=annotated_queries,
            reference_store=st.session_state.get(get_state_key("reference_store"), {}),
            peptidoform_psms=st.session_state.get(get_state_key("peptidoform_psms"), {}),
            comparison_frame=st.session_state.get(get_state_key("comparison_frame")),
            selection=selection,
            plot_settings=plot_settings,
            callbacks=ExportCallbacks(
                on_step=st.caption,
                on_plot_progress=lambda done, total: tracker.update(done, total, units="plots"),
                on_warning=st.warning,
            ),
        )
        _show_export_summary(export_result, tracker)

    # Post-export: list saved files and open the plot gallery below
    _show_saved_files_expander(export_result)
    _render_plot_gallery(export_result.export_root)
    return export_result.export_root
