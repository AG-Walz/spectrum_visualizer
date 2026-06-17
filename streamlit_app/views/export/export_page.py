"""Rendered as the Export tab — format checkboxes, run button, history browser, and ZIP download.
Sections are @st.fragments; heavy work lives in export_actions.py and gallery.py.
"""
import io
import json
import logging
import shutil
import zipfile
from pathlib import Path

import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_param, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import apply_run_settings, restore_analysis_if_needed

from streamlit_app.views.export.gallery import _render_plot_gallery
from streamlit_app.views.export.export_actions import _run_export

params = page_setup()
logger = logging.getLogger(__name__)


def _list_export_runs(export_path: Path) -> list[tuple[Path, dict]]:
    """Scan export_path for past runs; return [(run_dir, run_settings_dict)] newest-first.
    """
    runs = []
    for settings_file in sorted(export_path.glob("*/run_settings.json"), reverse=True):
        try:
            with open(settings_file) as f:
                run_settings = json.load(f)
            runs.append((settings_file.parent, run_settings))
        except Exception:
            continue
    return runs


def _format_run_label(run_dir: Path, run_settings: dict) -> str:
    """Build the one-line expander label (name, date, peptide/file counts) for one historical run.
    """
    peptides = run_settings.get("peptides", [])
    idxml_count = len(run_settings.get("idxml_files", []))
    mzml_count = len(run_settings.get("mzml_files", []))
    created = (run_settings.get("_created_at") or "")[:10]
    return (
        f"**{run_dir.name}**"
        f" -- {created}"
        f" -- {len(peptides)} peptide(s)"
        f" -- {idxml_count} idXML / {mzml_count} mzML"
    )


def _render_run_metadata(run_settings: dict) -> None:
    """Render the metadata column (peptides, files, tolerance, timestamps) inside a history expander.
    """
    peptides = run_settings.get("peptides", [])
    peptide_preview = ", ".join(peptides[:3])
    if len(peptides) > 3:
        peptide_preview += f" ... +{len(peptides) - 3}"
    idxml_count = len(run_settings.get("idxml_files", []))
    mzml_count = len(run_settings.get("mzml_files", []))

    st.markdown(f"**Peptides:** {peptide_preview or '--'}")
    st.markdown(f"**Files:** {idxml_count} idXML, {mzml_count} mzML")
    st.markdown(f"**Reference substring:** `{run_settings.get('reference_substring') or '--'}`")
    st.markdown(f"**Ion types:** {', '.join(run_settings.get('ion_types', []))}")
    st.markdown(f"**m/z tolerance:** {run_settings.get('mz_tolerance', '--')} Da")
    created_full = (run_settings.get("_created_at") or "--")[:19]
    last_run_full = (run_settings.get("_last_run") or "--")[:19]
    st.caption(f"Created: {created_full}  |  Last run: {last_run_full}")


def _render_run_action_buttons(run_dir: Path, run_settings: dict) -> None:
    """Render Reload/Delete buttons for one historical run, with a two-step confirm before delete.
    """
    # Reload: rehydrate this run's parameters into session state
    if st.button("Reload settings", key=f"reload_{run_dir.name}"):
        apply_run_settings(run_settings)
        st.success("Settings loaded -- switch to Setup page to review.")
        st.rerun()

    # Two-step delete: confirmation flag is set by the initial click, then this branch fires
    confirm_key = f"confirm_delete_{run_dir.name}"
    if st.session_state.get(confirm_key):
        st.warning("Delete all files in this run?")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Yes, delete", key=f"yes_del_{run_dir.name}", type="primary"):
                shutil.rmtree(run_dir)
                st.session_state.pop(confirm_key, None)
                st.rerun()
        with c2:
            if st.button("Cancel", key=f"cancel_del_{run_dir.name}"):
                st.session_state.pop(confirm_key, None)
                st.rerun()
    # Initial click: arm the confirmation flag and rerun so the branch above shows next
    else:
        if st.button("Delete run", key=f"delete_{run_dir.name}"):
            st.session_state[confirm_key] = True
            st.rerun()


def _render_history_run(run_dir: Path, run_settings: dict) -> None:
    """Render the full expander (metadata, actions, embedded gallery) for one historical run.
    """
    with st.expander(_format_run_label(run_dir, run_settings), expanded=False):
        col_meta, col_actions = st.columns([3, 1])
        with col_meta:
            _render_run_metadata(run_settings)
        with col_actions:
            _render_run_action_buttons(run_dir, run_settings)
        _render_plot_gallery(run_dir)


@st.fragment
def export_history_fragment(export_dir: str):
    """Browse past export runs in export_dir, with reload-settings and delete actions per run.
    """
    st.subheader("Export History")

    # Empty states: directory missing or no completed runs inside
    export_path = Path(export_dir)
    if not export_path.exists():
        st.caption(f"No exports directory found at `{export_dir}`.")
        return

    runs = _list_export_runs(export_path)
    if not runs:
        st.caption(f"No completed export runs found in `{export_dir}`.")
        return

    # Render: one expander per run, newest first
    st.caption(f"{len(runs)} run(s) in `{export_dir}`")
    for run_dir, run_settings in runs:
        _render_history_run(run_dir, run_settings)


@st.fragment
def zip_download_fragment(export_root: Path) -> None:
    """Render a download button that zips export_root on click 
    (deferred so the zip isn't built on every rerun).
    """
    def _make_zip() -> bytes:
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for file in export_root.rglob("*"):
                if not file.is_file():
                    continue
                zip_file.write(file, file.relative_to(export_root.parent))
        return buffer.getvalue()

    st.download_button(
        label="Download Export as ZIP",
        data=_make_zip,
        file_name=f"{export_root.name}.zip",
        mime="application/zip",
    )


@st.fragment
def export_controls_fragment():
    """Render the export-directory input, format checkboxes, and the "Export All" trigger button.
    """
    # Guard: Load & Annotate must have populated session state first
    if not st.session_state.get(get_state_key("data_ready")):
        st.info("Load data first (Load & Annotate page).")
        return

    # Export directory: default to <workspace>/exports when unset
    workspace = st.session_state.get("workspace", None)
    export_fallback = str(Path(workspace) / "exports") if workspace else "./exports"
    export_default = params.get("export_directory", "") or export_fallback
    st.text_input(
        "Export Directory",
        value=export_default,
        key="export_directory",
    )

    # Format checkboxes: dataframes/summary in column 1, plots in column 2
    export_col1, export_col2 = st.columns(2)
    with export_col1:
        export_parquet = st.checkbox("Parquet DataFrames", value=True,
                                     key=f"{STATE_PREFIX}_exp_parquet")
        export_csv = st.checkbox("Comparison Summary (CSV/XLSX)", value=True,
                                 key=f"{STATE_PREFIX}_exp_csv")
        export_html = st.checkbox("Interactive Plots (HTML)", value=False,
                                  key=f"{STATE_PREFIX}_exp_html")
    with export_col2:
        export_mirror = st.checkbox("Mirror Plots (PNG/SVG)", value=True,
                                    key=f"{STATE_PREFIX}_exp_mirror")
        export_summary = st.checkbox("Summary Analysis Plots", value=True,
                                     key=f"{STATE_PREFIX}_exp_summary")

    # Trigger: dispatch _run_export with the chosen flags
    if st.button("Export All", key=f"{STATE_PREFIX}_export_btn", type="primary"):
        _run_export(
            params,
            get_param(params, "export_directory", "") or export_fallback,
            export_parquet, export_csv, export_mirror,
            export_summary, export_html,
        )


def render():
    """Page entry point: rebuild analysis cache, render the three fragments, then persist params.
    """
    # Restore cached analysis (if any) and render the page header
    restore_analysis_if_needed()
    st.title("Export")
    st.caption("Export analysis results, plots, and data files")

    # Section 1: format selection + Export button
    export_controls_fragment()

    # Resolve the export directory once, then reuse it for history + ZIP
    workspace = st.session_state.get("workspace", None)
    export_fallback = str(Path(workspace) / "exports") if workspace else "./exports"
    export_dir = get_param(params, "export_directory", "") or export_fallback

    # Section 2: history browser over previous runs
    st.markdown("---")
    export_history_fragment(export_dir)

    # Section 3: ZIP download (only when the directory has content)
    export_path = Path(export_dir)
    if export_path.exists() and any(export_path.rglob("*")):
        st.markdown("---")
        st.subheader("Download")
        zip_download_fragment(export_path)

    save_params(params)


render()
