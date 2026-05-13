"""
Sidebar, workspace switcher, and page-setup entry point.

``page_setup()`` is called at the top of every view; ``render_sidebar``
renders the shared left-hand controls; ``_render_workspace_switcher`` is the
workspace selector within the sidebar.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import streamlit as st

from core.common.defaults import DEFAULTS

from streamlit_app.utils.state.state_core import PROJECT_ROOT, STATE_PREFIX, _persisted_key
from streamlit_app.utils.state.state_params import _PARAM_KEYS, load_params, save_params
from streamlit_app.utils.state.state_workspace import (
    _get_workspaces_dir,
    is_safe_workspace_name,
)


def _flush_and_switch_workspace(params: dict) -> None:
    """Selectbox callback: persist outgoing params, clear stale state, switch.

    Keeping this at module level (instead of as a nested closure) makes the
    workspace-switcher easy to read: each widget just dispatches to a named
    helper and the heavy state surgery lives here.
    """
    # Flush: write OUTGOING workspace's params.json before mutating state
    # Flush any in-flight session state to the OUTGOING workspace's
    # params.json BEFORE deleting keys. Without this, peptides the user
    # typed in a text_area whose blur event races with the dropdown
    # on_change can be lost on switch.
    save_params(params)

    # Reset widgets: pop every canonical widget key + anything from params
    # Pop every known widget key — not just the ones that happened to
    # be on disk in ``params.json`` — so settings from the outgoing
    # workspace do not bleed into the incoming one's widgets (and from
    # there, back into its params.json via save_params).
    for key in (*_PARAM_KEYS, *params.keys()):
        st.session_state.pop(key, None)

    # Reset analysis: drop page-namespaced keys + _persisted backups
    # Clear page-namespaced keys so stale analysis state (and the
    # ``_persisted`` widget-value backups) from the previous workspace
    # are not shown in the new one.
    stale_keys = [key for key in st.session_state if key.startswith(STATE_PREFIX)]
    for key in stale_keys:
        st.session_state.pop(key, None)

    # Switch & rerun: update query param and halt the in-flight rerun
    st.query_params.workspace = st.session_state["chosen-workspace"]
    # Halt this rerun immediately to avoid the script continuing with
    # session_state.workspace still pointing at the OLD workspace.
    st.rerun()


def _render_workspace_select(workspaces_dir: Path, params: dict) -> str:
    """Render the existing-workspace dropdown; returns the current workspace name."""
    workspaces_dir.mkdir(parents=True, exist_ok=True)
    options = sorted(entry.name for entry in workspaces_dir.iterdir() if entry.is_dir())
    if not options:
        options = ["default"]

    current_workspace = str(st.session_state.workspace.stem)
    index = options.index(current_workspace) if current_workspace in options else 0
    st.selectbox(
        "choose existing workspace",
        options, index=index,
        on_change=_flush_and_switch_workspace, args=(params,),
        key="chosen-workspace",
    )
    return current_workspace


def _render_create_workspace(workspaces_dir: Path) -> None:
    """Render the 'New workspace name' input and Create button."""
    new_workspace_name = st.text_input("New workspace name", "", key="new-workspace-name").strip()
    if not st.button("**Create Workspace**"):
        return
    if new_workspace_name and is_safe_workspace_name(new_workspace_name):
        new_path = Path(workspaces_dir, new_workspace_name)
        new_path.mkdir(parents=True, exist_ok=True)
        st.session_state.workspace = new_path
        st.query_params.workspace = new_workspace_name
        st.rerun()
    else:
        st.warning("Please enter a valid workspace name.")


def _render_delete_workspace(workspaces_dir: Path, current_workspace: str) -> None:
    """Render the delete-with-confirmation flow for the current workspace."""
    if current_workspace == "default":
        return
    confirm_key = "_confirm_delete_workspace"
    if not st.session_state.get(confirm_key):
        if st.button("Delete this workspace"):
            st.session_state[confirm_key] = True
            st.rerun()
        return

    st.warning(f"Delete workspace **{current_workspace}** and all its data?")
    confirm_column, cancel_column = st.columns(2)
    with confirm_column:
        if st.button("Yes, delete", type="primary", key="workspace-delete-confirm-yes"):
            workspace_path = Path(workspaces_dir, current_workspace)
            if workspace_path.exists():
                shutil.rmtree(workspace_path)
            st.session_state.workspace = Path(workspaces_dir, "default")
            st.query_params.workspace = "default"
            st.session_state.pop(confirm_key, None)
            st.rerun()
    with cancel_column:
        if st.button("Cancel", key="workspace-delete-confirm-no"):
            st.session_state.pop(confirm_key, None)
            st.rerun()


def _render_workspace_switcher(params: dict) -> None:
    """Render workspace selector, create, and delete controls."""
    workspaces_dir = _get_workspaces_dir()
    with st.expander("**Workspaces**"):
        current_workspace = _render_workspace_select(workspaces_dir, params)
        _render_create_workspace(workspaces_dir)
        _render_delete_workspace(workspaces_dir, current_workspace)


def page_setup(page: str = "") -> dict[str, Any]:
    """Set up page config, workspace, and sidebar. Call at the top of every page."""
    # Settings: lazy-load settings.json into session_state on first call
    if "settings" not in st.session_state:
        settings_path = PROJECT_ROOT / "streamlit_app" / "config" / "settings.json"
        with open(settings_path, "r") as handle:
            st.session_state.settings = json.load(handle)

    # Page config: title, favicon, wide layout, sidebar default
    logo_path = PROJECT_ROOT / "assets" / "openms_transparent_bg_logo.svg"
    st.set_page_config(
        page_title=st.session_state.settings.get("app_name", "MSpecViz"),
        page_icon=str(logo_path) if logo_path.exists() else None,
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    # Workspace resolution: switch from query param or fall back to default
    if "workspace" not in st.session_state or (
        "workspace" in st.query_params
        and st.query_params.workspace != st.session_state.workspace.name
    ):
        st.cache_data.clear()
        st.cache_resource.clear()
        st.session_state.location = "local"

        workspaces_dir = _get_workspaces_dir()
        if st.session_state.settings.get("enable_workspaces", False):
            if "workspace" in st.query_params:
                requested_workspace = st.query_params.workspace
                if not is_safe_workspace_name(requested_workspace):
                    st.session_state.workspace = Path(workspaces_dir, "default")
                    st.query_params.workspace = "default"
                else:
                    st.session_state.workspace = Path(workspaces_dir, requested_workspace)
            else:
                st.session_state.workspace = Path(workspaces_dir, "default")
                st.query_params.workspace = "default"
        else:
            st.session_state.workspace = Path(workspaces_dir, "default")

    # Finalize: sync query param, ensure workspace dir exists, render sidebar
    if (
        "workspace" not in st.query_params
        and st.session_state.settings.get("enable_workspaces", False)
    ):
        st.query_params.workspace = st.session_state.workspace.name

    st.session_state.workspace.mkdir(parents=True, exist_ok=True)

    return render_sidebar(page)


def render_sidebar(page: str = "") -> dict[str, Any]:
    """Render sidebar with workspace switcher, settings, and version info."""
    # Load params for the active workspace
    params = load_params()

    with st.sidebar:
        # Workspace switcher: only when workspaces are enabled
        if st.session_state.settings.get("enable_workspaces", False):
            _render_workspace_switcher(params)

        # Settings expander: image export format selector, persisted across pages
        with st.expander("**Settings**"):
            image_formats = ["svg", "png", "jpeg", "webp"]
            persisted = _persisted_key("image_format")
            # Seed session_state BEFORE the widget is created so Streamlit
            # uses the restored value on page navigation. Reading from the
            # non-widget persisted backup survives widget garbage collection.
            if "image_format" not in st.session_state:
                restored = st.session_state.get(
                    persisted, params.get("image_format", DEFAULTS.image_format)
                )
                if restored not in image_formats:
                    restored = DEFAULTS.image_format
                st.session_state["image_format"] = restored
            selected_format = st.selectbox(
                "image export format",
                image_formats,
                key="image_format",
            )
            st.session_state[persisted] = selected_format
            params["image_format"] = selected_format
            save_params(params)

        # Version footer: branded app name + version chip
        version_info = st.session_state.settings.get("version", "0.1.0")
        app_name = st.session_state.settings.get("app_name", "MSpecViz")
        st.markdown(
            '<style>.version-box{border:1px solid #a4a5ad;padding:10px;'
            'border-radius:.5rem;text-align:center;display:flex;'
            'justify-content:center;align-items:center}</style>',
            unsafe_allow_html=True,
        )
        st.markdown(
            f'<div class="version-box">{app_name}<br>Version: {version_info}</div>',
            unsafe_allow_html=True,
        )

    return params
