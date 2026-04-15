"""
Workspace management — switcher, create/delete, safe-name validation.

The switcher is the only Streamlit widget in this module; everything else
is pure filesystem logic so it can be unit-tested without a Streamlit app.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import streamlit as st

from streamlit_app.utils.state.core import PROJECT_ROOT, STATE_PREFIX


def is_safe_workspace_name(name: str) -> bool:
    """Check if a workspace name is safe (no path traversal characters)."""
    if not name:
        return False
    return "/" not in name and "\\" not in name and name not in ("..", ".")


def _get_workspaces_dir() -> Path:
    """Return the workspaces directory from settings.

    Absolute paths in ``workspaces_dir`` are honored as-is. Relative paths
    are anchored to the project root (not the CWD) so workspaces always
    land at ``<project-root>/workspaces/`` regardless of where
    ``streamlit run`` is invoked from.
    """
    configured = Path(st.session_state.settings.get("workspaces_dir", "workspaces"))
    if configured.is_absolute():
        return configured
    return PROJECT_ROOT / configured


def _render_workspace_switcher(params: dict) -> None:
    """Render workspace selector, create, and delete controls."""
    # Local import avoids a circular reference (params -> workspace -> params).
    from streamlit_app.utils.state.params import save_params

    workspaces_dir = _get_workspaces_dir()

    with st.expander("**Workspaces**"):
        def change_workspace():
            # Flush any in-flight session state to the OUTGOING workspace's
            # params.json BEFORE deleting keys. Without this, peptides the
            # user typed in a text_area whose blur event races with the
            # dropdown on_change can be lost on switch.
            save_params(params)

            for key in params:
                st.session_state.pop(key, None)
            # Clear page-namespaced keys so stale analysis state from
            # the previous workspace is not shown in the new one.
            stale_keys = [
                k for k in st.session_state
                if k.startswith(STATE_PREFIX)
            ]
            for key in stale_keys:
                st.session_state.pop(key, None)
            st.query_params.workspace = st.session_state["chosen-workspace"]

        workspaces_dir.mkdir(parents=True, exist_ok=True)
        options = sorted(
            f.name for f in workspaces_dir.iterdir() if f.is_dir()
        )
        if not options:
            options = ["default"]

        current_workspace = str(st.session_state.workspace.stem)
        index = options.index(current_workspace) if current_workspace in options else 0

        st.selectbox(
            "choose existing workspace",
            options, index=index,
            on_change=change_workspace, key="chosen-workspace",
        )

        new_ws_name = st.text_input("New workspace name", "", key="new-workspace-name").strip()
        if st.button("**Create Workspace**"):
            if new_ws_name and is_safe_workspace_name(new_ws_name):
                new_path = Path(workspaces_dir, new_ws_name)
                new_path.mkdir(parents=True, exist_ok=True)
                st.session_state.workspace = new_path
                st.query_params.workspace = new_ws_name
                st.rerun()
            else:
                st.warning("Please enter a valid workspace name.")

        if current_workspace != "default":
            confirm_key = "_confirm_delete_workspace"
            if st.session_state.get(confirm_key):
                st.warning(f"Delete workspace **{current_workspace}** and all its data?")
                yes_col, no_col = st.columns(2)
                with yes_col:
                    if st.button("Yes, delete", type="primary", key="ws-confirm-yes"):
                        ws_path = Path(workspaces_dir, current_workspace)
                        if ws_path.exists():
                            shutil.rmtree(ws_path)
                        st.session_state.workspace = Path(workspaces_dir, "default")
                        st.query_params.workspace = "default"
                        st.session_state.pop(confirm_key, None)
                        st.rerun()
                with no_col:
                    if st.button("Cancel", key="ws-confirm-no"):
                        st.session_state.pop(confirm_key, None)
                        st.rerun()
            else:
                if st.button("Delete this workspace"):
                    st.session_state[confirm_key] = True
                    st.rerun()
