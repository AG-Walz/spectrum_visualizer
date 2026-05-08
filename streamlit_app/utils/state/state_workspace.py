"""
Workspace management — safe-name validation and directory resolution.

Pure filesystem logic so this module can be imported without a running
Streamlit app. The workspace switcher widget lives in ``params.py``.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from streamlit_app.utils.state.state_core import PROJECT_ROOT


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
