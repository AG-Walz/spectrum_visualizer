"""Session-state, workspace, params, and sidebar helpers for the Streamlit app.

Combines four layers in one file:

- Namespaced `st.session_state` key helpers (`STATE_PREFIX`, `get_state_key`,
  `_persisted_key`, `get_param`, `PROJECT_ROOT`).
- Workspace directory helpers (`is_safe_workspace_name`, `_get_workspaces_dir`).
- Params load/save against each workspace's `params.json`
  (`load_params`, `save_params`).
- The sidebar + `page_setup` entry point called at the top of every view.
"""
from __future__ import annotations

import functools
import json
import logging
import os
import shutil
import tempfile
import traceback
from collections import deque
from pathlib import Path
from typing import Any

import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

try:
    import psutil
except ImportError:
    psutil = None

from core import __version__
from core.bug_report import build_issue_url, build_snapshot
from core.defaults import DEFAULTS


# ============================================================================
# Namespaced session_state key helpers
# ============================================================================

# Prefix on every app-owned session_state key, so they cannot collide with third-party widget keys.
STATE_PREFIX = "peptide_mirror"

# Anchor for CWD-relative paths (workspaces, assets, settings) so the app behaves the
# same no matter which directory streamlit run was launched from.
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def get_state_key(key: str) -> str:
    """Prefix key with the app namespace (e.g. "foo" -> "peptide_mirror.foo")."""
    return f"{STATE_PREFIX}.{key}"


def _persisted_key(key: str) -> str:
    """Return the <key>_persisted backup-key name — survives Streamlit's widget garbage collection on page navigation."""
    return get_state_key(key + "_persisted")


def get_param(params: dict, key: str, default: Any = None) -> Any:
    """Resolve a param with fallback order: live widget > _persisted backup > params.json > default."""
    if key in st.session_state:
        return st.session_state[key]
    persisted = _persisted_key(key)
    if persisted in st.session_state:
        return st.session_state[persisted]
    return params.get(key, default)


# ============================================================================
# Loaded-idXML PSMList cache
# ============================================================================
# The Setup picker (peptide counts) and Load & Annotate (target lookup) both
# consume a PSMList loaded from the selected idXML files. Caching that PSMList
# keyed on the file selection lets the second stage reuse the first stage's load
# instead of re-parsing the same files (see core.peptide_lookup.load_idxmls).

_LOADED_IDXML_PSMLIST_KEY = "_loaded_idxml_psmlist"
_LOADED_IDXML_FILES_KEY = "_loaded_idxml_files"


def cached_idxml_psmlist(idxml_files: list[str]):
    """Return the PSMList previously loaded from exactly idxml_files, or None when
    the cache is empty or was loaded from a different idXML selection.
    """
    if st.session_state.get(_LOADED_IDXML_FILES_KEY) == idxml_files:
        return st.session_state.get(_LOADED_IDXML_PSMLIST_KEY)
    return None


def cache_idxml_psmlist(idxml_files: list[str], psm_list) -> None:
    """Stash a freshly loaded PSMList keyed on the idXML file selection it came from."""
    st.session_state[_LOADED_IDXML_PSMLIST_KEY] = psm_list
    st.session_state[_LOADED_IDXML_FILES_KEY] = list(idxml_files)


# ============================================================================
# Workspace directory helpers
# ============================================================================

def is_safe_workspace_name(name: str) -> bool:
    """Reject empty names and any name containing path-traversal characters (/, \\, ., ..)."""
    if not name:
        return False
    return "/" not in name and "\\" not in name and name not in ("..", ".")


def _get_workspaces_dir() -> Path:
    """Return the workspaces root from settings.workspaces_dir — relative paths are anchored to the project root, not CWD."""
    configured = Path(st.session_state.settings.get("workspaces_dir", "workspaces"))
    if configured.is_absolute():
        return configured
    return PROJECT_ROOT / configured


# ============================================================================
# Params load/save (per-workspace params.json)
# ============================================================================

# Canonical list of widget keys save_params writes to params.json. Listing every key
# explicitly (not just those already on disk) lets first-ever saves capture widget values, and
# lets _flush_and_switch_workspace pop them when switching workspaces.
_PARAM_KEYS: tuple[str, ...] = (
    "selected_idxmls", "selected_mzmls",
    "idxml_browser_dir", "mzml_browser_dir",
    "target_peptides", "reference_substring",
    "ion_types", "mz_tolerance",
    "add_neutral_losses",
    "intensity_mode", "workers",
    "koina_enabled", "koina_model",
    "koina_collision_energy", "koina_instrument_type",
    "image_format",
    "summary_metrics", "summary_plot_types", "summary_sort",
)


def load_params(default: bool = False) -> dict[str, Any]:
    """Load the current workspace's params.json — empty dict if workspaces disabled or no file yet."""
    # Short-circuit: workspaces disabled or caller asked for defaults
    if default or not st.session_state.settings.get("enable_workspaces", False):
        return {}

    # Read: load params.json from the current workspace
    path = Path(st.session_state.workspace, "params.json")
    if not path.exists():
        return {}

    with open(path, "r", encoding="utf-8") as handle:
        raw_params = json.load(handle)

    # Migrate legacy hyphenated keys (reference-substring) to snake_case. Pre-migration
    # workspaces would otherwise hold stale hyphen keys that save_params cannot sync.
    migrated: dict[str, Any] = {key: value for key, value in raw_params.items() if "-" not in key}
    for key, value in raw_params.items():
        if "-" in key:
            migrated.setdefault(key.replace("-", "_"), value)
    return migrated


def save_params(params: dict[str, Any]) -> None:
    """Pull live widget values (or _persisted backups) into params and write atomically to params.json."""
    # Short-circuit: workspaces disabled means nothing to persist
    if not st.session_state.settings.get("enable_workspaces", False):
        return

    # Sync: pull live widget values (or _persisted backups) into the params dict.
    # Union with _PARAM_KEYS so first-ever saves capture widget values not yet on disk.
    tracked_keys = list(dict.fromkeys((*_PARAM_KEYS, *params.keys())))
    for key in tracked_keys:
        if key in st.session_state:
            params[key] = st.session_state[key]
        else:
            # _persisted backup (set by setup.py fragments) — survives widget garbage collection on page navigation.
            persisted = _persisted_key(key)
            if persisted in st.session_state:
                params[key] = st.session_state[persisted]

    # Atomic write via unique staging file. A single rerun can call save_params from the sidebar
    # plus several setup fragments; sharing one params.json.tmp path between them — or between
    # browser tabs on the same workspace — made the second mover's os.replace fail.
    params_path = Path(st.session_state.workspace, "params.json")
    params_path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, staging_name = tempfile.mkstemp(
        prefix=".params-", suffix=".json.tmp", dir=params_path.parent,
    )
    staging_path = Path(staging_name)
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as outfile:
            json.dump(params, outfile, indent=4)
        os.replace(staging_path, params_path)
    except Exception:
        if staging_path.exists():
            try:
                staging_path.unlink()
            except OSError:
                pass
        raise


# ============================================================================
# Sidebar + page_setup entry point
# ============================================================================

def _flush_and_switch_workspace(params: dict) -> None:
    """Workspace-dropdown on_change: save the outgoing workspace, wipe stale state, then rerun on the new one."""
    # Flush: write the OUTGOING workspace's params.json before mutating state. Without this,
    # peptides typed in a text_area whose blur event races with the dropdown on_change are lost.
    save_params(params)

    # Reset widgets: pop every canonical widget key (not just what was on disk) so settings from
    # the outgoing workspace don't bleed into the incoming one — and from there, back into its
    # params.json via save_params.
    for key in (*_PARAM_KEYS, *params.keys()):
        st.session_state.pop(key, None)

    # Reset analysis: clear page-namespaced keys and _persisted widget backups so stale
    # analysis state from the previous workspace doesn't show up in the new one.
    stale_keys = [key for key in st.session_state if key.startswith(STATE_PREFIX)]
    for key in stale_keys:
        st.session_state.pop(key, None)

    # Switch & rerun. Halt this rerun immediately so the script doesn't continue with
    # session_state.workspace still pointing at the OLD workspace.
    st.query_params.workspace = st.session_state["chosen-workspace"]
    st.rerun()


def _render_workspace_select(workspaces_dir: Path, params: dict) -> str:
    """Render the existing-workspace dropdown and return the current workspace name."""
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
    """Render the "New workspace name" input and the Create button."""
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
    """Render the two-step delete-with-confirmation flow for the current workspace."""
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
    """Render the Workspaces expander: selector, create, and delete controls."""
    workspaces_dir = _get_workspaces_dir()
    with st.expander("**Workspaces**"):
        current_workspace = _render_workspace_select(workspaces_dir, params)
        _render_create_workspace(workspaces_dir)
        _render_delete_workspace(workspaces_dir, current_workspace)


@st.fragment(run_every=5)
def _render_resource_monitor() -> None:
    """Point-in-time system CPU/RAM pinned into the top-right header strip; self-refreshes every 5s.
    Silent when psutil is absent (optional dependency)."""
    if psutil is None:
        return
    memory = psutil.virtual_memory()
    st.markdown(
        f"""
        <div class="mspecviz-resource-monitor">🖥️ CPU {psutil.cpu_percent(interval=0.1):.0f}% &middot; 🧠 RAM {memory.percent:.0f}% ({memory.used / 1e9:.1f}/{memory.total / 1e9:.1f} GB)</div>
        <style>
        .mspecviz-resource-monitor {{
            position: fixed;
            /* Sit just below Streamlit's 3.75rem header instead of inside it, so the chip
               never shares the toolbar row with the "Running" status widget (whose width
               varies with state/locale) — see issue #103. */
            top: 3.9rem;
            right: 1rem;
            z-index: 1000000;
            font-size: 0.8rem;
            white-space: nowrap;
            color: rgba(128, 128, 128, 0.95);
            pointer-events: none;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def page_setup(page: str = "") -> dict[str, Any]:
    """Initialise settings, page config, workspace, and sidebar — call once at the top of every view."""
    # Terminal-error capture: attach the session log collector once per process so
    # the bug report can surface errors logged anywhere this session.
    _install_terminal_log_collector()

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

    # Workspace resolution: honour ?workspace=… query param, else fall back to "default"
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

    _render_resource_monitor()
    return render_sidebar(page)


def capture_last_error(func):
    """Decorator: stash an exception's traceback for the bug report, then re-raise
    (Streamlit's RerunException is a BaseException, so it slips through)."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            st.session_state[get_state_key("last_error")] = traceback.format_exc()
            raise
    return wrapper


# Session terminal-error log: a root-logger handler stashes ERROR+ records into the
# active session's ring buffer, so the bug report surfaces every error logged to the
# terminal this session — not just the last decorated-action exception.
TERMINAL_LOG_KEY = "terminal_log"
_TERMINAL_LOG_MAXLEN = 10        # keep the most recent errors; tracebacks are bulky
_TERMINAL_LOG_RECORD_CHARS = 2000  # cap each record so one trace can't dominate the URL


class _SessionLogCollector(logging.Handler):
    """Route ERROR+ log records into the active Streamlit session's ring buffer.
    Records logged outside a script run (startup, ProcessPool workers) have no
    session context attached, so they are dropped rather than misfiled."""
    def emit(self, record: logging.LogRecord) -> None:
        if get_script_run_ctx() is None:
            return  # no active session on this thread — nothing to attach to
        try:
            message = self.format(record)[:_TERMINAL_LOG_RECORD_CHARS]
            key = get_state_key(TERMINAL_LOG_KEY)
            buffer = st.session_state.get(key)
            if buffer is None:
                buffer = deque(maxlen=_TERMINAL_LOG_MAXLEN)
                st.session_state[key] = buffer
            buffer.append(message)
        except Exception:
            pass  # logging the error must never itself raise


def _install_terminal_log_collector() -> None:
    """Attach the session log collector to the root logger exactly once per process."""
    root = logging.getLogger()
    if any(isinstance(handler, _SessionLogCollector) for handler in root.handlers):
        return
    handler = _SessionLogCollector()
    handler.setLevel(logging.ERROR)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    root.addHandler(handler)


def _render_bug_report(page: str) -> None:
    """Render the 'Report a bug / request a feature' expander, pre-filling a GitHub
    issue with a tiny snapshot: workspace (a pointer to the server-side params.json),
    page, the last captured error, and any errors logged to the terminal this session."""
    repo = st.session_state.settings.get("bug_report_repo", "AG-Walz/spectrum_visualizer")
    workspace = st.session_state.workspace.name
    error = st.session_state.get(get_state_key("last_error"))
    terminal_log = list(st.session_state.get(get_state_key(TERMINAL_LOG_KEY), []))

    bug_url, inlined = build_issue_url(
        repo, "bug_report.yml",
        build_snapshot(workspace, page, error, terminal_log),
    )
    feature_url, _ = build_issue_url(
        repo, "feature_request.yml",
        build_snapshot(workspace, page),
    )

    with st.expander("🐞 **Report a bug / request a feature**"):
        st.caption("Opens a pre-filled GitHub issue — use a desktop browser (the mobile app drops pre-filled fields).")
        st.link_button("🐞 Report a bug", bug_url, use_container_width=True)
        st.link_button("💡 Request a feature", feature_url, use_container_width=True)
        if not inlined:
            st.warning("Some captured context was too large to embed — please paste any terminal errors into the issue.")


def render_sidebar(page: str = "") -> dict[str, Any]:
    """Render the workspace switcher, Settings expander, and version footer; return the active workspace's params."""
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
            # Seed session_state from the _persisted backup BEFORE creating the widget — otherwise
            # Streamlit's default ("svg") wins over the user's saved choice after page navigation.
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

        # Bug report: inline a compact snapshot into a pre-filled GitHub issue
        _render_bug_report(page)

        # Version footer: branded app name + version chip (version from core, the
        # single source of truth — see core.__version__)
        version_info = __version__
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
