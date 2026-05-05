"""
Directory-browser widget used on the Setup page.

Self-contained: takes a workspace ``params`` dict + state key; writes selected
paths back to ``st.session_state[state_key]``. Callers render inside their
column layout. Extracted from ``views.setup`` to keep that file under the
500-line review threshold.
"""
from pathlib import Path

import streamlit as st

from streamlit_app.utils.page_state import list_files, list_subdirectories


def _browser_keys(key_prefix: str) -> dict:
    """Bundle the per-browser session-state keys derived from ``key_prefix``."""
    return {
        "dir": f"{key_prefix}_dir",
        "nav": f"{key_prefix}_nav_to",
        "select_all": f"{key_prefix}_select_all",
        "recursive": f"{key_prefix}_recursive",
        "path_input": f"{key_prefix}_path_input",
    }


def _init_browser_state(keys: dict, state_key: str, default_dir: str, params: dict) -> None:
    """Seed session state for the selection list and current directory."""
    if state_key not in st.session_state:
        st.session_state[state_key] = params.get(state_key, [])
    if keys["dir"] not in st.session_state:
        st.session_state[keys["dir"]] = default_dir


def _apply_deferred_navigation(keys: dict) -> None:
    """Consume any pending Up/folder-button navigation request."""
    nav_key = keys["nav"]
    if not st.session_state.get(nav_key):
        return
    st.session_state[keys["dir"]] = st.session_state[nav_key]
    st.session_state[nav_key] = None
    for k in (keys["select_all"], keys["recursive"]):
        st.session_state.pop(k, None)


def _render_path_bar(key_prefix: str, keys: dict, current_dir: str) -> None:
    """Render the editable path input plus the parent-directory button."""
    col1, col2 = st.columns([6, 1])
    with col1:
        # Sync the input widget to current_dir when navigation just changed it,
        # without rebuilding the widget (rebuilding loses focus/caret state).
        path_input_key = keys["path_input"]
        if st.session_state.get(path_input_key) != current_dir:
            st.session_state[path_input_key] = current_dir
        new_path = st.text_input(
            "Path", key=path_input_key,
            label_visibility="collapsed", placeholder="Enter path...",
        )
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[keys["dir"]] = new_path
            for k in (keys["select_all"], keys["recursive"]):
                st.session_state.pop(k, None)
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("up", key=f"{key_prefix}_up_btn", help="Parent directory",
                     disabled=(parent == current_dir)):
            st.session_state[keys["nav"]] = parent
            st.rerun()


def _render_subfolder_grid(key_prefix: str, keys: dict, current_dir: str) -> None:
    """Render up to 16 sub-folders as a 4-column button grid."""
    subdirs = [Path(path) for path in list_subdirectories(current_dir)]
    if not subdirs:
        return

    st.markdown(
        "<style>.folder-nav button{font-size:0.8rem !important;"
        "padding:0.15rem 0.4rem !important;min-height:0 !important}</style>",
        unsafe_allow_html=True,
    )
    cols_per_row = 4
    max_folders = 16
    for row_start in range(0, min(len(subdirs), max_folders), cols_per_row):
        cols = st.columns(cols_per_row)
        for col_offset, col in enumerate(cols):
            folder_index = row_start + col_offset
            if folder_index >= len(subdirs) or folder_index >= max_folders:
                continue
            folder = subdirs[folder_index]
            with col:
                if st.button(folder.name, key=f"{key_prefix}_folder_{folder_index}",
                             help=str(folder)):
                    st.session_state[keys["nav"]] = str(folder)
                    st.rerun()
    if len(subdirs) > max_folders:
        st.caption(f"... and {len(subdirs) - max_folders} more folders")


def _list_browser_files(current_dir: str, extensions: list, recursive: bool):
    """Return ``(file_labels, label_to_path)`` for the active directory."""
    all_files = [Path(path) for path in list_files(current_dir, tuple(extensions), recursive)]
    current_path = Path(current_dir)
    file_labels = ([str(f.relative_to(current_path)) for f in all_files]
                   if recursive else [f.name for f in all_files])
    label_to_path = {lbl: str(f) for lbl, f in zip(file_labels, all_files)}
    return file_labels, label_to_path


def _selections_outside_scope(state_key: str, file_labels, label_to_path,
                              current_path: Path, recursive: bool) -> list:
    """Return previously-selected paths that fall outside the current view."""
    current_scope_paths = {label_to_path[lbl] for lbl in file_labels}
    selections = st.session_state.get(state_key, [])
    if recursive:
        return [p for p in selections if p not in current_scope_paths]
    return [p for p in selections if Path(p).parent != current_path]


def _render_file_picker(key_prefix: str, keys: dict, state_key: str,
                       current_dir: str, extensions: list, recursive: bool) -> None:
    """Render the filter input, multiselect, and persist the selection list."""
    file_labels, label_to_path = _list_browser_files(current_dir, extensions, recursive)
    filter_text = st.text_input("Filter files", key=f"{key_prefix}_filter",
                                placeholder="Type to filter filenames...")
    filtered_labels = (
        [lbl for lbl in file_labels if filter_text.lower() in lbl.lower()]
        if filter_text else file_labels
    )

    if not filtered_labels:
        if filter_text:
            st.info(f"No files match filter '{filter_text}'")
        else:
            st.info(f"No {extensions[0]} files found")
        return

    all_selected_paths = set(st.session_state.get(state_key, []))
    other_dir_selections = _selections_outside_scope(
        state_key, file_labels, label_to_path, Path(current_dir), recursive,
    )

    valid_defaults = [lbl for lbl in filtered_labels if label_to_path[lbl] in all_selected_paths]
    all_in_scope_selected = set(valid_defaults) == set(filtered_labels) and len(valid_defaults) > 0

    select_all = st.checkbox(
        f"Select all {len(filtered_labels)} files"
        + (f" (filtered from {len(file_labels)})" if filter_text else ""),
        value=all_in_scope_selected, key=keys["select_all"],
    )

    if select_all:
        selected_labels = filtered_labels
    else:
        selected_labels = st.multiselect(
            "Select files" + (" (recursive):" if recursive else ":"),
            options=filtered_labels,
            default=valid_defaults if not all_in_scope_selected else [],
            key=f"{key_prefix}_files",
        )

    current_selected_paths = [label_to_path[lbl] for lbl in selected_labels]
    merged_selections = list(dict.fromkeys(other_dir_selections + current_selected_paths))
    st.session_state[state_key] = merged_selections


def _render_clear_all(key_prefix: str, keys: dict, state_key: str) -> None:
    """Render the 'Clear all selections' button if anything is selected."""
    all_selected = st.session_state.get(state_key, [])
    if not all_selected:
        return
    if st.button(f"Clear all selections ({len(all_selected)} files)",
                 key=f"{key_prefix}_clear_all", width="stretch"):
        st.session_state[state_key] = []
        st.session_state.pop(keys["select_all"], None)
        st.rerun()


def render_directory_browser(
    key_prefix: str, extensions: list, state_key: str,
    default_dir: str = "/", params: dict | None = None,
):
    """Render a directory browser + file picker.

    Args:
        key_prefix: Unique prefix for widget keys (lets multiple browsers coexist).
        extensions: File extensions to list (e.g. ['.idXML', '.idxml']).
        state_key: ``st.session_state`` key that holds the selected paths list.
        default_dir: Initial directory when state_key is unset.
        params: Workspace params dict (for initial selections + default dir).
    """
    params = params or {}
    keys = _browser_keys(key_prefix)

    _init_browser_state(keys, state_key, default_dir, params)
    _apply_deferred_navigation(keys)

    current_dir = st.session_state[keys["dir"]]
    _render_path_bar(key_prefix, keys, current_dir)
    _render_subfolder_grid(key_prefix, keys, current_dir)

    recursive = st.checkbox("Include subdirectories", value=False,
                            key=keys["recursive"],
                            help="Scan all subdirectories for matching files")
    _render_file_picker(key_prefix, keys, state_key, current_dir, extensions, recursive)
    _render_clear_all(key_prefix, keys, state_key)