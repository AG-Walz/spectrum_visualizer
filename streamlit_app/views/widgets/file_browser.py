"""Directory-browser + file-picker widget for the Setup page: breadcrumbs, sub-folder grid, glob-filtered multiselect — selections persist in st.session_state[state_key]."""
from pathlib import Path

import streamlit as st

from streamlit_app.utils.page_state import list_files, list_subdirectories


def _browser_keys(key_prefix: str) -> dict:
    """Return the per-browser session_state widget keys derived from key_prefix."""
    return {
        "dir": f"{key_prefix}_dir",
        "nav": f"{key_prefix}_nav_to",
        "select_all": f"{key_prefix}_select_all",
        "recursive": f"{key_prefix}_recursive",
        "path_input": f"{key_prefix}_path_input",
    }


def _init_browser_state(keys: dict, state_key: str, default_dir: str, params: dict) -> None:
    """Seed session_state slots (selection list, current directory, path input) on first run."""
    if state_key not in st.session_state:
        st.session_state[state_key] = params.get(state_key, [])
    if keys["dir"] not in st.session_state:
        st.session_state[keys["dir"]] = default_dir
    if keys["path_input"] not in st.session_state:
        st.session_state[keys["path_input"]] = st.session_state[keys["dir"]]


def _apply_deferred_navigation(keys: dict) -> None:
    """Consume any pending Up/folder-button navigation request before widgets render this run."""
    # Guard: bail when no navigation was requested this run
    nav_key = keys["nav"]
    if not st.session_state.get(nav_key):
        return
    # Apply: update current directory + visible path input
    new_dir = st.session_state[nav_key]
    st.session_state[keys["dir"]] = new_dir
    # Sync the path input to the new directory only on real navigation events
    # (Up/folder buttons) — not on every rerun — so in-progress typing in the
    # path input is preserved.
    st.session_state[keys["path_input"]] = new_dir
    # Cleanup: clear nav slot + stale select-all/recursive widget state
    st.session_state[nav_key] = None
    for k in (keys["select_all"], keys["recursive"]):
        st.session_state.pop(k, None)


def _render_path_bar(key_prefix: str, keys: dict, current_dir: str) -> None:
    """Render the editable path input + the Up (parent directory) button."""
    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path", key=keys["path_input"],
            label_visibility="collapsed", placeholder="Enter path...",
        )
        # Resolve to absolute so the Up button (disabled when parent == current_dir)
        # works after typing a relative path like "." or "subdir/". Route through
        # the deferred-nav slot — Streamlit forbids writing a widget-bound key
        # after the widget instantiates this run.
        if new_path and new_path != current_dir:
            resolved = str(Path(new_path).expanduser().resolve())
            if Path(resolved).is_dir() and resolved != current_dir:
                st.session_state[keys["nav"]] = resolved
                st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("up", key=f"{key_prefix}_up_btn", help="Parent directory",
                     disabled=(parent == current_dir)):
            st.session_state[keys["nav"]] = parent
            st.rerun()


def _render_subfolder_grid(key_prefix: str, keys: dict, current_dir: str) -> None:
    """Render up to 16 sub-folders as a 4-column grid of navigation buttons."""
    # List sub-folders; bail early when the directory has none
    subdirs = [Path(path) for path in list_subdirectories(current_dir)]
    if not subdirs:
        return

    # Inject compact-button styling for the folder grid
    st.markdown(
        "<style>.folder-nav button{font-size:0.8rem !important;"
        "padding:0.15rem 0.4rem !important;min-height:0 !important}</style>",
        unsafe_allow_html=True,
    )
    # Render the 4×4 button grid; each button defers navigation via the nav slot
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
    # Overflow note: indicate how many folders were truncated
    if len(subdirs) > max_folders:
        st.caption(f"... and {len(subdirs) - max_folders} more folders")


def _list_browser_files(current_dir: str, extensions: list, recursive: bool):
    """List files in current_dir matching extensions; return (labels, label→path)."""
    all_files = [Path(path) for path in list_files(current_dir, tuple(extensions), recursive)]
    current_path = Path(current_dir)
    file_labels = ([str(f.relative_to(current_path)) for f in all_files]
                   if recursive else [f.name for f in all_files])
    label_to_path = {lbl: str(f) for lbl, f in zip(file_labels, all_files)}
    return file_labels, label_to_path


def _selections_outside_scope(state_key: str, file_labels, label_to_path,
                              current_path: Path, recursive: bool) -> list:
    """Return previously-selected paths that fall outside the current directory view (so they survive nav)."""
    current_scope_paths = {label_to_path[lbl] for lbl in file_labels}
    selections = st.session_state.get(state_key, [])
    if recursive:
        return [p for p in selections if p not in current_scope_paths]
    return [p for p in selections if Path(p).parent != current_path]


def _render_file_picker(key_prefix: str, keys: dict, state_key: str,
                       current_dir: str, extensions: list, recursive: bool) -> None:
    """Render the filter input + select-all + multiselect and persist picks to session_state[state_key]."""
    # List + filter files in the current view by the substring filter input
    file_labels, label_to_path = _list_browser_files(current_dir, extensions, recursive)
    filter_text = st.text_input("Filter files", key=f"{key_prefix}_filter",
                                placeholder="Type to filter filenames...")
    filtered_labels = (
        [lbl for lbl in file_labels if filter_text.lower() in lbl.lower()]
        if filter_text else file_labels
    )

    # Empty-state: show info text and bail before rendering pickers
    if not filtered_labels:
        if filter_text:
            st.info(f"No files match filter '{filter_text}'")
        else:
            st.info(f"No {extensions[0]} files found")
        return

    # Compute defaults: preserve out-of-scope picks, mark in-scope ones as selected
    all_selected_paths = set(st.session_state.get(state_key, []))
    other_dir_selections = _selections_outside_scope(
        state_key, file_labels, label_to_path, Path(current_dir), recursive,
    )

    valid_defaults = [lbl for lbl in filtered_labels if label_to_path[lbl] in all_selected_paths]
    all_in_scope_selected = set(valid_defaults) == set(filtered_labels) and len(valid_defaults) > 0

    # Render select-all checkbox + the multiselect (or skip multiselect when all selected)
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

    # Persist: merge with out-of-scope selections, dedupe while keeping order
    current_selected_paths = [label_to_path[lbl] for lbl in selected_labels]
    merged_selections = list(dict.fromkeys(other_dir_selections + current_selected_paths))
    st.session_state[state_key] = merged_selections


def _render_clear_all(key_prefix: str, keys: dict, state_key: str) -> None:
    """Render the 'Clear all selections' button when at least one file is selected."""
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
    """Render the directory browser + file picker on the Setup page; persists selections to session_state[state_key].

    key_prefix namespaces widget keys; extensions filters listed files; params seeds initial dir + selection."""
    # Bootstrap: initialise state + consume any pending navigation request
    params = params or {}
    keys = _browser_keys(key_prefix)

    _init_browser_state(keys, state_key, default_dir, params)
    _apply_deferred_navigation(keys)

    # Render navigation widgets: path bar + sub-folder grid
    current_dir = st.session_state[keys["dir"]]
    _render_path_bar(key_prefix, keys, current_dir)
    _render_subfolder_grid(key_prefix, keys, current_dir)

    # Render file picker (with recursive toggle) and the clear-all action
    recursive = st.checkbox("Include subdirectories", value=False,
                            key=keys["recursive"],
                            help="Scan all subdirectories for matching files")
    _render_file_picker(key_prefix, keys, state_key, current_dir, extensions, recursive)
    _render_clear_all(key_prefix, keys, state_key)