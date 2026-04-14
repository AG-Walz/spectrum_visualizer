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
    dir_key = f"{key_prefix}_dir"
    nav_key = f"{key_prefix}_nav_to"
    select_all_key = f"{key_prefix}_select_all"
    recursive_key = f"{key_prefix}_recursive"

    if state_key not in st.session_state:
        st.session_state[state_key] = params.get(state_key, [])

    if dir_key not in st.session_state:
        st.session_state[dir_key] = default_dir

    # Deferred navigation (from subfolder buttons / up button)
    if nav_key in st.session_state and st.session_state[nav_key]:
        st.session_state[dir_key] = st.session_state[nav_key]
        st.session_state[nav_key] = None
        for k in (select_all_key, recursive_key):
            st.session_state.pop(k, None)

    current_dir = st.session_state[dir_key]

    # --- Path bar + up button ---
    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path", value=current_dir,
            key=f"{key_prefix}_path_input_{abs(hash(current_dir))}",
            label_visibility="collapsed", placeholder="Enter path...",
        )
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[dir_key] = new_path
            for k in (select_all_key, recursive_key):
                st.session_state.pop(k, None)
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("up", key=f"{key_prefix}_up_btn", help="Parent directory",
                     disabled=(parent == current_dir)):
            st.session_state[nav_key] = parent
            st.rerun()

    # --- Subfolder navigation grid ---
    subdirs = [Path(d) for d in list_subdirectories(current_dir)]
    if subdirs:
        st.markdown(
            "<style>.folder-nav button{font-size:0.8rem !important;"
            "padding:0.15rem 0.4rem !important;min-height:0 !important}</style>",
            unsafe_allow_html=True,
        )
        cols_per_row = 4
        max_folders = 16
        for i in range(0, min(len(subdirs), max_folders), cols_per_row):
            cols = st.columns(cols_per_row)
            for j, col in enumerate(cols):
                idx = i + j
                if idx < len(subdirs) and idx < max_folders:
                    folder = subdirs[idx]
                    with col:
                        if st.button(folder.name, key=f"{key_prefix}_folder_{idx}",
                                     help=str(folder)):
                            st.session_state[nav_key] = str(folder)
                            st.rerun()
        if len(subdirs) > max_folders:
            st.caption(f"... and {len(subdirs) - max_folders} more folders")

    # --- File list + filter + multiselect ---
    recursive = st.checkbox("Include subdirectories", value=False,
                            key=recursive_key,
                            help="Scan all subdirectories for matching files")

    all_files = [Path(f) for f in list_files(current_dir, tuple(extensions), recursive)]
    current_path = Path(current_dir)
    file_labels = ([str(f.relative_to(current_path)) for f in all_files]
                   if recursive else [f.name for f in all_files])
    label_to_path = {lbl: str(f) for lbl, f in zip(file_labels, all_files)}

    filter_text = st.text_input("Filter files", key=f"{key_prefix}_filter",
                                placeholder="Type to filter filenames...")
    filtered_labels = (
        [lbl for lbl in file_labels if filter_text.lower() in lbl.lower()]
        if filter_text else file_labels
    )

    if filtered_labels:
        all_selected_paths = set(st.session_state.get(state_key, []))
        current_scope_paths = {label_to_path[lbl] for lbl in file_labels}
        if recursive:
            other_dir_selections = [p for p in st.session_state.get(state_key, [])
                                    if p not in current_scope_paths]
        else:
            other_dir_selections = [p for p in st.session_state.get(state_key, [])
                                    if Path(p).parent != current_path]

        valid_defaults = [lbl for lbl in filtered_labels if label_to_path[lbl] in all_selected_paths]
        all_in_scope_selected = set(valid_defaults) == set(filtered_labels) and len(valid_defaults) > 0

        select_all = st.checkbox(
            f"Select all {len(filtered_labels)} files"
            + (f" (filtered from {len(file_labels)})" if filter_text else ""),
            value=all_in_scope_selected, key=select_all_key,
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

    else:
        if filter_text:
            st.info(f"No files match filter '{filter_text}'")
        else:
            st.info(f"No {extensions[0]} files found")

    all_selected = st.session_state.get(state_key, [])
    if all_selected:
        if st.button(f"Clear all selections ({len(all_selected)} files)",
                     key=f"{key_prefix}_clear_all", width="stretch"):
            st.session_state[state_key] = []
            st.session_state.pop(select_all_key, None)
            st.rerun()
