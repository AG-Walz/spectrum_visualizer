"""Shared file selection UI components for workspaces.

Provides reusable directory browsing and file selection widgets to avoid code duplication
between Spectrum Analysis and Dataset Filter workspaces.
"""
import streamlit as st
from pathlib import Path
import os
from typing import List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def get_files_with_substring(directory_path: str, substring: str, extensions: List[str]) -> List[str]:
    """Get all files in directory and subdirectories that contain substring and match extensions.
    
    Args:
        directory_path: Root directory to search
        substring: Filter files by this substring (empty string = all files)
        extensions: List of file extensions to include (e.g., ['.mzML', '.idXML'])
        
    Returns:
        Sorted list of file paths matching criteria
    """
    try:
        path = Path(directory_path)
        if not path.exists():
            return []
        
        files = []
        for file_path in path.rglob("*"):
            if file_path.is_file():
                if substring.lower() in file_path.name.lower() or substring == "":
                    if file_path.suffix.lower() in [ext.lower() for ext in extensions]:
                        files.append(str(file_path))
        
        return sorted(files)
    except Exception as e:
        logger.error(f"Error accessing directory {directory_path}: {e}")
        return []


def _get_subdirectories(directory: str) -> List[Path]:
    """Get sorted list of subdirectories (cached for performance)."""
    try:
        path = Path(directory)
        if not path.exists() or not path.is_dir():
            return []
        subdirs = [d for d in path.iterdir() if d.is_dir()]
        return sorted(subdirs, key=lambda x: x.name.lower())
    except Exception as e:
        logger.error(f"Error listing subdirectories: {e}")
        return []


def render_directory_browser(
    session_key: str,
    default_path: str = "/mnt/volume/",
    show_subdirs: bool = True,
    compact: bool = False,
    allow_manual_input: bool = True
) -> str:
    """Render a directory browser widget.
    
    Args:
        session_key: Unique key for session state (e.g., 'filter_directory', 'mzml_directory')
        default_path: Starting directory path
        show_subdirs: Whether to show subdirectory navigation buttons
        compact: If True, uses simpler layout without parent/refresh buttons
        allow_manual_input: If True, allows user to type or paste directory path
        
    Returns:
        Currently selected directory path
    """
    # Initialize session state
    if session_key not in st.session_state:
        st.session_state[session_key] = default_path
    
    current_dir = st.session_state[session_key]
    
    # Display current directory with back button
    col1, col2 = st.columns([5, 1])
    with col1:
        if allow_manual_input:
            # Allow user to type or paste path
            manual_path = st.text_input(
                "Current Directory (editable - type or paste path)",
                value=current_dir,
                key=f"{session_key}_manual_input",
                help="You can type or paste a directory path here"
            )
            # Update session state if user changed the path
            if manual_path != current_dir and manual_path:
                if Path(manual_path).exists() and Path(manual_path).is_dir():
                    st.session_state[session_key] = manual_path
                    st.rerun()
                elif manual_path:  # User entered something but it's invalid
                    st.warning(f"⚠️ Directory does not exist: {manual_path}")
        else:
            # Non-editable display
            st.text_input(
                "Current Directory",
                value=current_dir,
                key=f"{session_key}_display",
                disabled=True,
                help="Currently selected directory"
            )
    with col2:
        # Add spacing to align button with text input (label + input height)
        st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
        parent_dir = str(Path(current_dir).parent)
        # Only show back button if we're not at root and parent exists
        if parent_dir != current_dir and Path(parent_dir).exists():
            if st.button("⬆️ Back", key=f"{session_key}_back_btn", use_container_width=True, help="Go to parent directory"):
                st.session_state[session_key] = parent_dir
                st.rerun()
    
    # Show subdirectories as clickable buttons in a compact grid
    if show_subdirs:
        try:
            if os.path.exists(current_dir) and os.path.isdir(current_dir):
                subdirs = _get_subdirectories(current_dir)
                
                if subdirs:
                    st.markdown("**Subdirectories:**")
                    
                    # Display folders as compact buttons in a grid (4 columns)
                    num_cols = 4
                    for i in range(0, len(subdirs), num_cols):
                        cols = st.columns(num_cols)
                        for j, col in enumerate(cols):
                            idx = i + j
                            if idx < len(subdirs):
                                subdir = subdirs[idx]
                                with col:
                                    # Truncate long names - doubled from 15 to 30 characters
                                    display_name = subdir.name if len(subdir.name) <= 30 else subdir.name[:27] + "..."
                                    if st.button(
                                        f"📁 {display_name}",
                                        key=f"{session_key}_dir_{idx}",
                                        use_container_width=True,
                                        help=subdir.name
                                    ):
                                        st.session_state[session_key] = str(subdir)
                                        st.rerun()
                else:
                    st.info("No subdirectories found")
        except Exception as e:
            st.error(f"Error accessing directory: {e}")
    
    return current_dir


def render_file_filter_and_selection(
    directory: str,
    file_type: str,
    extensions: List[str],
    session_key: str,
    default_substring: str = "",
    include_subdirs: bool = True,
    show_select_all: bool = True,
    multiselect: bool = True
) -> Tuple[List[str], str, bool]:
    """Render file filtering and selection UI.
    
    Args:
        directory: Directory to search for files
        file_type: Display name for file type (e.g., "mzML", "idXML")
        extensions: List of file extensions (e.g., ['.mzML'], ['.idXML', '.mzid'])
        session_key: Unique key for widgets
        default_substring: Default value for substring filter
        include_subdirs: Whether to search subdirectories
        show_select_all: Whether to show "select all" checkbox
        multiselect: If True, allows multiple file selection; if False, single select
        
    Returns:
        Tuple of (selected_files, substring_filter, include_subdirs_flag)
    """
    st.markdown(f"**{file_type} Files**")
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        substring = st.text_input(
            f"Filter {file_type} files",
            value=default_substring,
            key=f"{session_key}_substring",
            placeholder=f"e.g., HNSCC, Sample1",
            help="Filter files by substring in filename"
        )
    
    with col2:
        include_subdirs = st.checkbox(
            "Subdirs",
            value=include_subdirs,
            key=f"{session_key}_subdirs",
            help="Include subdirectories in search"
        )
    
    # Find files
    if include_subdirs:
        files_found = get_files_with_substring(directory, substring, extensions)
    else:
        # Only current directory
        files_found = sorted([
            str(f) for f in Path(directory).iterdir()
            if f.is_file() and f.suffix.lower() in [ext.lower() for ext in extensions]
            and (substring.lower() in f.name.lower() or substring == "")
        ])
    
    if not files_found:
        st.warning(f"No {file_type} files found")
        return [], substring, include_subdirs
    
    # Show file count
    count_msg = f"Found {len(files_found)} {file_type} file(s)"
    if len(files_found) >= 10000:
        st.warning(f"⚠️ {count_msg} (search limit reached)")
    else:
        st.info(count_msg)
    
    selected_files = []
    
    # Selection UI
    if show_select_all:
        if st.checkbox(f"Select all {len(files_found)} files", key=f"{session_key}_select_all"):
            selected_files = files_found
        elif multiselect:
            # Show multiselect for manual selection
            from shared_imports import create_file_options_with_sizes, map_display_names_to_paths
            file_options = create_file_options_with_sizes(files_found)
            selected_display = st.multiselect(
                f"Select {file_type} files",
                options=file_options,
                key=f"{session_key}_multiselect",
                help=f"Choose specific {file_type} files"
            )
            selected_files = map_display_names_to_paths(selected_display, files_found)
    elif multiselect:
        # Only multiselect without select all
        from shared_imports import create_file_options_with_sizes, map_display_names_to_paths
        file_options = create_file_options_with_sizes(files_found)
        selected_display = st.multiselect(
            f"Select {file_type} files",
            options=file_options,
            key=f"{session_key}_multiselect",
            help=f"Choose specific {file_type} files"
        )
        selected_files = map_display_names_to_paths(selected_display, files_found)
    else:
        # Single select radio buttons
        file_names = [Path(f).name for f in files_found]
        selected_name = st.radio(
            f"Select {file_type} file",
            options=file_names,
            key=f"{session_key}_radio"
        )
        if selected_name:
            selected_files = [f for f in files_found if Path(f).name == selected_name]
    
    return selected_files, substring, include_subdirs


def render_folder_multi_selector(
    session_key: str,
    default_path: str = "/mnt/volume/"
) -> List[str]:
    """Render folder selection UI allowing multiple folders.
    
    Args:
        session_key: Unique key for session state
        default_path: Starting directory path
        
    Returns:
        List of selected folder paths
    """
    # Initialize session state
    if f"{session_key}_folders" not in st.session_state:
        st.session_state[f"{session_key}_folders"] = []
    
    # Directory browser for adding folders
    st.markdown("**Browse and Add Folders:**")
    current_folder = render_directory_browser(
        session_key=f"{session_key}_browser",
        default_path=default_path,
        show_subdirs=True,
        compact=True
    )
    
    # Add folder button
    if current_folder and st.button("➕ Add Folder", key=f"{session_key}_add_btn", use_container_width=True):
        if current_folder not in st.session_state[f"{session_key}_folders"]:
            st.session_state[f"{session_key}_folders"].append(current_folder)
            st.rerun()
    
    # Show selected folders
    if st.session_state[f"{session_key}_folders"]:
        st.markdown(f"**Selected Folders ({len(st.session_state[f'{session_key}_folders'])}):**")
        
        folders_to_remove = []
        for idx, folder in enumerate(st.session_state[f"{session_key}_folders"]):
            col_folder, col_remove = st.columns([6, 1])
            with col_folder:
                st.text(f"📁 {folder}")
            with col_remove:
                if st.button("❌", key=f"{session_key}_remove_{idx}", help="Remove this folder"):
                    folders_to_remove.append(folder)
        
        # Remove folders after iteration
        for folder in folders_to_remove:
            st.session_state[f"{session_key}_folders"].remove(folder)
            st.rerun()
    
    return st.session_state[f"{session_key}_folders"]