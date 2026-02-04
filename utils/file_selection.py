"""
File selection utilities for directory browsing with substring filtering.
"""
import streamlit as st
from pathlib import Path
from typing import List, Optional, Tuple


def get_files_in_directory(
    directory: str,
    extensions: List[str],
    substring_filter: str = ""
) -> List[Path]:
    """
    Get files in a directory matching extensions and optional substring filter.

    Args:
        directory: Path to directory to search
        extensions: List of extensions to match (e.g., ['.mzML', '.mzml'])
        substring_filter: Optional substring to filter filenames

    Returns:
        List of matching file paths sorted by name
    """
    dir_path = Path(directory)
    if not dir_path.exists() or not dir_path.is_dir():
        return []

    files = []
    for ext in extensions:
        files.extend(dir_path.glob(f"*{ext}"))

    # Apply substring filter if provided
    if substring_filter:
        filter_lower = substring_filter.lower()
        files = [f for f in files if filter_lower in f.name.lower()]

    return sorted(files, key=lambda x: x.name)


def get_subdirectories(directory: str) -> List[Path]:
    """Get immediate subdirectories of a directory."""
    dir_path = Path(directory)
    if not dir_path.exists() or not dir_path.is_dir():
        return []
    return sorted([d for d in dir_path.iterdir() if d.is_dir()], key=lambda x: x.name)


def directory_browser(
    label: str,
    key_prefix: str,
    base_path: str = "/mnt/volume",
    extensions: List[str] = None,
    help_text: str = None
) -> Tuple[str, List[Path]]:
    """
    Render a directory browser with navigation and file filtering.

    Args:
        label: Label for the browser section
        key_prefix: Unique prefix for session state keys
        base_path: Starting directory path
        extensions: File extensions to show (if None, shows directories only)
        help_text: Optional help text to display

    Returns:
        Tuple of (current_directory, selected_files)
    """
    # Initialize session state
    current_dir_key = f"{key_prefix}_current_dir"
    if current_dir_key not in st.session_state:
        st.session_state[current_dir_key] = base_path

    current_dir = st.session_state[current_dir_key]

    st.subheader(label)
    if help_text:
        st.caption(help_text)

    # Path display and navigation
    col1, col2 = st.columns([4, 1])
    with col1:
        st.text_input(
            "Current path",
            value=current_dir,
            key=f"{key_prefix}_path_display",
            disabled=True,
            label_visibility="collapsed"
        )
    with col2:
        if st.button("⬆️ Up", key=f"{key_prefix}_up_btn"):
            parent = str(Path(current_dir).parent)
            st.session_state[current_dir_key] = parent
            st.rerun()

    # Subdirectory navigation
    subdirs = get_subdirectories(current_dir)
    if subdirs:
        subdir_names = ["(Select folder)"] + [d.name for d in subdirs]
        selected_subdir = st.selectbox(
            "Navigate to folder",
            options=subdir_names,
            key=f"{key_prefix}_subdir_select",
            label_visibility="collapsed"
        )
        if selected_subdir != "(Select folder)":
            new_path = str(Path(current_dir) / selected_subdir)
            st.session_state[current_dir_key] = new_path
            st.rerun()

    selected_files = []

    # File selection if extensions provided
    if extensions:
        # Substring filter
        filter_text = st.text_input(
            "Filter by filename",
            key=f"{key_prefix}_filter",
            placeholder="Type to filter files..."
        )

        # Get and display files
        files = get_files_in_directory(current_dir, extensions, filter_text)

        if files:
            # Multiselect for files
            file_names = [f.name for f in files]
            selected_names = st.multiselect(
                f"Select files ({len(files)} found)",
                options=file_names,
                key=f"{key_prefix}_file_select"
            )
            selected_files = [f for f in files if f.name in selected_names]

            if selected_files:
                st.success(f"{len(selected_files)} file(s) selected")
        else:
            st.info(f"No {', '.join(extensions)} files found in this directory")

    return current_dir, selected_files


def reference_filter_widget(
    files: List[Path],
    key_prefix: str
) -> Tuple[List[Path], List[Path]]:
    """
    Widget to separate files into reference and query based on substring.

    Args:
        files: List of file paths to categorize
        key_prefix: Unique prefix for session state keys

    Returns:
        Tuple of (reference_files, query_files)
    """
    if not files:
        return [], []

    st.markdown("**Reference File Selection**")
    st.caption("Enter a substring to identify reference files (e.g., 'Ref', 'Standard', 'Synpep')")

    ref_substring = st.text_input(
        "Reference identifier",
        key=f"{key_prefix}_ref_substring",
        placeholder="e.g., Synpep"
    )

    reference_files = []
    query_files = []

    for f in files:
        if ref_substring and ref_substring.lower() in f.name.lower():
            reference_files.append(f)
        else:
            query_files.append(f)

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**Reference files ({len(reference_files)})**")
        for f in reference_files[:5]:  # Show first 5
            st.caption(f"📊 {f.name}")
        if len(reference_files) > 5:
            st.caption(f"... and {len(reference_files) - 5} more")

    with col2:
        st.markdown(f"**Query files ({len(query_files)})**")
        for f in query_files[:5]:  # Show first 5
            st.caption(f"📄 {f.name}")
        if len(query_files) > 5:
            st.caption(f"... and {len(query_files) - 5} more")

    return reference_files, query_files
