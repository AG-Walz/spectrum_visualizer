"""
Peptide Mirror Analysis Workspace

Interactive workspace for comparing query spectra against reference spectra.
"""
import streamlit as st
from pathlib import Path
import json
import logging

logger = logging.getLogger(__name__)

# State key prefix for this workspace
STATE_PREFIX = "peptide_mirror"


def get_state_key(key: str) -> str:
    """Get namespaced session state key."""
    return f"{STATE_PREFIX}.{key}"


def init_session_state():
    """Initialize session state with defaults."""
    defaults = {
        'idxml_dir': '/mnt/volume',
        'mzml_dir': '/mnt/volume',
        'selected_idxmls': [],
        'selected_mzmls': [],
        'target_peptides': '',
        'reference_substring': '',
        'ion_types': ['b', 'y'],
        'mz_tolerance': 0.02,
        'search_results': None,
        'loaded_variants': None,
        'current_step': 1,
    }

    for key, default in defaults.items():
        state_key = get_state_key(key)
        if state_key not in st.session_state:
            st.session_state[state_key] = default


def render_file_selection():
    """Render file selection step."""
    st.header("Step 1: File Selection")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(
            key_prefix=f"{STATE_PREFIX}_idxml",
            extensions=['.idXML', '.idxml'],
            state_key='selected_idxmls'
        )

    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(
            key_prefix=f"{STATE_PREFIX}_mzml",
            extensions=['.mzML', '.mzml'],
            state_key='selected_mzmls'
        )

    # Show selection summary
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    if idxmls and mzmls:
        st.success(f"Selected {len(idxmls)} idXML file(s) and {len(mzmls)} mzML file(s)")
        return True
    else:
        if not idxmls:
            st.warning("Please select at least one idXML file")
        if not mzmls:
            st.warning("Please select at least one mzML file")
        return False


def render_directory_browser(key_prefix: str, extensions: list, state_key: str):
    """Render a directory browser with file selection."""
    # State keys
    dir_key = f"{key_prefix}_dir"
    nav_key = f"{key_prefix}_nav_to"  # Pending navigation target
    select_all_key = f"{key_prefix}_select_all"

    # Initialize directory state
    if dir_key not in st.session_state:
        st.session_state[dir_key] = "/mnt/volume"

    # Handle pending navigation (set by button callbacks)
    if nav_key in st.session_state and st.session_state[nav_key]:
        st.session_state[dir_key] = st.session_state[nav_key]
        st.session_state[nav_key] = None
        # Reset select all when changing directory
        if select_all_key in st.session_state:
            del st.session_state[select_all_key]

    current_dir = st.session_state[dir_key]

    # Path input with up button
    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path",
            value=current_dir,
            key=f"{key_prefix}_path_input_{hash(current_dir) % 10000}",  # Dynamic key forces refresh
            label_visibility="collapsed",
            placeholder="Enter path..."
        )
        # Handle manual path entry
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[dir_key] = new_path
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("⬆️", key=f"{key_prefix}_up_btn", help="Parent directory", disabled=(parent == current_dir)):
            st.session_state[nav_key] = parent
            st.rerun()

    # Show clickable subdirectories
    try:
        subdirs = sorted([d for d in Path(current_dir).iterdir() if d.is_dir()])
        if subdirs:
            st.caption("Click folder to navigate:")
            # Display folders in a grid
            cols_per_row = 3
            max_folders = 15
            for i in range(0, min(len(subdirs), max_folders), cols_per_row):
                cols = st.columns(cols_per_row)
                for j, col in enumerate(cols):
                    idx = i + j
                    if idx < len(subdirs) and idx < max_folders:
                        folder = subdirs[idx]
                        folder_path = str(folder)
                        with col:
                            if st.button(
                                f"📁 {folder.name}",
                                key=f"{key_prefix}_folder_{idx}",
                                width="stretch"
                            ):
                                st.session_state[nav_key] = folder_path
                                st.rerun()
            if len(subdirs) > max_folders:
                st.caption(f"... and {len(subdirs) - max_folders} more folders")
    except PermissionError:
        st.warning("Permission denied")
    except Exception as e:
        st.error(f"Cannot access: {e}")

    st.markdown("---")

    # Get matching files
    try:
        all_files = []
        for ext in extensions:
            all_files.extend(Path(current_dir).glob(f"*{ext}"))
        all_files = sorted(all_files, key=lambda x: x.name)

        # Filter input
        filter_text = st.text_input(
            "Filter files",
            key=f"{key_prefix}_filter",
            placeholder="Type to filter filenames..."
        )

        # Apply filter
        if filter_text:
            filter_lower = filter_text.lower()
            filtered_files = [f for f in all_files if filter_lower in f.name.lower()]
        else:
            filtered_files = all_files

        if filtered_files:
            file_names = [f.name for f in filtered_files]

            # Get currently selected files
            current_selected = st.session_state.get(get_state_key(state_key), [])
            current_names = [Path(f).name for f in current_selected if Path(f).exists()]
            valid_defaults = [n for n in current_names if n in file_names]

            # Check if all are currently selected
            all_selected = set(valid_defaults) == set(file_names) and len(valid_defaults) > 0

            # Select All checkbox
            select_all = st.checkbox(
                f"Select all {len(filtered_files)} files" + (f" (filtered from {len(all_files)})" if filter_text else ""),
                value=all_selected,
                key=select_all_key
            )

            if select_all:
                # All files selected
                selected_names = file_names
            else:
                # Manual selection via multiselect
                selected_names = st.multiselect(
                    "Or select specific files:",
                    options=file_names,
                    default=valid_defaults if not all_selected else [],
                    key=f"{key_prefix}_files"
                )

            # Update state with full paths
            selected_paths = [str(f) for f in filtered_files if f.name in selected_names]
            st.session_state[get_state_key(state_key)] = selected_paths

            if selected_paths:
                st.success(f"✅ {len(selected_paths)} file(s) selected")
        else:
            if filter_text:
                st.info(f"No files match filter '{filter_text}'")
            else:
                st.info(f"No {extensions[0]} files in this directory")

    except Exception as e:
        st.error(f"Error listing files: {e}")


def render_peptide_input():
    """Render peptide input step."""
    st.header("Step 2: Target Peptides")

    col1, col2 = st.columns([2, 1])

    with col1:
        peptide_text = st.text_area(
            "Enter peptide sequences (one per line)",
            value=st.session_state.get(get_state_key('target_peptides'), ''),
            height=150,
            key=f"{STATE_PREFIX}_peptide_input",
            placeholder="PEPTIDESEQ\nANOTHERPEPTIDE\n..."
        )
        st.session_state[get_state_key('target_peptides')] = peptide_text

    with col2:
        st.markdown("**Tips:**")
        st.markdown("""
        - Enter base sequences (no modifications)
        - One peptide per line
        - Modified variants will be found automatically
        - Case insensitive
        """)

        # Count entered peptides
        peptides = [p.strip() for p in peptide_text.split('\n') if p.strip()]
        if peptides:
            st.success(f"{len(peptides)} peptide(s) entered")

            # Preview button
            if st.button("🔍 Preview in idXML", key=f"{STATE_PREFIX}_preview_btn"):
                with st.spinner("Scanning idXML files..."):
                    preview_peptides_in_idxmls()

    return len(peptides) > 0


def preview_peptides_in_idxmls():
    """Preview which peptides are found in selected idXML files."""
    from core.peptide_search import get_unique_peptides_from_idxmls

    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    if not idxmls:
        st.warning("No idXML files selected")
        return

    # Get sample of peptides from files
    sample_peptides = get_unique_peptides_from_idxmls(idxmls, limit=50)

    if sample_peptides:
        st.info(f"Sample peptides from idXML files (first 50):")
        # Display in columns
        cols = st.columns(3)
        for i, pep in enumerate(sample_peptides[:30]):
            cols[i % 3].caption(pep)


def render_reference_selection():
    """Render reference file selection step."""
    st.header("Step 3: Reference Identification")

    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    if not mzmls:
        st.warning("Please select mzML files first")
        return False

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")

    ref_substring = st.text_input(
        "Reference identifier",
        value=st.session_state.get(get_state_key('reference_substring'), ''),
        key=f"{STATE_PREFIX}_ref_input",
        placeholder="e.g., Synpep"
    )
    st.session_state[get_state_key('reference_substring')] = ref_substring

    # Show categorization preview
    if ref_substring:
        reference_files = []
        query_files = []

        for mzml in mzmls:
            fname = Path(mzml).name
            if ref_substring.lower() in fname.lower():
                reference_files.append(fname)
            else:
                query_files.append(fname)

        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Reference files ({len(reference_files)})**")
            for f in reference_files[:5]:
                st.caption(f"📊 {f}")
            if len(reference_files) > 5:
                st.caption(f"... and {len(reference_files) - 5} more")

        with col2:
            st.markdown(f"**Query files ({len(query_files)})**")
            for f in query_files[:5]:
                st.caption(f"📄 {f}")
            if len(query_files) > 5:
                st.caption(f"... and {len(query_files) - 5} more")

        if reference_files and query_files:
            st.success("Reference/Query separation configured")
            return True
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")

    else:
        st.info("Without a reference identifier, spectra will be compared query-to-query")
        return True

    return False


def render_settings():
    """Render analysis settings."""
    st.header("Step 4: Settings")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Ion Types")
        all_ions = ['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium']
        current_ions = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])

        selected_ions = st.multiselect(
            "Select ion types to annotate",
            options=all_ions,
            default=[i for i in current_ions if i in all_ions],
            key=f"{STATE_PREFIX}_ion_select"
        )
        st.session_state[get_state_key('ion_types')] = selected_ions

    with col2:
        st.subheader("m/z Tolerance")
        current_tol = st.session_state.get(get_state_key('mz_tolerance'), 0.02)

        tolerance = st.number_input(
            "Tolerance (Da)",
            min_value=0.001,
            max_value=1.0,
            value=float(current_tol),
            step=0.005,
            format="%.3f",
            key=f"{STATE_PREFIX}_tol_input"
        )
        st.session_state[get_state_key('mz_tolerance')] = tolerance

    return True


def render_search():
    """Render peptide search step."""
    st.header("Step 5: Search Peptides")

    # Check prerequisites - read from both possible state keys for robustness
    # Files may be stored directly or via the widget
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])

    # Peptides: try the synced key first, then the widget key directly
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    if not peptide_text:
        # Fallback: read directly from text_area widget key
        peptide_text = st.session_state.get(f"{STATE_PREFIX}_peptide_input", '')
    peptides = [p.strip().upper() for p in peptide_text.split('\n') if p.strip()]

    # Debug info
    with st.expander("Debug: Current State", expanded=False):
        st.write(f"idXML files: {len(idxmls)}")
        st.write(f"Peptide text length: {len(peptide_text)}")
        st.write(f"Peptides parsed: {peptides}")
        if idxmls:
            st.write(f"First idXML: {idxmls[0]}")

    if not idxmls:
        st.warning("Please select idXML files first (Files tab)")
        return

    if not peptides:
        st.warning("Please enter target peptides first (Peptides tab)")
        return

    # Show what will be searched
    st.info(f"Ready to search {len(peptides)} peptide(s) in {len(idxmls)} idXML file(s)")

    if st.button("🔍 Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(idxmls, peptides)

    # Display search results
    search_results = st.session_state.get(get_state_key('search_results'))
    if search_results:
        st.markdown("---")
        st.subheader("Search Results")
        display_search_results(search_results)

        total_matches = sum(len(vg.matches) for vg in search_results.values())
        st.success(f"Found {total_matches} PSMs across {len(search_results)} peptide variants. Go to **Analysis** tab to load spectra.")
    else:
        st.caption("Search results will appear here after running the search.")


def render_analysis():
    """Render spectrum loading and analysis step."""
    st.header("Step 6: Load & Analyze Spectra")

    # Check prerequisites
    search_results = st.session_state.get(get_state_key('search_results'))
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    # Reference substring: try synced key, then widget key
    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')
    if not ref_substring:
        ref_substring = st.session_state.get(f"{STATE_PREFIX}_ref_input", '')

    if not search_results:
        st.warning("Please run peptide search first (Search tab)")
        return

    if not mzmls:
        st.warning("Please select mzML files first (Files tab)")
        return

    # Calculate expected spectra count
    total_psms = sum(len(vg.matches) for vg in search_results.values())
    loaded_variants = st.session_state.get(get_state_key('loaded_variants'))

    # Show loading status
    col1, col2 = st.columns([2, 1])

    with col1:
        if loaded_variants:
            total_loaded = sum(
                len(vg.query_spectra) + len(vg.reference_spectra)
                for vg in loaded_variants.values()
            )
            st.metric("Spectra Loaded", f"{total_loaded} / {total_psms} PSMs")
        else:
            st.metric("Spectra Loaded", f"0 / {total_psms} PSMs")

    with col2:
        # Reference info
        if ref_substring:
            ref_count = sum(1 for m in mzmls if ref_substring.lower() in Path(m).name.lower())
            st.metric("Reference Files", f"{ref_count} / {len(mzmls)}")

    st.markdown("---")

    # Loading button and progress
    if st.button("📂 Load Spectra from mzML", key=f"{STATE_PREFIX}_load_btn", type="primary"):
        run_spectrum_loading_with_progress()

    # Show loading progress placeholder
    progress_container = st.container()

    # Display loaded spectra summary
    if loaded_variants:
        st.markdown("---")
        st.subheader("Loaded Spectra Summary")
        display_loaded_spectra(loaded_variants)

        total_query = sum(len(vg.query_spectra) for vg in loaded_variants.values())
        total_ref = sum(len(vg.reference_spectra) for vg in loaded_variants.values())

        if total_query > 0 or total_ref > 0:
            st.success(f"✅ Ready for visualization! {total_query} query + {total_ref} reference spectra loaded. Go to **Results** tab.")
        else:
            st.warning("No spectra could be loaded. Check that mzML files match the idXML identifications.")
    else:
        st.caption("Click 'Load Spectra' to load matched spectra from mzML files.")


def run_peptide_search(idxmls: list, peptides: list):
    """Execute peptide search."""
    from core.peptide_search import search_peptides_in_idxmls

    progress_bar = st.progress(0, "Starting search...")
    status_text = st.empty()
    error_container = st.container()
    error_messages = []

    def progress_callback(current, total, message):
        if total > 0:
            progress_bar.progress(current / total, message)
            status_text.text(message)

    def error_callback(error_msg: str):
        """Handle errors immediately."""
        error_messages.append(error_msg)
        with error_container:
            st.warning(f"⚠️ {error_msg}")

    try:
        results = search_peptides_in_idxmls(
            idxmls, peptides,
            progress_callback=progress_callback,
            error_callback=error_callback
        )

        st.session_state[get_state_key('search_results')] = results
        progress_bar.progress(1.0, "Search complete!")
        status_text.text("Search complete!")

        if results:
            total_matches = sum(len(vg.matches) for vg in results.values())

            # Count resolved vs unresolved mzML mappings
            resolved = sum(1 for vg in results.values() for m in vg.matches if m.mzml_filename)
            unresolved = total_matches - resolved

            st.success(f"Found {total_matches} PSMs across {len(results)} peptide variants")

            if unresolved > 0:
                st.warning(f"{unresolved} PSMs could not be mapped to mzML files (missing id_merge_index)")

            # Show error summary
            if error_messages:
                with st.expander(f"⚠️ {len(error_messages)} warnings during search", expanded=False):
                    for msg in error_messages[:10]:
                        st.caption(msg)
                    if len(error_messages) > 10:
                        st.caption(f"... and {len(error_messages) - 10} more")
        else:
            st.warning("No matches found for the specified peptides")

    except Exception as e:
        st.error(f"Search failed: {e}")
        logger.exception("Peptide search error")


def run_spectrum_loading_with_progress():
    """Load spectra for search results with detailed progress tracking."""
    from core.spectrum_loader import load_matched_spectra, LoadingProgress

    search_results = st.session_state.get(get_state_key('search_results'))
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')

    if not search_results:
        st.warning("No search results to load")
        return

    mzml_paths = [Path(m) for m in mzmls]
    total_psms = sum(len(vg.matches) for vg in search_results.values())

    # Progress display
    progress_bar = st.progress(0)
    status_text = st.empty()

    # Stats display in columns
    stats_cols = st.columns(4)
    files_metric = stats_cols[0].empty()
    loaded_metric = stats_cols[1].empty()
    failed_metric = stats_cols[2].empty()
    current_metric = stats_cols[3].empty()

    # Error display area
    error_container = st.container()
    error_messages = []

    def progress_callback(progress: LoadingProgress):
        """Handle progress updates from spectrum loader."""
        # Update progress bar
        progress_bar.progress(progress.progress_fraction)

        # Update metrics
        files_metric.metric(
            "Files",
            f"{progress.files_done} / {progress.total_files}"
        )
        loaded_metric.metric(
            "Loaded",
            f"{progress.spectra_loaded} / {progress.spectra_expected}"
        )
        failed_metric.metric(
            "Not Found",
            f"{progress.spectra_not_found}"
        )

        # Current file info
        if progress.current_file:
            file_name = progress.current_file[:20] + "..." if len(progress.current_file) > 20 else progress.current_file
            current_metric.metric(
                "Current File",
                file_name,
                delta=f"{progress.current_file_spectra}/{progress.current_file_total}"
            )

        # Status text
        status_text.text(
            f"Loading {progress.current_file}: {progress.current_file_spectra}/{progress.current_file_total} spectra"
        )

    def error_callback(error_msg: str):
        """Handle errors immediately."""
        error_messages.append(error_msg)
        # Display errors as they come in
        with error_container:
            st.warning(f"⚠️ {error_msg}")

    try:
        # Make a copy of search results for loading
        import copy
        variant_groups = copy.deepcopy(search_results)

        status_text.text("Starting spectrum loading...")
        files_metric.metric("Files", f"0 / ?")
        loaded_metric.metric("Loaded", f"0 / {total_psms}")
        failed_metric.metric("Not Found", "0")
        current_metric.metric("Current File", "Starting...")

        loaded = load_matched_spectra(
            variant_groups, mzml_paths,
            reference_substring=ref_substring,
            progress_callback=progress_callback,
            error_callback=error_callback
        )

        st.session_state[get_state_key('loaded_variants')] = loaded

        # Final counts
        total_query = sum(len(vg.query_spectra) for vg in loaded.values())
        total_ref = sum(len(vg.reference_spectra) for vg in loaded.values())
        total_loaded = total_query + total_ref
        total_failed = total_psms - total_loaded

        progress_bar.progress(1.0)
        status_text.text("Loading complete!")
        files_metric.metric("Files", f"✅ Complete")
        loaded_metric.metric("Loaded", f"{total_loaded}")
        failed_metric.metric("Not Found", f"{total_failed}")
        current_metric.metric("Current File", "✅ Done")

        # Summary
        if total_loaded > 0:
            st.success(f"Loaded {total_query} query + {total_ref} reference spectra")
        else:
            st.error("No spectra could be loaded. Check mzML files match the identifications.")

        # Show error summary if any
        if error_messages:
            with st.expander(f"⚠️ {len(error_messages)} warnings during loading", expanded=False):
                for msg in error_messages[:20]:  # Limit display
                    st.caption(msg)
                if len(error_messages) > 20:
                    st.caption(f"... and {len(error_messages) - 20} more")

    except Exception as e:
        st.error(f"Loading failed: {e}")
        status_text.text(f"Error: {e}")
        logger.exception("Spectrum loading error")


def display_search_results(results: dict):
    """Display search results table."""
    from core.peptide_search import summarize_search_results

    summary_df = summarize_search_results(results)
    if not summary_df.empty:
        # Rename columns for display
        display_df = summary_df.rename(columns={
            'base_sequence': 'Peptide',
            'modification': 'Modification',
            'psm_count': 'PSMs',
            'idxml_files': 'idXML Files',
            'mzml_resolved': 'mzML Resolved',
            'unresolved': 'Unresolved'
        })
        # Drop display_label from table view
        display_df = display_df.drop(columns=['display_label'], errors='ignore')

        st.dataframe(
            display_df,
            width="stretch",
            hide_index=True,
            column_config={
                "Unresolved": st.column_config.NumberColumn(
                    "Unresolved",
                    help="PSMs that could not be mapped to mzML files"
                )
            }
        )


def display_loaded_spectra(loaded_variants: dict):
    """Display summary of loaded spectra."""
    rows = []
    for variant_key, vg in loaded_variants.items():
        rows.append({
            'Peptide': vg.display_label,
            'Query Spectra': len(vg.query_spectra),
            'Reference Spectra': len(vg.reference_spectra),
            'Total': len(vg.query_spectra) + len(vg.reference_spectra)
        })

    if rows:
        import pandas as pd
        df = pd.DataFrame(rows)
        st.dataframe(df, width="stretch", hide_index=True)


def render_results_gallery():
    """Render results gallery with mirror plots."""
    st.header("Step 7: Results & Export")

    loaded_variants = st.session_state.get(get_state_key('loaded_variants'))
    if not loaded_variants:
        st.info("Load spectra first to view results")
        return

    # Filter to variants with spectra
    variants_with_spectra = {
        k: v for k, v in loaded_variants.items()
        if v.query_spectra or v.reference_spectra
    }

    if not variants_with_spectra:
        st.warning("No spectra loaded. Check that mzML files match the idXML identifications.")
        return

    # Create tabs for each variant
    tab_names = list(variants_with_spectra.keys())
    if len(tab_names) > 10:
        tab_names = tab_names[:10]
        st.info(f"Showing first 10 of {len(variants_with_spectra)} variants")

    tabs = st.tabs([variants_with_spectra[k].display_label for k in tab_names])

    for tab, variant_key in zip(tabs, tab_names):
        with tab:
            render_variant_plots(variants_with_spectra[variant_key])

    # Export section
    st.markdown("---")
    render_export_section(loaded_variants)


def render_variant_plots(vg):
    """Render mirror plots for a variant group."""
    from core.peak_annotation import prepare_spectrum_for_plot
    from core.mirror_plot import create_mirror_plot, figure_to_bytes

    ion_types = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])
    mz_tolerance = st.session_state.get(get_state_key('mz_tolerance'), 0.02)

    query_spectra = vg.query_spectra
    ref_spectra = vg.reference_spectra

    st.caption(f"{len(query_spectra)} query, {len(ref_spectra)} reference spectra")

    # Determine pairs to plot
    pairs = []
    if query_spectra and ref_spectra:
        for q in query_spectra:
            for r in ref_spectra:
                pairs.append((q, r, "Query vs Reference"))
    elif len(query_spectra) >= 2:
        first = query_spectra[0]
        for other in query_spectra[1:]:
            pairs.append((other, first, "Query vs Query"))

    if not pairs:
        st.info("Need at least one query and one reference spectrum, or two query spectra")
        return

    # Limit display
    max_plots = 6
    if len(pairs) > max_plots:
        st.info(f"Showing first {max_plots} of {len(pairs)} comparisons")
        pairs = pairs[:max_plots]

    # Display plots in columns
    cols = st.columns(min(2, len(pairs)))

    for i, (upper, lower, comparison_type) in enumerate(pairs):
        with cols[i % 2]:
            try:
                # Prepare data
                upper_df = prepare_spectrum_for_plot(
                    upper.mz, upper.intensity,
                    upper.peptidoform, upper.precursor_charge,
                    ion_types=ion_types, mz_tolerance=mz_tolerance
                )
                lower_df = prepare_spectrum_for_plot(
                    lower.mz, lower.intensity,
                    lower.peptidoform, lower.precursor_charge,
                    ion_types=ion_types, mz_tolerance=mz_tolerance
                )

                # Create plot
                title = f"{vg.base_sequence}"
                fig = create_mirror_plot(
                    upper_df, lower_df,
                    title=title,
                    upper_label=Path(upper.source_file).stem[:20],
                    lower_label=Path(lower.source_file).stem[:20],
                    figsize=(8, 6),
                    dpi=100
                )

                # Display
                st.pyplot(fig)

                # Download button
                png_bytes = figure_to_bytes(fig, 'png')
                st.download_button(
                    "📥 PNG",
                    data=png_bytes,
                    file_name=f"{vg.base_sequence}_{i+1}.png",
                    mime="image/png",
                    key=f"download_{vg.base_sequence}_{i}"
                )

                # Cleanup
                import matplotlib.pyplot as plt
                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating plot: {e}")


def render_export_section(loaded_variants: dict):
    """Render batch export controls."""
    st.subheader("Batch Export")

    col1, col2 = st.columns([2, 1])

    with col1:
        export_dir = st.text_input(
            "Export directory",
            value=str(Path.home() / "MSpecViz_exports"),
            key=f"{STATE_PREFIX}_export_dir"
        )

    with col2:
        formats = st.multiselect(
            "Formats",
            options=['png', 'svg'],
            default=['png'],
            key=f"{STATE_PREFIX}_export_formats"
        )

    if st.button("📤 Export All Plots", key=f"{STATE_PREFIX}_export_all_btn", type="primary"):
        run_batch_export(loaded_variants, export_dir, formats)


def run_batch_export(loaded_variants: dict, export_dir: str, formats: list):
    """Run batch export of all plots."""
    from core.export import export_mirror_plots

    ion_types = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])
    mz_tolerance = st.session_state.get(get_state_key('mz_tolerance'), 0.02)

    progress_bar = st.progress(0, "Starting export...")

    def progress_callback(current, total, message):
        if total > 0:
            progress_bar.progress(current / total, message)

    try:
        exported = export_mirror_plots(
            loaded_variants,
            export_dir,
            formats=formats,
            ion_types=ion_types,
            mz_tolerance=mz_tolerance,
            progress_callback=progress_callback
        )

        progress_bar.progress(1.0, "Export complete!")

        total_files = sum(len(files) for files in exported.values())
        st.success(f"Exported {total_files} files to {export_dir}")

        # Show export location
        for fmt, files in exported.items():
            if files:
                st.caption(f"{fmt.upper()}: {len(files)} files")

    except Exception as e:
        st.error(f"Export failed: {e}")
        logger.exception("Export error")


def render():
    """Main render function for the workspace."""
    init_session_state()

    st.title("Peptide Mirror Analysis")
    st.caption("Compare experimental spectra against reference spectra for target peptides")

    # Workflow steps as tabs
    steps = [
        ("1. Files", render_file_selection),
        ("2. Peptides", render_peptide_input),
        ("3. Reference", render_reference_selection),
        ("4. Settings", render_settings),
        ("5. Search", render_search),
        ("6. Analysis", render_analysis),
        ("7. Results", render_results_gallery),
    ]

    # Create tabs for steps
    tab_names = [s[0] for s in steps]
    tabs = st.tabs(tab_names)

    for tab, (name, render_func) in zip(tabs, steps):
        with tab:
            render_func()
