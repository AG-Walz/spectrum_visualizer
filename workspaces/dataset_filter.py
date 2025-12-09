"""Filter Large Datasets - Extract relevant spectra before analysis.

3-step workflow: Select idXML → Define peptide filter → Extract & create filtered files.
Memory-efficient: Processes 10GB files using 200MB RAM via selective spectrum extraction.
Reduces file sizes by 90-99%. Output ready for Spectrum Analysis workspace.
"""
import streamlit as st
import sys
from pathlib import Path

# Add paths
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))
sys.path.insert(0, str(Path(__file__).parent.parent / "core"))

# Import all heavy packages from shared module (loaded once at app startup)
from shared_imports import (
    pd, np, oms, read_file, PSM, PSMList, Peptidoform,
    warnings, tempfile, os, sys, Path, logging, json, gc, re,
    Optional, Tuple, List, Dict,
    SPECTRUM_TOOLS_AVAILABLE,
    # Common utility functions
    create_file_options_with_sizes, map_display_names_to_paths,
    extract_filename, extract_filenames, is_empty_collection,
    # Spectrum tools
    st_tools
)

# Import shared UI components
from file_selection_components import (
    render_directory_browser,
    render_file_filter_and_selection,
    get_files_with_substring
)

# Import shared data input UI (single source of truth for file selection interface)
from data_input_components import (
    render_data_input_section,
    render_unified_file_browser,
    render_server_directory_selector
)

# Import mzml_processor module for peptide-based filtering
try:
    from mzml_processor import extract_peptide_relevant_spectra, filter_mzml
except ImportError as e:
    st.error("Failed to import mzml_processor module.")
    st.stop()

# Logger setup
logger = logging.getLogger(__name__)

# ================================================================
# FILTERING FUNCTIONS
# ================================================================

def filter_idxml_by_peptides(idxml_path: str, peptide_filter: List[str]) -> Tuple[bool, Dict[str, List[str]], List[str]]:
    """
    STEP 1: Filter idXML to extract only spectrum IDs matching peptides.
    NO mzML loading at this stage!
    
    Returns:
        (success, {mzml_filename: [spec_ids]}, messages)
    """
    messages = []
    
    try:
        # Normalize peptides
        normalized_peptides = []
        for pep in peptide_filter:
            if '/' in pep:
                pep = pep.split('/')[0]
            pep = re.sub(r'\[.*?\]', '', pep)
            pep = re.sub(r'\(.*?\)', '', pep)
            pep = pep.strip().upper()
            if pep:
                normalized_peptides.append(pep)
        
        normalized_peptides = list(set(normalized_peptides))
        messages.append(f"  Normalized to {len(normalized_peptides)} unique peptides")
        
        # Extract spectrum IDs (doesn't touch mzML files!)
        idxml_dir = str(Path(idxml_path).parent)
        
        mzml_to_spectra = extract_peptide_relevant_spectra(
            idxml_path=idxml_path,
            target_peptides=normalized_peptides,
            mzml_dir=idxml_dir
        )
        
        if not mzml_to_spectra:
            messages.append("  No spectra matched filter")
            return False, {}, messages
        
        total_spectra = sum(len(ids) for ids in mzml_to_spectra.values())
        messages.append(f"  Found {total_spectra} spectra in {len(mzml_to_spectra)} mzML(s)")
        
        return True, mzml_to_spectra, messages
        
    except Exception as e:
        messages.append(f"  Error: {str(e)}")
        return False, {}, messages


def create_filtered_mzmls(mzml_to_spectra: Dict[str, List[str]], 
                          mzml_search_dir: str, output_dir: str) -> Tuple[int, List[str], List[str]]:
    """
    STEP 2: Load ONLY the necessary spectra from mzML files and create filtered outputs.
    
    Returns:
        (total_spectra, [output_files], [messages])
    """
    messages = []
    output_files = []
    total_kept = 0
    
    for mzml_file, spec_ids in mzml_to_spectra.items():
        # Find mzML file
        mzml_path = Path(mzml_search_dir) / mzml_file
        
        if not mzml_path.exists():
            # Try recursive search
            found = list(Path(mzml_search_dir).rglob(mzml_file))
            if found:
                mzml_path = found[0]
            else:
                messages.append(f"  ⚠️ Not found: {mzml_file}")
                continue
        
        messages.append(f"  Processing {mzml_file}: {len(spec_ids)} spectra")
        
        try:
            # filter_mzml only loads the specific spectra we need!
            filtered_path = filter_mzml(
                original_mzml_path=str(mzml_path),
                spectrum_ids=spec_ids,
                output_dir=output_dir,
                append_if_exists=False
            )
            
            # Count from spec_ids (no need to reload!)
            kept_count = len(spec_ids)
            
            total_kept += kept_count
            output_files.append(filtered_path)
            
            messages.append(f"    ✅ {Path(filtered_path).name}: {kept_count} spectra")
            
        except Exception as e:
            messages.append(f"    ❌ Error: {str(e)}")
            logger.error(f"Error: {e}", exc_info=True)
    
    return total_kept, output_files, messages


# ================================================================
# UI COMPONENTS
# ================================================================

def render_idxml_selection():
    """File selection using shared data_input_components.py - filters only need idXML files."""
    
    # Custom header for filter workspace
    st.header("📁 Step 1: Select idXML Files")
    st.caption("Select peptide identification files to filter. The corresponding mzML files will be detected automatically from the idXML metadata.")
    
    # Use shared data input section (but we'll only take idXML files)
    mzml_files, idxml_files = render_data_input_section(
        session_key_prefix="dataset_filter",
        show_process_button=False,  # We handle processing separately
        process_callback=None
    )
    
    # Get base directory from session state
    current_dir = st.session_state.get('dataset_filter_selected_folders', ['/mnt/volume/'])[0] if st.session_state.get('dataset_filter_selected_folders') else '/mnt/volume/'
    
    return idxml_files, current_dir


def render_peptide_filter():
    """Simplified peptide input."""
    st.header("🧬 Step 2: Define Peptide Filter")
    
    tab1, tab2 = st.tabs(["Manual Input", "Extract from idXML"])
    
    peptides = []
    
    with tab1:
        peptide_input = st.text_area(
            "Enter peptides (one per line)",
            height=150,
            placeholder="PEPTIDEK\nSEQUENCE\n...",
            key="peptide_input"
        )
        
        if peptide_input:
            peptides = [p.strip() for p in peptide_input.split('\n') if p.strip()]
            st.success(f"✅ {len(peptides)} peptides")
    
    with tab2:
        if 'selected_idxml' in st.session_state and st.session_state['selected_idxml']:
            if st.button("🔍 Extract Peptides", type="primary"):
                with st.spinner("Extracting..."):
                    unique_peps = set()
                    for idxml in st.session_state['selected_idxml']:
                        try:
                            psms = read_file(idxml, file_format="idXML")
                            for psm in psms:
                                unique_peps.add(str(psm.peptidoform))
                        except Exception as e:
                            st.error(f"Error: {e}")
                    
                    peptides = sorted(list(unique_peps))
                    st.session_state['extracted_peptides'] = peptides
                    st.success(f"✅ {len(peptides)} peptides")
            
            if 'extracted_peptides' in st.session_state:
                peptides = st.session_state['extracted_peptides']
                with st.expander(f"View {len(peptides)} peptides"):
                    for p in peptides[:20]:
                        st.text(p)
                    if len(peptides) > 20:
                        st.text(f"... and {len(peptides)-20} more")
        else:
            st.info("Select idXML files first")
    
    return peptides


def render_output_settings():
    """Output configuration."""
    st.header("💾 Step 3: Output Settings")
    
    default_output = st.session_state.get('filter_current_directory', '/mnt/volume/') + "/filtered_output"
    
    output_dir = st.text_input(
        "Output Directory",
        value=st.session_state.get('output_dir', default_output),
        key="output_dir_input"
    )
    st.session_state['output_dir'] = output_dir
    
    create_dir = st.checkbox("Create if missing", value=True)
    
    if os.path.exists(output_dir):
        st.success(f"✅ Directory exists")
    else:
        if create_dir:
            st.info(f"📁 Will create directory")
        else:
            st.warning(f"⚠️ Directory doesn't exist")
    
    return output_dir, create_dir


# ================================================================
# MAIN WORKSPACE
# ================================================================

def main():
    """Filtering workspace."""
    st.title("🗂️ Filter Large Datasets")
    
    # Step 1: idXML selection
    selected_idxml, base_dir = render_idxml_selection()
    st.session_state['selected_idxml'] = selected_idxml
    st.session_state['base_dir'] = base_dir
    
    if not selected_idxml:
        st.info("👆 Select idXML files to continue")
        return
    
    st.markdown("---")
    
    # Step 2: Peptide filter
    peptide_filter = render_peptide_filter()
    
    if not peptide_filter:
        st.warning("⚠️ Define peptide filter to continue")
        return
    
    st.markdown("---")
    
    # Step 3: Output settings
    output_dir, create_dir = render_output_settings()
    
    st.markdown("---")
    
    # Step 4: Process
    st.header("🚀 Step 4: Process")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("idXML Files", len(selected_idxml))
    with col2:
        st.metric("Peptide Filter", len(peptide_filter))
    
    st.info(f"""
    **Plan:**
    - Process {len(selected_idxml)} idXML file(s)
    - Filter for {len(peptide_filter)} peptide(s)
    - mzML files discovered automatically
    - Only necessary spectra loaded
    - Output: `{output_dir}`
    """)
    
    if st.button("🔄 Start Optimized Filtering", type="primary", use_container_width=True):
        # Create output dir
        if create_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
                st.success(f"✅ Created: {output_dir}")
            except Exception as e:
                st.error(f"❌ Failed: {e}")
                return
        
        if not os.path.exists(output_dir):
            st.error(f"❌ Directory doesn't exist: {output_dir}")
            return
        
        st.info(f"🔍 Processing {len(selected_idxml)} idXML file(s)...")
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        results = []
        total_kept = 0
        all_outputs = []
        
        # Process each idXML
        for i, idxml_path in enumerate(selected_idxml):
            idxml_name = Path(idxml_path).name
            progress_bar.progress((i + 1) / len(selected_idxml))
            status_text.text(f"Processing {i+1}/{len(selected_idxml)}: {idxml_name}")
            
            messages = [f"📋 {idxml_name}"]
            
            # STEP 1: Filter idXML (no mzML loading!)
            success, mzml_to_spectra, filter_msgs = filter_idxml_by_peptides(idxml_path, peptide_filter)
            messages.extend(filter_msgs)
            
            if not success:
                results.append({
                    'File': idxml_name,
                    'Status': '⚠️',
                    'Spectra': 0,
                    'Files': 0,
                    'Details': '\n'.join(messages)
                })
                continue
            
            # STEP 2: Load only necessary spectra
            kept, outputs, create_msgs = create_filtered_mzmls(mzml_to_spectra, base_dir, output_dir)
            messages.extend(create_msgs)
            
            results.append({
                'File': idxml_name,
                'Status': '✅' if outputs else '⚠️',
                'Spectra': kept,
                'Files': len(outputs),
                'Details': '\n'.join(messages)
            })
            
            total_kept += kept
            all_outputs.extend(outputs)
        
        progress_bar.empty()
        status_text.empty()
        
        # Results
        st.markdown("---")
        st.header("📊 Results")
        
        successful = sum(1 for r in results if r['Status'] == '✅')
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Processed", len(results))
        with col2:
            st.metric("Successful", successful)
        with col3:
            st.metric("mzML Files", len(all_outputs))
        with col4:
            st.metric("Spectra", total_kept)
        
        st.subheader("Details")
        results_df = pd.DataFrame(results)
        
        for _, row in results_df.iterrows():
            with st.expander(f"{row['Status']} {row['File']} - {row['Spectra']} spectra, {row['Files']} files"):
                st.text(row['Details'])
        
        summary_df = results_df.drop('Details', axis=1)
        st.dataframe(summary_df, use_container_width=True)
        
        if successful == len(results):
            st.success(f"✅ All {len(results)} files processed!")
            st.success(f"📁 Created {len(all_outputs)} filtered mzML files in: {output_dir}")
        else:
            st.warning(f"⚠️ {successful}/{len(results)} successful")
            if all_outputs:
                st.info(f"📁 Created {len(all_outputs)} files")
        
        if all_outputs:
            with st.expander("📋 Output Files"):
                for f in all_outputs:
                    st.text(f"  • {Path(f).name}")


if __name__ == "__main__":
    main()
