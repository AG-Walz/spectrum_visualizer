import streamlit as st
import pandas as pd
from pathlib import Path
import tempfile
import os
import sys
from typing import Optional, Tuple, List, Dict
import logging
import json
import gc
import pyopenms as oms
from psm_utils import Peptidoform, PSM, PSMList

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure environment
os.environ.setdefault("BOKEH_ALLOW_WS_ORIGIN", "*")

def page_setup():
    """Common page setup function."""
    # Load settings
    if "settings" not in st.session_state:
        try:
            with open("settings.json", "r") as f:
                st.session_state.settings = json.load(f)
        except FileNotFoundError:
            st.session_state.settings = {
                "app-name": "immunuXplore",
                "version": "0.0.1"
            }

page_setup()

# Import spectrum tools
bin_dir = Path(__file__).parent.parent / "bin"
if str(bin_dir) not in sys.path:
    sys.path.insert(0, str(bin_dir))

try:
    import spectrum_tools as st_tools
except ImportError:
    try:
        from bin import spectrum_tools as st_tools
    except ImportError:
        st.error("❌ Could not import spectrum_tools module. Please check installation.")
        st.stop()

# ================================================================
# CORE DATA FUNCTIONS
# ================================================================

def clear_memory_cache():
    """Clear memory and caches to prevent OOM errors."""
    try:
        # Clear Streamlit cache
        st.cache_data.clear()
        
        # Force garbage collection
        gc.collect()
        
        # Remove the problematic sys.intern.clear() - sys.intern is a function, not an object
        # sys.intern.clear()  # This line causes the error - commented out
        
        logger.info("Memory cache cleared")
    except Exception as e:
        logger.warning(f"Error clearing memory cache: {e}")

@st.cache_data
def save_uploaded_file(uploaded_file, suffix: Optional[str] = None) -> str:
    """Save uploaded file to temporary directory and return path with original filename."""
    temp_dir = tempfile.mkdtemp()
    
    # Use original filename if available, otherwise fall back to generic name
    if hasattr(uploaded_file, 'name') and uploaded_file.name:
        original_name = uploaded_file.name
        # Ensure the suffix matches if provided
        if suffix and not original_name.endswith(suffix):
            file_name = f"{os.path.splitext(original_name)[0]}{suffix}"
        else:
            file_name = original_name
    else:
        # Fallback for files without names
        file_name = f"uploaded_file{suffix or ''}"
    
    file_path = os.path.join(temp_dir, file_name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getvalue())
    return file_path

# Remove @st.cache_data decorator to prevent memory accumulation
# @st.cache_data
def load_and_validate_files_multi(mzml_paths: List[str], idxml_path: str) -> Tuple[bool, str, Optional[pd.DataFrame], Optional[Dict]]:
    """Load and validate multiple mzML files with idXML. Returns (success, message, mapping_df, validation_info)."""
    try:
        logger.info(f"Validating files: {len(mzml_paths)} mzML files, {idxml_path}")
        
        # Extract expected mzML files from idXML
        expected_files = st_tools.extract_mzml_files_from_idxml(idxml_path)
        uploaded_files = [Path(p).name for p in mzml_paths]
        
        # Validate file mapping
        matched_files, missing_files, extra_files = st_tools.validate_mzml_files_mapping(uploaded_files, expected_files)
        
        validation_info = {
            "expected_files": expected_files,
            "uploaded_files": uploaded_files,
            "matched_files": matched_files,
            "missing_files": missing_files,
            "extra_files": extra_files
        }
        
        # Test loading the identification file
        protein_ids, peptide_ids = st_tools.load_identifications(idxml_path)
        
        if len(peptide_ids) == 0:
            return False, "ID file contains no peptide identifications", None, validation_info
        
        # Memory-efficient validation: only test one file at a time
        total_spectra = 0
        valid_files = 0
        
        for i, mzml_path in enumerate(mzml_paths):
            try:
                logger.info(f"Validating file {i+1}/{len(mzml_paths)}: {Path(mzml_path).name}")
                exp = st_tools.load_experiment(mzml_path)
                file_spectra = exp.size()
                total_spectra += file_spectra
                valid_files += 1
                logger.info(f"File {Path(mzml_path).name}: {file_spectra} spectra")
                
                # Clear memory immediately after validation
                del exp
                gc.collect()
                
            except Exception as e:
                logger.warning(f"Could not load {mzml_path}: {e}")
        
        if valid_files == 0:
            return False, "No mzML files could be loaded", None, validation_info
        
        if total_spectra == 0:
            return False, "No mzML files contain spectra", None, validation_info
        
        logger.info(f"Validation complete: {valid_files}/{len(mzml_paths)} files valid, {total_spectra} total spectra")
        
        # Extract mapping using memory-efficient multi-file function
        mapping_df = st_tools.extract_spectrum_mapping_multi(mzml_paths, idxml_path)
        
        if mapping_df.empty:
            return False, "No matching spectra found between mzML and ID files", None, validation_info
        
        success_msg = f"Successfully loaded {len(mapping_df)} spectrum-peptide mappings from {len(mzml_paths)} files"
        logger.info(success_msg)
        
        # Force memory cleanup before returning
        clear_memory_cache()
        
        return True, success_msg, st_tools.ensure_streamlit_compatible_dataframe(mapping_df), validation_info
        
    except Exception as e:
        logger.error(f"File validation error: {e}")
        return False, f"Error loading files: {str(e)}", None, None

@st.cache_data
def load_and_validate_files(mzml_path: str, idxml_path: str) -> Tuple[bool, str, Optional[pd.DataFrame]]:
    """Load and validate single mzML file with idXML. Returns (success, message, mapping_df)."""
    success, message, mapping_df, _ = load_and_validate_files_multi([mzml_path], idxml_path)
    return success, message, mapping_df

def get_available_models() -> Dict[str, str]:
    """Get available Koina models from settings or spectrum_tools."""
    try:
        # Try to get from settings first
        settings = st.session_state.get("settings", {})
        if "models" in settings:
            available_models = settings["models"]["available"]
            default_model = settings["models"]["default"]
            
            # Build display dictionary with default first
            ordered_models = {}
            
            # Add default model first
            if default_model in available_models:
                ordered_models[f"{default_model} (Default)"] = default_model
            
            # Add remaining models
            for model in available_models:
                if model != default_model:
                    ordered_models[model] = model
            
            return ordered_models
        
        # Fallback to spectrum_tools
        basic_models = st_tools.get_available_models()
        
        # Ensure default is first
        ordered_models = {}
        default_model = settings.get("models", {}).get("default", "ms2pip_timsTOF2024")
        
        if default_model in basic_models:
            display_name = f"{default_model} - {basic_models[default_model]}"
            ordered_models[display_name] = default_model
        
        # Add remaining models
        for model_name, description in basic_models.items():
            if model_name != default_model:
                display_name = f"{model_name} - {description}"
                ordered_models[display_name] = model_name
        
        return ordered_models
        
    except Exception as e:
        logger.warning(f"Could not get models: {e}")
        return {
            "ms2pip_timsTOF2024 (Default)": "ms2pip_timsTOF2024",
            "ms2pip_HCD2021": "ms2pip_HCD2021",
            "Prosit_2019_intensity": "Prosit_2019_intensity",
        }

# ================================================================
# UI COMPONENTS
# ================================================================

def render_data_input_section():
    """Render the data input section."""
    st.header("📁 Data Input")
    
    # File upload section
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📊 Mass Spectrometry Data")
        
        # Get max file size from settings
        max_size_gb = st.session_state.settings.get("performance", {}).get("max_file_size_gb", 2)
        
        # Multi-file uploader for mzML files
        mzml_files = st.file_uploader(
            "Select mzML files",
            type=["mzml", "mzML"],
            accept_multiple_files=True,
            help=f"Mass spectrometry data files containing MS/MS spectra (up to {max_size_gb} GB each)"
        )
        
        if mzml_files:
            total_size_mb = sum(f.size for f in mzml_files) / 1024 / 1024
            st.success(f"✅ Loaded {len(mzml_files)} mzML files")
            
            # Show file details
            with st.expander("📋 File Details"):
                for i, file in enumerate(mzml_files):
                    file_size_mb = file.size / 1024 / 1024
                    st.write(f"**{i+1}.** {file.name} ({file_size_mb:.1f} MB)")
            
            st.info(f"Total size: {total_size_mb:.1f} MB")
            
            # Show warning if files are approaching the limit
            if total_size_mb > (max_size_gb * 1024 * 0.8):  # 80% of limit
                st.warning(f"⚠️ Large files detected. Maximum supported size per file is {max_size_gb} GB.")
    
    with col2:
        st.subheader("🔍 Identification Data")
        idxml_file = st.file_uploader(
            "Select identification file", 
            type=["idxml", "idXML", "mzid", "mzidentml"],
            help=f"Peptide identification file from database search (up to {max_size_gb} GB)"
        )
        
        if idxml_file:
            file_size_mb = idxml_file.size / 1024 / 1024
            st.success(f"✅ Loaded: {idxml_file.name}")
            st.info(f"File size: {file_size_mb:.1f} MB")
    
    # Process files if both are uploaded
    if mzml_files and idxml_file is not None:
        # Process files button
        if st.button("🔄 Process Files", type="primary"):
            with st.spinner("🔄 Processing uploaded files..."):
                # Save all mzML files with original names preserved
                mzml_paths = []
                for mzml_file in mzml_files:
                    mzml_path = save_uploaded_file(mzml_file)  # No suffix - preserve original name
                    mzml_paths.append(mzml_path)
                
                # Save idXML file with original name preserved
                idxml_path = save_uploaded_file(idxml_file)  # No suffix - preserve original name
            
            # Validate files
            success, message, mapping_df, validation_info = load_and_validate_files_multi(mzml_paths, idxml_path)
            
            if success and mapping_df is not None:
                st.success(f"✅ {message}")
                
                # Store in session state
                st.session_state.update({
                    'files_loaded': True,
                    'mzml_paths': mzml_paths,
                    'idxml_path': idxml_path,
                    'mapping_df': mapping_df,
                    'validation_info': validation_info
                })
                
                # Show file mapping validation
                if validation_info:
                    st.markdown("---")
                    render_file_mapping_validation(validation_info)
                
                st.rerun()
            else:
                st.error(f"❌ {message}")
                
                # Show validation information even if loading failed
                if validation_info:
                    st.markdown("---")
                    render_file_mapping_validation(validation_info)
                
                with st.expander("💡 Troubleshooting"):
                    st.markdown("""
                    **Common Issues:**
                    - Ensure mzML files contain MS2 spectra
                    - Check ID file has peptide identifications
                    - Verify file names match exactly (case-sensitive)
                    - Confirm all files are from same experiment
                    - Check file formats are supported
                    """)
    
    return mzml_files, idxml_file

def render_analysis_settings():
    """Render analysis settings section."""
    st.header("🤖 Analysis Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("AI Model Selection")
        available_models = get_available_models()
        model_names = list(available_models.keys())
        
        selected_model_display = st.selectbox(
            "Choose prediction model",
            options=model_names,
            index=0,
            help="AI model for spectrum prediction via Koina API"
        )
        
        selected_model = available_models[selected_model_display]
        st.info(f"Selected: {selected_model}")
    
    with col2:
        st.subheader("Advanced Options")
        
        debug_mode = st.checkbox(
            "Enable Debug Mode",
            value=st.session_state.settings.get("ui", {}).get("debug_mode", False),
            help="Show detailed processing information"
        )
        
        resource_monitor = st.checkbox(
            "Resource Monitoring",
            value=st.session_state.settings.get("performance", {}).get("enable_monitoring", True),
            help="Display CPU and memory usage"
        )
        
        if resource_monitor:
            st.info("💡 Resource monitoring enabled - check Settings to configure")
            try:
                import psutil
                cpu_percent = psutil.cpu_percent(interval=1)
                memory_percent = psutil.virtual_memory().percent
                
                col2a, col2b = st.columns(2)
                with col2a:
                    st.metric("CPU Usage", f"{cpu_percent:.1f}%")
                with col2b:
                    st.metric("Memory Usage", f"{memory_percent:.1f}%")
            except ImportError:
                st.info("Install psutil for resource monitoring: `pip install psutil`")
    
    return selected_model, debug_mode

def render_spectrum_selector(mapping_df: pd.DataFrame) -> Optional[int]:
    """Render spectrum selection interface."""
    st.header("🔍 Spectrum Selection")
    
    # Summary statistics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Spectra", len(mapping_df))
    
    with col2:
        unique_peptides = mapping_df["peptide_sequence"].nunique() if "peptide_sequence" in mapping_df.columns else 0
        st.metric("Unique Peptides", unique_peptides)
    
    with col3:
        if "charge" in mapping_df.columns:
            charge_range = f"{mapping_df['charge'].min()}+ to {mapping_df['charge'].max()}+"
            st.metric("Charge Range", charge_range)
        else:
            st.metric("Charge Range", "N/A")
    
    with col4:
        if "source_file" in mapping_df.columns:
            unique_files = mapping_df["source_file"].nunique()
            st.metric("Source Files", unique_files)
        elif "score" in mapping_df.columns:
            avg_score = mapping_df["score"].mean()
            st.metric("Avg Score", f"{avg_score:.3f}")
        else:
            st.metric("Avg Score", "N/A")
    
    # Data table
    st.subheader("📋 Available Spectra")
    
    # Create display dataframe
    display_df = mapping_df.copy()
    display_columns = []
    
    # Add spectrum index
    if "spectrum_index" in display_df.columns:
        display_df["Index"] = display_df["spectrum_index"].astype(int)
        display_columns.append("Index")
    
    # Add peptide sequence with proper modification formatting
    if "peptide_sequence" in display_df.columns:
        display_df["Peptide"] = display_df["peptide_sequence"].apply(
            lambda x: st_tools.convert_to_mass_shift_format(str(x))
        )
        display_columns.append("Peptide")
    
    # Add additional columns
    for col, display_name in [
        ("charge", "Charge"),
        ("rt", "RT (min)"),
        ("rt_spec", "RT (min)"),
        ("precursor_mz", "Precursor m/z"),
        ("precursor_mz_spec", "Precursor m/z"),
        ("score", "Score"),
        ("source_file", "Source File")
    ]:
        if col in display_df.columns and display_name not in display_columns:
            if "RT" in display_name:
                display_df[display_name] = (display_df[col] / 60).round(2)
            elif "m/z" in display_name:
                display_df[display_name] = display_df[col].round(4)
            elif col == "score":
                display_df[display_name] = display_df[col].round(3)
            else:
                display_df[display_name] = display_df[col]
            display_columns.append(display_name)
    
    # Show data table
    st.dataframe(
        display_df[display_columns] if display_columns else display_df,
        use_container_width=True,
        height=400
    )
    
    # Spectrum selection
    available_indices = sorted(display_df["spectrum_index"].unique()) if "spectrum_index" in display_df.columns else list(range(len(display_df)))
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        selected_idx = st.selectbox(
            "Select spectrum index for visualization",
            options=available_indices,
            index=0,
            help=f"Choose from {len(available_indices)} spectra"
        )
        spectrum_index = int(selected_idx)
    
    with col2:
        if st.button("🔄 Refresh Data", help="Reload spectrum data"):
            st.rerun()
    
    return spectrum_index

def render_spectrum_info(mapping_df: pd.DataFrame, spectrum_index: int):
    """Render spectrum information section."""
    st.header("📊 Spectrum Information")
    
    # Get spectrum data
    if "spectrum_index" in mapping_df.columns:
        matching_rows = mapping_df[mapping_df["spectrum_index"] == spectrum_index]
        if matching_rows.empty:
            st.error(f"Spectrum index {spectrum_index} not found")
            return
        row = matching_rows.iloc[0]
    else:
        if spectrum_index >= len(mapping_df):
            st.error(f"Invalid spectrum index: {spectrum_index}")
            return
        row = mapping_df.iloc[spectrum_index]
    
    # Display key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        peptide_seq = st_tools.convert_to_mass_shift_format(str(row.get("peptide_sequence", "N/A")))
        st.metric("Peptide", peptide_seq)
    
    with col2:
        charge = int(row.get("charge", 0))
        st.metric("Charge", f"{charge}+")
    
    with col3:
        rt_value = row.get("rt", row.get("rt_spec", 0))
        rt = rt_value / 60
        st.metric("RT", f"{rt:.2f} min")
    
    with col4:
        if "source_file" in row.index:
            source_file = row["source_file"]
            st.metric("Source File", source_file)
        else:
            precursor_mz = row.get("precursor_mz", row.get("precursor_mz_spec", 0))
            st.metric("Precursor m/z", f"{precursor_mz:.4f}")
    
    # Additional info row if we have source file
    if "source_file" in row.index:
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            precursor_mz = row.get("precursor_mz", row.get("precursor_mz_spec", 0))
            st.metric("Precursor m/z", f"{precursor_mz:.4f}")
        with col2:
            if "score" in row.index:
                st.metric("Score", f"{row['score']:.3f}")
        with col3:
            pass
        with col4:
            pass
    
    # Detailed information
    with st.expander("📋 Detailed Properties"):
        details = {}
        for col in mapping_df.columns:
            if col in row.index:
                details[col.replace("_", " ").title()] = row[col]
        
        detail_df = pd.DataFrame(list(details.items()), columns=["Property", "Value"])
        st.dataframe(detail_df, use_container_width=True, hide_index=True)

def render_visualization_section(
    mzml_path: str,
    idxml_path: str, 
    spectrum_index: int,
    model_name: str,
    mapping_df: pd.DataFrame,
    debug_mode: bool = False
):
    """Render visualization section."""
    st.header("📈 Spectrum Visualization")
    
    # Get spectrum data
    if "spectrum_index" in mapping_df.columns:
        matching_rows = mapping_df[mapping_df["spectrum_index"] == spectrum_index]
        if matching_rows.empty:
            st.error(f"Spectrum index {spectrum_index} not found")
            return
        row = matching_rows.iloc[0]
    else:
        if spectrum_index >= len(mapping_df):
            st.error(f"Invalid spectrum index: {spectrum_index}")
            return
        row = mapping_df.iloc[spectrum_index]
    
    # Generate visualization
    with st.spinner("🔄 Generating spectrum visualization..."):
        try:
            # Debug information
            if debug_mode:
                with st.expander("🐛 Debug Information"):
                    st.json({
                        "spectrum_index": spectrum_index,
                        "model_name": model_name,
                        "peptide_sequence": row.get("peptide_sequence", "Unknown"),
                        "charge": row.get("charge", "Unknown"),
                        "files": {"mzml": mzml_path, "idxml": idxml_path}
                    })
            
            # Generate plot
            plot = st_tools.visualize_spectrum(
                mzml_path=mzml_path,
                id_path=idxml_path,
                spectrum_index=spectrum_index,
                model_name=model_name,
                width=1200,
                height=700,
                debug=debug_mode
            )
            
            if plot is None:
                st.error("❌ Failed to generate visualization")
                return
            
            # Display plot
            plot_displayed = display_plot_with_fallbacks(plot, debug_mode)
            
            if plot_displayed:
                st.success("✅ Visualization generated successfully!")
                
                # Show plot information
                with st.expander("ℹ️ Visualization Guide"):
                    st.markdown(f"""
                    **Experimental Spectrum** (top): Real peaks from mass spectrometer
                    **Predicted Spectrum** (bottom): AI-generated prediction using {model_name}
                    
                    **Ion Color Coding:**
                    - 🔵 **Blue**: b-ions (N-terminal fragments)
                    - 🟠 **Orange**: y-ions (C-terminal fragments)
                    - 🟢 **Green**: Other fragment types
                    
                    **Intensity Matching:**
                    - Darker colors: Higher confidence/intensity
                    - Peak annotations show ion type and mass
                    
                    **Peptide**: {row.get('peptide_sequence', 'Unknown')}
                    **Model**: {model_name}
                    """)
            else:
                st.error("❌ Visualization failed - check debug information")
                
        except Exception as e:
            logger.error(f"Visualization error: {e}")
            st.error(f"❌ Error: {str(e)}")
            
            if debug_mode:
                st.exception(e)

def display_plot_with_fallbacks(plot, debug_mode: bool) -> bool:
    """Display plot with multiple fallback methods."""
    plot_displayed = False
    
    # Method 1: streamlit-bokeh component
    try:
        from streamlit_bokeh import st_bokeh
        st_bokeh(plot, use_container_width=True)
        plot_displayed = True
        if debug_mode:
            st.success("✅ Using streamlit-bokeh component")
    except ImportError:
        if debug_mode:
            st.info("streamlit-bokeh not available, trying standard method")
    except Exception as e:
        if debug_mode:
            st.warning(f"streamlit-bokeh failed: {e}")
    
    # Method 2: Standard Streamlit
    if not plot_displayed:
        try:
            st.bokeh_chart(plot, use_container_width=True)
            plot_displayed = True
            if debug_mode:
                st.success("✅ Using standard bokeh chart")
        except Exception as e:
            if debug_mode:
                st.error(f"Standard bokeh failed: {e}")
    
    # Method 3: HTML export
    if not plot_displayed:
        try:
            from bokeh.embed import file_html
            from bokeh.resources import CDN
            html = file_html(plot, CDN, "Spectrum Visualization")
            st.components.v1.html(html, height=750, scrolling=True)
            plot_displayed = True
            if debug_mode:
                st.success("✅ Using HTML export")
        except Exception as e:
            if debug_mode:
                st.error(f"HTML export failed: {e}")
    
    return plot_displayed

def render_file_mapping_validation(validation_info: Dict):
    """Render file mapping validation information."""
    st.header("🔗 File Mapping Validation")
    
    if not validation_info:
        st.info("No validation information available")
        return
    
    expected_files = validation_info.get("expected_files", [])
    uploaded_files = validation_info.get("uploaded_files", [])
    matched_files = validation_info.get("matched_files", [])
    missing_files = validation_info.get("missing_files", [])
    extra_files = validation_info.get("extra_files", [])
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Expected Files", len(expected_files))
    
    with col2:
        st.metric("Uploaded Files", len(uploaded_files))
    
    with col3:
        st.metric("Matched Files", len(matched_files))
        
    with col4:
        if missing_files:
            st.metric("Missing Files", len(missing_files), delta=f"-{len(missing_files)}")
        else:
            st.metric("Missing Files", 0)
    
    # Detailed breakdown
    if matched_files:
        st.success(f"✅ **{len(matched_files)} files matched successfully:**")
        for file in matched_files:
            st.write(f"• {file}")
    
    if missing_files:
        st.error(f"❌ **{len(missing_files)} files are missing:**")
        for file in missing_files:
            st.write(f"• {file}")
        
        with st.expander("💡 Tips for Missing Files"):
            st.markdown("""
            **Common solutions:**
            - Check file names match exactly (case-sensitive)
            - Ensure all files from the same experiment are uploaded
            - Verify idXML file contains correct file references
            - Try uploading files in the correct order
            """)
    
    if extra_files:
        st.warning(f"⚠️ **{len(extra_files)} unexpected files uploaded:**")
        for file in extra_files:
            st.write(f"• {file}")
        
        st.info("These files may not be referenced in the identification file.")
    
    # Show expected files if validation fails
    if expected_files and (missing_files or not matched_files):
        with st.expander("📋 Expected Files from ID File"):
            st.markdown("**The identification file expects these mzML files:**")
            for file in expected_files:
                status = "✅" if file in matched_files else "❌"
                st.write(f"{status} {file}")

# ================================================================
# MAIN WORKSPACE
# ================================================================

def main():
    """Main workspace function - simplified based on SpecViz.ipynb workflow."""
    st.title("🔬 Spectrum Analysis Workspace")
    st.markdown("Upload your MS/MS data and analyze spectra with AI-powered predictions.")
    
    st.markdown("---")  # Separator
    
    # Step 1: File Upload (Keep existing functionality)
    render_data_input_section()
    
    # Only proceed if files are loaded
    if not st.session_state.get('files_loaded', False):
        st.info("👆 Please upload your mzML and idXML files to get started.")
        return
    
    # Step 2: Simplified Analysis Settings (based on notebook parameters)
    render_simplified_analysis_settings()
    
    # Step 3: Automatic Analysis and Visualization
    if st.session_state.get('analysis_ready', False):
        render_automated_analysis()

def render_simplified_analysis_settings():
    """Simplified analysis settings based on SpecViz.ipynb parameters."""
    st.header("⚙️ Analysis Settings")
    
    # Get available peptides from the PSM data
    try:
        mapping_df = st.session_state.get('mapping_df', pd.DataFrame())
        if mapping_df.empty:
            st.error("No PSM data available. Please check your files.")
            return
        
        # Extract unique unmodified peptide sequences (like in SpecViz.py line 84)
        # Load PSM data to extract proper peptide sequences
        idxml_path = st.session_state.get('idxml_path')
        if not idxml_path:
            st.error("No idXML file path found in session state.")
            return
        
        # Load PSM list and extract unmodified sequences
        from psm_utils.io import read_file
        from psm_utils import Peptidoform
        
        psm_list = read_file(idxml_path, file_format="idXML")
        
        # Extract unmodified sequences using the same method as SpecViz.py
        unmodified_sequences = set()
        for psm in psm_list:
            try:
                unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
                unmodified_sequences.add(unmodified_peptide)
            except Exception as e:
                logger.warning(f"Could not extract sequence from {psm.peptidoform}: {e}")
                continue
        
        unique_peptides = sorted(list(unmodified_sequences))
        
        if not unique_peptides:
            st.error("No valid peptide sequences found in the PSM data.")
            return
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Peptide selection
            selected_peptide = st.selectbox(
                "🧬 Select Peptide Sequence",
                unique_peptides,
                help="Choose the unmodified peptide sequence to analyze"
            )
            
            # Model selection
            available_models = st_tools.get_available_models()
            selected_model = st.selectbox(
                "🤖 Prediction Model",
                list(available_models.keys()),
                index=0,  # Default to first model
                help="AI model for fragment ion prediction"
            )
        
        with col2:
            # Maximum fragment charge
            max_frag_charge = st.number_input(
                "⚡ Max Fragment Charge",
                min_value=1,
                max_value=5,
                value=2,
                help="Maximum charge state for fragment ions"
            )
            
            # m/z tolerance
            mz_tolerance = st.number_input(
                "📏 m/z Tolerance (Da)",
                min_value=0.01,
                max_value=2.0,
                value=0.5,
                step=0.01,
                help="Mass tolerance for peak alignment"
            )
        
        # Show peptide statistics by counting PSMs with matching unmodified sequence
        matching_psm_count = 0
        for psm in psm_list:
            try:
                unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
                if unmodified_peptide == selected_peptide:
                    matching_psm_count += 1
            except:
                continue
        
        st.info(f"📊 Found **{matching_psm_count} PSMs** for peptide **{selected_peptide}**")
        
        # Start analysis button
        if st.button("🚀 Start Analysis", type="primary"):
            # Store settings in session state
            st.session_state.update({
                'selected_peptide': selected_peptide,
                'selected_model': available_models[selected_model],  # Get model ID
                'max_frag_charge': max_frag_charge,
                'mz_tolerance': mz_tolerance,
                'analysis_ready': True
            })
            st.rerun()
    
    except Exception as e:
        st.error(f"Error in analysis settings: {e}")
        logger.error(f"Analysis settings error: {e}", exc_info=True)

def render_automated_analysis():
    """Automated analysis following SpecViz.ipynb workflow."""
    st.header("📈 Analysis Results")
    
    try:
        # Get parameters from session state
        selected_peptide = st.session_state['selected_peptide']
        selected_model = st.session_state['selected_model']
        max_frag_charge = st.session_state['max_frag_charge']
        mz_tolerance = st.session_state['mz_tolerance']
        
        # Get file paths
        if 'mzml_paths' in st.session_state:
            mzml_paths = st.session_state['mzml_paths']
        else:
            mzml_paths = [st.session_state['mzml_path']]
        idxml_path = st.session_state['idxml_path']
        
        with st.spinner("🔄 Running analysis workflow..."):
            # Step 1: Load experiments (following notebook exactly)
            st.write("📁 Loading mzML experiments...")
            experiments = {}
            for mzml_file in mzml_paths:
                exp = oms.MSExperiment()
                oms.MzMLFile().load(mzml_file, exp)
                experiments[Path(mzml_file).name] = exp
            
            # Step 2: Create spectrum lookups (following notebook exactly)
            st.write("🔍 Creating spectrum lookups...")
            lookup_map = {}
            for exp_name, exp in experiments.items():
                spectrum_lookup = oms.SpectrumLookup()
                spectrum_lookup.readSpectra(exp, "scan=(?<SCAN>\\d+)")
                lookup_map[Path(exp_name).stem] = spectrum_lookup
            
            # Step 3: Load PSM list and filter by peptide (following notebook exactly)
            st.write("📋 Loading and filtering PSM data...")
            from psm_utils.io import read_file
            from psm_utils import Peptidoform
            psm_list = read_file(idxml_path, file_format="idXML")
            
            # Filter PSMs by peptide (exactly like notebook)
            filtered_psms = PSMList(psm_list=[])
            for psm in psm_list:
                unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
                if unmodified_peptide == selected_peptide:
                    filtered_psms.append(psm)
            
            if len(filtered_psms) == 0:
                st.error(f"No PSMs found for peptide: {selected_peptide}")
                return
            
            # Step 4: Attach spectra to PSMs (exactly like notebook)
            st.write("🔗 Attaching experimental spectra...")
            updated_psms = []
            for psm in filtered_psms:
                if not (hasattr(psm, "spectrum") and psm.spectrum is not None):
                    exp_name = f"{psm.run}.mzML"
                    if exp_name in experiments and psm.run in lookup_map:
                        try:
                            idx = lookup_map[psm.run].findByNativeID(psm.spectrum_id)
                            psm.spectrum = experiments[exp_name][idx]
                        except:
                            pass  # spectrum remains None
                updated_psms.append(psm)
            
            target_psms = PSMList(psm_list=updated_psms)
            
            # Filter PSMs that have spectra attached
            psms_with_spectra = []
            for psm in target_psms:
                if hasattr(psm, 'spectrum') and psm.spectrum is not None:
                    psms_with_spectra.append(psm)
            
            if len(psms_with_spectra) == 0:
                st.error("No experimental spectra could be loaded for the selected PSMs.")
                return
            
            st.success(f"✅ Loaded spectra for {len(psms_with_spectra)} out of {len(filtered_psms)} PSMs")
            
            # Step 5: Process alignments (exactly like notebook)
            st.write("⚖️ Processing spectrum alignments...")
            aligned_rows = []
            
            for psm in psms_with_spectra:
                try:
                    # Generate theoretical spectrum (exactly like notebook)
                    gen = oms.TheoreticalSpectrumGenerator()
                    params = gen.getParameters()
                    for ion in ("a","b","y","c","x","z"): 
                        params.setValue(f"add_{ion}_ions","true")
                    params.setValue("add_metainfo","true")
                    gen.setParameters(params)
                    theo_spec = oms.MSSpectrum()
                    gen.getSpectrum(theo_spec, oms.AASequence(str(psm.peptidoform).split('/')[0].replace('[','(').replace(']',')')), 
                                    1, max_frag_charge)
                    if theo_spec.size() > 1: 
                        theo_spec.sortByPosition()
                    
                    # Convert to dataframe (exactly like notebook)
                    if theo_spec is None or theo_spec.size() == 0: 
                        continue
                    sda, ida = theo_spec.getStringDataArrays(), theo_spec.getIntegerDataArrays()
                    names = [n.decode() for n in sda[0]] if sda else [""] * theo_spec.size()
                    charges = list(ida[0]) if ida else [None] * theo_spec.size()
                    theo_df = pd.DataFrame({"mz":[p.getMZ() for p in theo_spec], "intensity":[p.getIntensity() for p in theo_spec], 
                                          "annotation":names, "charge":charges})
                    theo_df = theo_df.query("annotation != ''")
                    
                    # Align spectra (exactly like notebook)
                    if psm.spectrum is None or theo_spec is None or psm.spectrum.size() == 0 or theo_spec.size() == 0:
                        continue
                    if psm.spectrum.size() > 1: 
                        psm.spectrum.sortByPosition()
                    if theo_spec.size() > 1: 
                        theo_spec.sortByPosition()
                    
                    try:
                        aligner = oms.SpectrumAlignment()
                        prm = aligner.getParameters()
                        prm.setValue("tolerance", mz_tolerance)
                        prm.setValue("is_relative_tolerance","false")
                        aligner.setParameters(prm)
                        pairs = []
                        aligner.getSpectrumAlignment(pairs, psm.spectrum, theo_spec)
                    except:
                        pairs = []
                    
                    # Annotate peaks (exactly like notebook)
                    if psm.spectrum is None or not pairs:
                        continue
                    theo_idxed = theo_df.reset_index().rename(columns={"index":"theo_idx"})
                    rows = [{"exp_idx":e, "theo_idx":t} for e,t in pairs]
                    aligned_df = pd.DataFrame(rows).merge(theo_idxed, on="theo_idx", how="left")
                    aligned_df["exp_mz"] = aligned_df.exp_idx.map(lambda i: psm.spectrum[i].getMZ())
                    aligned_df["exp_intensity"] = aligned_df.exp_idx.map(lambda i: psm.spectrum[i].getIntensity())
                    
                    if not aligned_df.empty:
                        aligned_df["psm_id"] = psm.spectrum_id
                        aligned_rows.append(aligned_df.assign(run=psm.run))
                except Exception as e:
                    logger.warning(f"Error processing PSM {psm.spectrum_id}: {e}")
                    continue
            
            if not aligned_rows:
                st.error("No successful alignments could be generated.")
                return
            
            # Combine all alignments
            aligned_df = pd.concat(aligned_rows, ignore_index=True, sort=False)
            st.success(f"✅ Generated {len(aligned_df)} peak alignments")
        
        # Step 6: Generate plots (following notebook exactly)
        st.subheader("📊 Experimental vs Predicted Spectra")
        
        exp_pred_plots = []
        for psm in psms_with_spectra:
            try:
                plot = create_exp_vs_pred_plot_streamlit(psm, aligned_df, selected_model)
                if plot:
                    exp_pred_plots.append(plot)
            except Exception as e:
                logger.error(f"Error creating exp vs pred plot for {psm.spectrum_id}: {e}")
        
        # Display exp vs pred plots
        if exp_pred_plots:
            for i, plot in enumerate(exp_pred_plots):
                st.bokeh_chart(plot, use_container_width=True)
        else:
            st.warning("No experimental vs predicted plots could be generated.")
        
        # Step 7: Generate replicate comparison plots (if multiple PSMs)
        exp_exp_plots = []
        if len(psms_with_spectra) >= 2:
            st.subheader("🔄 Replicate Comparisons")
            
            import itertools
            psm_combinations = list(itertools.combinations(psms_with_spectra, 2))
            
            for psm1, psm2 in psm_combinations:
                try:
                    plot = create_exp_vs_exp_plot_streamlit(psm1, psm2, aligned_df)
                    if plot:
                        exp_exp_plots.append(plot)
                except Exception as e:
                    logger.error(f"Error creating exp vs exp plot: {e}")
            
            # Display exp vs exp plots
            if exp_exp_plots:
                for plot in exp_exp_plots:
                    st.bokeh_chart(plot, use_container_width=True)
            else:
                st.warning("No replicate comparison plots could be generated.")
        
        # Analysis summary
        st.subheader("📋 Analysis Summary")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Target Peptide", selected_peptide)
        with col2:
            st.metric("PSMs with Spectra", len(psms_with_spectra))
        with col3:
            st.metric("Peak Alignments", len(aligned_df))
        with col4:
            st.metric("Total Plots", len(exp_pred_plots) + len(exp_exp_plots))
        
        # Option to restart analysis
        if st.button("🔄 Analyze Different Peptide"):
            st.session_state['analysis_ready'] = False
            st.rerun()
    
    except Exception as e:
        st.error(f"Error in automated analysis: {e}")
        logger.error(f"Automated analysis error: {e}", exc_info=True)

def create_exp_vs_pred_plot_streamlit(psm, aligned_df, model):
    """Create experimental vs predicted plot for Streamlit (based on notebook)."""
    try:
        # Extract experimental spectrum data
        psm_aligned = aligned_df[aligned_df.psm_id == psm.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        exp_df = st_tools.extract_annotations(psm.spectrum, psm_aligned)
        
        # Get predicted data
        pred_df = st_tools.predict_koina_ions(str(psm.peptidoform), psm.peptidoform.precursor_charge, model)
        
        # Create mirror plot
        title = f"Exp vs Pred: {psm.run} ({psm.peptidoform})"
        plot = st_tools.create_mirror_plot(
            exp_df, pred_df, 
            title=title,
            upper_colors=st_tools.EXP_COLORS,
            lower_colors=st_tools.PRED_COLORS,
            width=800, height=600
        )
        return plot
    except Exception as e:
        logger.error(f"Error creating exp vs pred plot: {e}")
        return None

def create_exp_vs_exp_plot_streamlit(psm1, psm2, aligned_df):
    """Create experimental vs experimental plot for Streamlit (based on notebook)."""
    try:
        # Extract data for both PSMs
        psm1_aligned = aligned_df[aligned_df.psm_id == psm1.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        exp_df1 = st_tools.extract_annotations(psm1.spectrum, psm1_aligned)
        
        psm2_aligned = aligned_df[aligned_df.psm_id == psm2.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        exp_df2 = st_tools.extract_annotations(psm2.spectrum, psm2_aligned)
        
        # Create mirror plot
        title = f"Exp vs Exp: {psm1.run} ({psm1.peptidoform}) vs {psm2.run} ({psm2.peptidoform})"
        plot = st_tools.create_mirror_plot(
            exp_df1, exp_df2,
            title=title,
            upper_colors=st_tools.EXP_COLORS,
            lower_colors=st_tools.EXP_COLORS,
            width=800, height=600
        )
        return plot
    except Exception as e:
        logger.error(f"Error creating exp vs exp plot: {e}")
        return None

if __name__ == "__main__":
    main()