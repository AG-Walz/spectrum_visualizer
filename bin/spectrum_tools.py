"""
spectrum_tools.py - Compatibility layer for existing Streamlit app
"""

import logging
from typing import List, Dict, Optional, Any, Tuple
import pyopenms as oms
import pandas as pd
from pathlib import Path
import tempfile
import itertools

# Import all functions from SpecViz
from .SpecViz import (
    filter_psms_by_peptide,
    attach_spectrum_to_psm,
    generate_theoretical_spectrum,
    align_spectra,
    annotate_peaks,
    extract_annotations,
    load_mzml_experiments,
    create_spectrum_lookup,
    convert_spectrum_to_df,
    predict_koina_ions,
    create_mirror_plot,
    INCLUDED_IONS,
    EXP_COLORS,
    PRED_COLORS
)

logger = logging.getLogger(__name__)

# Default parameters - inherit from SpecViz
DEFAULT_MODEL = "ms2pip_timsTOF2024"

# ================================================================
# STREAMLIT-SPECIFIC UTILITY FUNCTIONS
# ================================================================

def get_available_models() -> Dict[str, str]:
    """Get available models for Streamlit interface."""
    return {
        "MS2PIP TimsTOF 2024": "ms2pip_timsTOF2024",
        "Prosit Transformer": "Prosit_2023_intensity_timsTOF",
        "AlphaPept": "AlphaPept_ms2_generic"
    }


# ================================================================
# LEGACY COMPATIBILITY FUNCTIONS
# ================================================================
def create_comparative_spectrum_mirror_plot(upper_data, lower_data, title):
    """Create comparative spectrum mirror plot."""
    try:
        from bokeh.plotting import figure
        from bokeh.models import HoverTool, ColumnDataSource, Range1d
        # Setup plot with custom hover tool
        hover = HoverTool(tooltips=[
            ("m/z", "@mz{0.000}"),
            ("Intensity", "@intensity{0.000}"),
            ("Annotation", "@annotation")
        ])
        p = figure(width=800, height=500, title=title,
                   tools=["pan", "wheel_zoom", "box_zoom", "reset", "save", hover],
                   x_axis_label="m/z", y_axis_label="Relative Intensity")
        # Process upper spectrum
        if not upper_data.empty:
            upper_data = upper_data.copy()
            max_intensity = upper_data["intensity"].max()
            if max_intensity > 0:
                upper_data["intensity"] = upper_data["intensity"] / max_intensity
            upper_data["spectrum_type"] = "Experimental"
            # Create data source and plot upper spectrum
            upper_source = ColumnDataSource(upper_data)
            p.segment(x0='mz', y0=0, x1='mz', y1='intensity',
                     line_color="blue", line_width=2, alpha=0.8,
                     source=upper_source)
        # Process lower spectrum (mirrored)
        if not lower_data.empty:
            lower_data = lower_data.copy()
            max_intensity = lower_data["intensity"].max()
            if max_intensity > 0:
                lower_data["intensity"] = -lower_data["intensity"] / max_intensity
            lower_data["spectrum_type"] = "Predicted"
            # Create data source and plot lower spectrum
            lower_source = ColumnDataSource(lower_data)
            p.segment(x0='mz', y0=0, x1='mz', y1='intensity',
                     line_color="red", line_width=2, alpha=0.8,
                     source=lower_source)
        # Set ranges
        p.y_range = Range1d(-1.2 if not lower_data.empty else 0, 1.2)
        # Add center line - fix concat warning
        if not (upper_data.empty and lower_data.empty):
            mz_dfs = []
            if not upper_data.empty:
                mz_dfs.append(upper_data[["mz"]])
            if not lower_data.empty:
                mz_dfs.append(lower_data[["mz"]])
            if mz_dfs:
                all_mz = pd.concat(mz_dfs, sort=False).dropna()["mz"]
                if len(all_mz) > 0:
                    p.line([all_mz.min(), all_mz.max()], [0, 0],
                           line_color="gray", line_width=1, alpha=0.5)
        return p
    except Exception as e:
        logger.error(f"Error creating mirror plot: {e}")
        return None

def create_single_spectrum_plot(mzml_path: str, idxml_path: str, spectrum_index: int,
                               model: str = DEFAULT_MODEL) -> Optional[Any]:
    """Create a single spectrum plot for Streamlit integration."""
    try:
        from psm_utils.io import read_file
        
        # Load experiment
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        
        if spectrum_index >= exp.size():
            logger.error(f"Spectrum index {spectrum_index} out of range (max: {exp.size()-1})")
            return None
        
        # Get spectrum
        spectrum = exp[spectrum_index]
        
        # Load PSM data to find matching peptide
        psm_list = read_file(idxml_path, file_format="idXML")
        
        # Find PSM for this spectrum
        matching_psm = None
        for psm in psm_list:
            if hasattr(psm, 'spectrum_id'):
                try:
                    if int(psm.spectrum_id.split('scan=')[1]) == spectrum_index + 1:
                        matching_psm = psm
                        break
                except:
                    pass
        
        # Get experimental data
        exp_df = convert_spectrum_to_df(spectrum)
        if exp_df.empty:
            exp_df = pd.DataFrame({"mz": [0], "intensity": [0], "annotation": ["no data"]})
        else:
            exp_df["annotation"] = "unassigned"
        
        # Get predicted data if PSM found
        pred_df = pd.DataFrame(columns=['mz', 'intensity', 'annotation'])
        if matching_psm:
            try:
                pred_df = predict_koina_ions(str(matching_psm.peptidoform), 
                                           matching_psm.peptidoform.precursor_charge, model)
            except Exception as e:
                logger.error(f"Prediction failed: {e}")
        
        title = f"Spectrum {spectrum_index}"
        if matching_psm:
            title += f": {matching_psm.peptidoform.sequence}"
        
        return create_mirror_plot(exp_df, pred_df, title=title)
        
    except Exception as e:
        logger.error(f"Error creating spectrum plot: {e}")
        return None

logger = logging.getLogger(__name__)

# Legacy function names for backward compatibility
def load_experiments(file_paths: List[str]) -> Dict[str, oms.MSExperiment]:
    """Legacy wrapper for load_mzml_experiments."""
    return load_mzml_experiments(file_paths)

def create_lookups(experiments: Dict[str, oms.MSExperiment]) -> Dict[str, oms.SpectrumLookup]:
    """Legacy wrapper for create_spectrum_lookup."""
    return create_spectrum_lookup(experiments)

def get_spectrum_data(spectrum: oms.MSSpectrum) -> pd.DataFrame:
    """Legacy wrapper for convert_spectrum_to_df."""
    return convert_spectrum_to_df(spectrum)

def create_plot(exp_data: pd.DataFrame, pred_data: pd.DataFrame, title: str = "Spectrum Plot") -> Any:
    """Legacy wrapper for create_mirror_plot."""
    return create_mirror_plot(exp_data, pred_data, title=title)

def get_koina_predictions(peptide: str, charge: int, model: str = DEFAULT_MODEL) -> pd.DataFrame:
    """Legacy wrapper for predict_koina_ions."""
    return predict_koina_ions(peptide, charge, model)

# Additional utility functions that might be needed by the Streamlit app
def validate_mzml_file(file_path: str) -> bool:
    """Validate if a file is a valid mzML file."""
    try:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file_path, exp)
        return exp.size() > 0
    except Exception as e:
        logger.error(f"Invalid mzML file {file_path}: {e}")
        return False

def validate_idxml_file(file_path: str) -> bool:
    """Validate if a file is a valid idXML file."""
    try:
        from psm_utils.io import read_file
        psm_list = read_file(file_path, file_format="idXML")
        return len(psm_list) > 0
    except Exception as e:
        logger.error(f"Invalid idXML file {file_path}: {e}")
        return False

def get_file_info(file_path: str) -> Dict[str, Any]:
    """Get basic information about a data file."""
    try:
        path = Path(file_path)
        if path.suffix.lower() == '.mzml':
            exp = oms.MSExperiment()
            oms.MzMLFile().load(file_path, exp)
            return {
                "type": "mzML",
                "size": exp.size(),
                "valid": True,
                "spectra_count": exp.size()
            }
        elif path.suffix.lower() == '.idxml':
            from psm_utils.io import read_file
            psm_list = read_file(file_path, file_format="idXML")
            return {
                "type": "idXML",
                "size": len(psm_list),
                "valid": True,
                "psm_count": len(psm_list)
            }
        else:
            return {"type": "unknown", "valid": False}
    except Exception as e:
        logger.error(f"Error getting file info for {file_path}: {e}")
        return {"type": "error", "valid": False, "error": str(e)}

def get_available_models_list() -> List[str]:
    """Get list of available Koina models (simple list version)."""
    return [
        "ms2pip_timsTOF2024",
        "ms2pip_Bruker_timsTOF_DDA_Orbi_HCD",
        "ms2pip_TMT",
        "ms2pip_iTRAQ",
        "ms2pip_CID",
        "ms2pip_HCD",
        "ms2pip_ETD"
    ]

def ensure_streamlit_compatible_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure DataFrame is compatible with Streamlit/Arrow serialization.
    
    Args:
        df: Input DataFrame
        
    Returns:
        DataFrame with proper data types for Streamlit display
    """
    if df.empty:
        return df
        
    df_copy = df.copy()
    
    # Convert columns to appropriate types
    for col in df_copy.columns:
        if col in ['mz', 'intensity', 'exp_mz', 'exp_intensity', 'theo_mz', 'theo_intensity']:
            df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce').astype(float)
        elif col in ['annotation', 'psm_id', 'run']:
            df_copy[col] = df_copy[col].astype(str)
        elif col in ['exp_idx', 'theo_idx']:
            df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce').astype('Int64')
        elif col == 'Value':  # This is causing the Arrow error
            # Convert spectrum IDs that look like integers but are strings
            df_copy[col] = df_copy[col].astype(str)
    
    return df_copy

def get_spectrum_count(file_path: str) -> int:
    """Get number of spectra in an mzML file."""
    try:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file_path, exp)
        return exp.size()
    except Exception:
        return 0

def get_psm_count(file_path: str) -> int:
    """Get number of PSMs in an idXML file."""
    try:
        from psm_utils.io import read_file
        psm_list = read_file(file_path, file_format="idXML")
        return len(psm_list)
    except Exception:
        return 0

def extract_unique_peptides(file_path: str) -> List[str]:
    """Extract unique peptide sequences from idXML file."""
    try:
        from psm_utils.io import read_file
        from psm_utils import Peptidoform
        
        psm_list = read_file(file_path, file_format="idXML")
        unique_peptides = set()
        
        for psm in psm_list:
            try:
                peptide_seq = Peptidoform(str(psm.peptidoform)).sequence
                unique_peptides.add(peptide_seq)
            except Exception:
                continue
        
        return sorted(list(unique_peptides))
    except Exception as e:
        logger.error(f"Error extracting peptides from {file_path}: {e}")
        return []

# Additional functions required by the Streamlit app
def extract_mzml_files_from_idxml(idxml_path: str) -> List[str]:
    """Extract mzML file names referenced in idXML file."""
    try:
        from psm_utils.io import read_file
        psm_list = read_file(idxml_path, file_format="idXML")
        
        # Extract unique run names from PSMs
        run_names = set()
        for psm in psm_list:
            if hasattr(psm, 'run') and psm.run:
                run_names.add(f"{psm.run}.mzML")
        
        return sorted(list(run_names))
    except Exception as e:
        logger.error(f"Error extracting mzML files from idXML: {e}")
        return []

def validate_mzml_files_mapping(uploaded_files: List[str], expected_files: List[str]) -> tuple:
    """Validate mapping between uploaded and expected mzML files."""
    uploaded_set = set(uploaded_files)
    expected_set = set(expected_files)
    
    matched_files = list(uploaded_set.intersection(expected_set))
    missing_files = list(expected_set - uploaded_set)
    extra_files = list(uploaded_set - expected_set)
    
    return matched_files, missing_files, extra_files

def load_identifications(idxml_path: str) -> tuple:
    """Load protein and peptide identifications from idXML file."""
    try:
        from psm_utils.io import read_file
        psm_list = read_file(idxml_path, file_format="idXML")
        
        # Extract protein IDs (simplified - just unique proteins)
        protein_ids = set()
        peptide_ids = []
        
        for psm in psm_list:
            if hasattr(psm, 'peptidoform'):
                peptide_ids.append(psm)
                # Try to extract protein info if available
                if hasattr(psm, 'protein'):
                    protein_ids.add(psm.protein)
        
        return list(protein_ids), peptide_ids
    except Exception as e:
        logger.error(f"Error loading identifications: {e}")
        return [], []

def load_experiment(mzml_path: str) -> oms.MSExperiment:
    """Load single mzML experiment."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml_path, exp)
    return exp

def extract_spectrum_mapping_multi(mzml_paths: List[str], idxml_path: str) -> pd.DataFrame:
    """Extract spectrum-peptide mapping from multiple mzML files and idXML."""
    try:
        from psm_utils.io import read_file
        psm_list = read_file(idxml_path, file_format="idXML")
        
        # Load experiments to get spectrum information
        experiments = {}
        for mzml_path in mzml_paths:
            exp = load_experiment(mzml_path)
            experiments[Path(mzml_path).stem] = exp
        
        mappings = []
        for i, psm in enumerate(psm_list):
            try:
                # Ensure proper data types for Arrow compatibility
                spectrum_id = str(psm.spectrum_id) if hasattr(psm, 'spectrum_id') else f"scan_{i}"
                peptide_seq = str(psm.peptidoform) if hasattr(psm, 'peptidoform') else "Unknown"
                charge = int(psm.peptidoform.precursor_charge) if hasattr(psm, 'peptidoform') else 1
                run = str(psm.run) if hasattr(psm, 'run') else "unknown"
                score = float(psm.score) if hasattr(psm, 'score') else 0.0
                mz = float(psm.precursor_mz) if hasattr(psm, 'precursor_mz') else 0.0
                rt = float(psm.retention_time) if hasattr(psm, 'retention_time') else 0.0
                
                mapping = {
                    "index": i,
                    "spectrum_id": spectrum_id,
                    "peptide_sequence": peptide_seq,
                    "charge": charge,
                    "run": run,
                    "score": score,
                    "mz": mz,
                    "rt": rt
                }
                mappings.append(mapping)
            except Exception as e:
                logger.warning(f"Error processing PSM {i}: {e}")
                continue
        
        df = pd.DataFrame(mappings)
        return ensure_streamlit_compatible_dataframe(df)
    except Exception as e:
        logger.error(f"Error extracting spectrum mapping: {e}")
        return pd.DataFrame()

def convert_to_mass_shift_format(peptide_seq: str) -> str:
    """Convert peptide sequence to mass shift format."""
    # Simple conversion - in practice this might need more sophisticated handling
    return str(peptide_seq)

def visualize_spectrum(mzml_path: str, id_path: str, spectrum_index: int, 
                      model_name: str = DEFAULT_MODEL, width: int = 800, 
                      height: int = 500, debug: bool = False) -> Optional[Any]:
    """Create spectrum visualization for Streamlit app."""
    try:
        # Use the create_single_spectrum_plot function from SpecViz
        plot = create_single_spectrum_plot(mzml_path, id_path, spectrum_index, model_name)
        return plot
    except Exception as e:
        logger.error(f"Error in visualize_spectrum: {e}")
        if debug:
            raise
        return None

# Export functions that actually exist in this module
__all__ = [
    # Core functions from SpecViz (imported)
    'filter_psms_by_peptide', 'attach_spectrum_to_psm', 'generate_theoretical_spectrum',
    'align_spectra', 'annotate_peaks', 'extract_annotations',
    'load_mzml_experiments', 'create_spectrum_lookup', 'convert_spectrum_to_df',
    'predict_koina_ions', 'create_mirror_plot',
    
    # Legacy compatibility
    'load_experiments', 'create_lookups', 'get_spectrum_data', 'create_plot',
    'get_koina_predictions',
    
    # Streamlit app functions
    'extract_mzml_files_from_idxml', 'validate_mzml_files_mapping', 'load_identifications',
    'load_experiment', 'extract_spectrum_mapping_multi', 'convert_to_mass_shift_format',
    'visualize_spectrum', 'ensure_streamlit_compatible_dataframe',
    
    # Utility functions
    'validate_mzml_file', 'validate_idxml_file', 'get_file_info',
    'get_available_models', 'get_spectrum_count', 'get_psm_count',
    'extract_unique_peptides',
    
    # Constants from SpecViz (imported)
    'EXP_COLORS', 'PRED_COLORS', 'INCLUDED_IONS', 'DEFAULT_MODEL'
]