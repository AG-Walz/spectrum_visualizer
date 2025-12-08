"""Code deduplication utilities - Eliminate repetitive patterns.

Consolidates 200+ lines of duplicate code into reusable functions: progress tracking,
file display, experiment management, directory operations. Used by both workspaces.
Example: create_experiment_map() replaces 2 duplicate mzML loading implementations.
"""

import streamlit as st
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from typing import List, Dict, Tuple, Optional, Callable, Any
import logging

logger = logging.getLogger(__name__)


# ================================================================
# PROGRESS TRACKING UTILITIES
# ================================================================

def process_items_with_progress(
    items: List[Any],
    process_func: Callable,
    status_template: str = "Processing {name}...",
    get_name: Callable = lambda x: str(x),
    show_progress: bool = True
) -> List[Any]:
    """
    Process a list of items with progress tracking.
    
    Args:
        items: List of items to process
        process_func: Function to apply to each item
        status_template: Template string for status messages (use {name} placeholder)
        get_name: Function to extract display name from item
        show_progress: Whether to show progress bar
        
    Returns:
        List of processed results
    """
    if not items:
        return []
    
    results = []
    
    if show_progress:
        progress_bar = st.progress(0)
        status_text = st.empty()
    
    for i, item in enumerate(items):
        if show_progress:
            name = get_name(item)
            status_text.text(status_template.format(name=name, index=i+1, total=len(items)))
        
        result = process_func(item)
        results.append(result)
        
        if show_progress:
            progress_bar.progress((i + 1) / len(items))
    
    if show_progress:
        status_text.text("✅ Processing complete")
        progress_bar.progress(1.0)
    
    return results


# ================================================================
# FILE DISPLAY UTILITIES
# ================================================================

def display_file_list(
    files: List[str],
    status_type: str = "info",
    header: Optional[str] = None,
    show_bullet: bool = True,
    max_display: Optional[int] = None
) -> None:
    """
    Display a list of files with consistent formatting.
    
    Args:
        files: List of file paths or names
        status_type: One of "success", "error", "warning", "info"
        header: Optional header message
        show_bullet: Whether to show bullet points
        max_display: Maximum number of files to display (None = all)
    """
    if not files:
        return
    
    display_files = files[:max_display] if max_display else files
    remaining = len(files) - len(display_files) if max_display and len(files) > max_display else 0
    
    if header:
        status_func = getattr(st, status_type, st.info)
        status_func(header)
    
    bullet = "• " if show_bullet else ""
    for f in display_files:
        st.write(f"{bullet}{f}")
    
    if remaining > 0:
        st.caption(f"...and {remaining} more files")


def display_categorized_files(
    matched: List[str],
    missing: List[str],
    extra: List[str],
    matched_header: str = "✅ Matched files:",
    missing_header: str = "❌ Missing files:",
    extra_header: str = "⚠️ Extra files:"
) -> None:
    """Display files categorized into matched, missing, and extra."""
    if matched:
        display_file_list(matched, "success", f"**{len(matched)}** {matched_header}")
    
    if missing:
        display_file_list(missing, "error", f"**{len(missing)}** {missing_header}")
    
    if extra:
        display_file_list(extra, "warning", f"**{len(extra)}** {extra_header}")


# ================================================================
# DIRECTORY UTILITIES
# ================================================================

def validate_directory(
    directory: str,
    create_if_missing: bool = False,
    show_ui: bool = True
) -> Tuple[bool, Optional[str]]:
    """
    Validate that a directory exists and is accessible.
    
    Args:
        directory: Path to directory
        create_if_missing: Whether to create the directory if it doesn't exist
        show_ui: Whether to show UI messages
        
    Returns:
        (is_valid, error_message)
    """
    path = Path(directory)
    
    if not path.exists():
        if create_if_missing:
            try:
                path.mkdir(parents=True, exist_ok=True)
                if show_ui:
                    st.success(f"✅ Created directory: {directory}")
                return True, None
            except Exception as e:
                error_msg = f"Failed to create directory: {e}"
                if show_ui:
                    st.error(f"❌ {error_msg}")
                return False, error_msg
        else:
            error_msg = f"Directory not found: {directory}"
            if show_ui:
                st.error(f"❌ {error_msg}")
            return False, error_msg
    
    if not path.is_dir():
        error_msg = f"Path exists but is not a directory: {directory}"
        if show_ui:
            st.error(f"❌ {error_msg}")
        return False, error_msg
    
    if show_ui:
        st.success(f"📂 Current directory: {directory}")
    
    return True, None


# ================================================================
# EXPERIMENT/SPECTRUM LOOKUP UTILITIES
# ================================================================

def create_experiment_map(mzml_paths: List[str], show_ui: bool = True) -> Dict[str, Any]:
    """
    Load mzML files into experiments dictionary.
    
    Args:
        mzml_paths: List of paths to mzML files
        show_ui: Whether to show UI progress
        
    Returns:
        Dict mapping filename to MSExperiment object
    """
    from shared_imports import oms
    
    experiments = {}
    
    def load_experiment(mzml_file: str) -> Tuple[str, Any]:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_file, exp)
        name = Path(mzml_file).name
        logger.info(f"  Loaded {name}: {exp.getNrSpectra()} spectra")
        return name, exp
    
    if show_ui:
        results = process_items_with_progress(
            mzml_paths,
            load_experiment,
            "Loading {name}...",
            lambda x: Path(x).name
        )
    else:
        results = [load_experiment(path) for path in mzml_paths]
    
    experiments = dict(results)
    
    total_spectra = sum(exp.getNrSpectra() for exp in experiments.values())
    logger.info(f"✅ Total experiments loaded: {len(experiments)}")
    logger.info(f"✅ Total spectra available: {total_spectra}")
    
    if show_ui:
        st.info(f"📂 Loaded {len(experiments)} mzML file(s) with {total_spectra} total spectra")
    
    return experiments


def create_spectrum_lookups(experiments: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create spectrum lookup objects for each experiment with smart format detection.
    
    Args:
        experiments: Dict mapping filename to MSExperiment
        
    Returns:
        Dict mapping file stem to SpectrumLookup object
    """
    from shared_imports import oms
    
    lookup_map = {}
    
    for exp_name, exp in experiments.items():
        spectrum_lookup = oms.SpectrumLookup()
        
        if exp.getNrSpectra() > 0:
            first_spectrum = exp.getSpectrum(0)
            native_id = first_spectrum.getNativeID()
            
            # Determine the appropriate regex pattern based on the native ID format
            if 'scan=' in native_id:
                pattern = "scan=(?<SCAN>\\d+)"
                logger.info(f"{exp_name}: Using scan-based pattern")
            elif 'index=' in native_id:
                pattern = "index=(?<SCAN>\\d+)"
                logger.info(f"{exp_name}: Using index-based pattern")
            elif 'spectrum=' in native_id:
                pattern = "spectrum=(?<SCAN>\\d+)"
                logger.info(f"{exp_name}: Using spectrum-based pattern")
            else:
                # Fallback: try to extract any number
                pattern = "(?<SCAN>\\d+)"
                logger.warning(f"{exp_name}: Using fallback pattern for unknown format: {native_id}")
            
            try:
                spectrum_lookup.readSpectra(exp, pattern)
            except Exception as e:
                logger.error(f"Failed to read spectra for {exp_name} with pattern {pattern}: {e}")
                # Try absolute fallback
                try:
                    spectrum_lookup.readSpectra(exp, "(?<SCAN>\\d+)")
                except Exception as e2:
                    logger.error(f"Fallback also failed for {exp_name}: {e2}")
        else:
            logger.warning(f"{exp_name}: No spectra found in experiment")
        
        lookup_map[Path(exp_name).stem] = spectrum_lookup
    
    logger.info(f"✅ Created {len(lookup_map)} spectrum lookups")
    return lookup_map


def create_fuzzy_match_helpers(experiments: Dict[str, Any]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Create helper dictionaries for fuzzy filename matching.
    
    Args:
        experiments: Dict mapping filename to experiment
        
    Returns:
        (stem_to_name, base_to_name) mapping dictionaries
    """
    exp_stem_to_name = {Path(exp_name).stem: exp_name for exp_name in experiments.keys()}
    exp_base_to_name = {}
    
    for exp_name in experiments.keys():
        base = Path(exp_name).stem.split('.')[0]
        exp_base_to_name[base] = exp_name
    
    return exp_stem_to_name, exp_base_to_name


def find_matching_experiment(
    run_name: str,
    experiments: Dict[str, Any],
    exp_stem_to_name: Dict[str, str],
    exp_base_to_name: Dict[str, str]
) -> Tuple[Optional[str], Optional[str]]:
    """
    Find matching experiment using multiple strategies.
    
    Args:
        run_name: The run name from PSM
        experiments: Experiment dictionary
        exp_stem_to_name: Stem to name mapping
        exp_base_to_name: Base to name mapping
        
    Returns:
        (experiment_name, lookup_key) or (None, None) if not found
    """
    # Strategy 1: Direct match with .mzML extension
    if f"{run_name}.mzML" in experiments:
        return f"{run_name}.mzML", run_name
    
    # Strategy 2: Direct match without extension
    if run_name in experiments:
        return run_name, Path(run_name).stem
    
    # Strategy 3: Match by stem
    if run_name in exp_stem_to_name:
        exp_name = exp_stem_to_name[run_name]
        return exp_name, run_name
    
    # Strategy 4: Match stem of run_name
    run_stem = Path(run_name).stem
    if run_stem in exp_stem_to_name:
        exp_name = exp_stem_to_name[run_stem]
        return exp_name, run_stem
    
    # Strategy 5: Match by base name (before first dot)
    run_base = Path(run_name).stem.split('.')[0]
    if run_base in exp_base_to_name:
        exp_name = exp_base_to_name[run_base]
        return exp_name, Path(exp_name).stem
    
    return None, None


# ================================================================
# STATUS MESSAGE UTILITIES
# ================================================================

def show_count_message(
    count: int,
    item_name: str,
    status_type: str = "info",
    additional_info: str = ""
) -> None:
    """
    Display a count message with consistent formatting.
    
    Args:
        count: Number of items
        item_name: Name of items (will be pluralized if count != 1)
        status_type: "success", "error", "warning", or "info"
        additional_info: Additional information to append
    """
    plural = "s" if count != 1 else ""
    message = f"{count} {item_name}{plural}"
    
    if additional_info:
        message = f"{message} {additional_info}"
    
    status_func = getattr(st, status_type, st.info)
    status_func(message)


def show_file_stats(
    file_list: List[str],
    total_size_mb: Optional[float] = None,
    status_type: str = "success"
) -> None:
    """
    Display file selection statistics.
    
    Args:
        file_list: List of file paths
        total_size_mb: Optional total size in MB
        status_type: "success", "info", "warning"
    """
    count = len(file_list)
    message = f"✅ Selected {count} file{'s' if count != 1 else ''}"
    
    if total_size_mb is not None:
        message += f" ({total_size_mb:.1f} MB)"
    
    status_func = getattr(st, status_type, st.info)
    status_func(message)


# ================================================================
# PSM FILTERING UTILITIES
# ================================================================

def filter_psms_by_spectra(
    psm_list: List[Any],
    selected_spectra_indices: set,
    mapping_df: Any
) -> List[Any]:
    """
    Filter PSMs to only include those with selected spectrum indices.
    
    Args:
        psm_list: List of PSM objects
        selected_spectra_indices: Set of selected spectrum indices
        mapping_df: DataFrame containing spectrum mappings
        
    Returns:
        Filtered list of PSMs
    """
    filtered_psms = []
    
    for psm in psm_list:
        # Try to find spectrum in mapping
        spectrum_matches = mapping_df[
            (mapping_df['spectrum_id'] == psm.spectrum_id) &
            (mapping_df['peptidoform'] == str(psm.peptidoform))
        ]
        
        # If found, check if spectrum_index is in selected set
        if not spectrum_matches.empty:
            idx = spectrum_matches.iloc[0].get('spectrum_index')
            if idx in selected_spectra_indices:
                filtered_psms.append(psm)
    
    return filtered_psms


def normalize_peptide_sequence(peptide: str) -> str:
    """
    Normalize a peptide sequence by removing modifications and charges.
    
    Args:
        peptide: Raw peptide string
        
    Returns:
        Normalized peptide string (uppercase, no modifications)
    """
    import re
    
    # Remove charge if present
    if '/' in peptide:
        peptide = peptide.split('/')[0]
    
    # Remove modifications in brackets and parentheses
    peptide = re.sub(r'\[.*?\]', '', peptide)
    peptide = re.sub(r'\(.*?\)', '', peptide)
    
    # Strip whitespace and convert to uppercase
    return peptide.strip().upper()
