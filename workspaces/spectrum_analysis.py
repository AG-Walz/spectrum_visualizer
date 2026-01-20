"""Spectrum Analysis - Main interactive workspace for MS/MS visualization.

5-step pipeline: File upload → Peptide/spectrum selection → Analysis settings →
Automated analysis (Exp vs Predicted, All-vs-All, One-vs-All modes) → Export plots.
Memory-optimized: Only loads mzML files containing selected spectra. Supports Koina AI.
"""
import sys as _sys
from pathlib import Path as _Path

import streamlit as st

# Add paths
_sys.path.insert(0, str(_Path(__file__).parent.parent / "utils"))
_sys.path.insert(0, str(_Path(__file__).parent.parent / "core"))

# Import all heavy packages from shared module (loaded once at app startup)
# This replaces individual imports of pandas, numpy, pyopenms, psm_utils, etc.
from shared_imports import (
    pd, np, oms, read_file, PSM, PSMList, Peptidoform,
    warnings, psutil, tempfile, os, sys, Path, logging, json, gc, socket,
    Optional, Tuple, List, Dict, Any, Union, itertools, re, StringIO,
    PSUTIL_AVAILABLE, BOKEH_AVAILABLE, TABULATE_AVAILABLE,
    SPECTRUM_TOOLS_AVAILABLE,
    figure, file_html, CDN,
    # Common utility functions
    create_file_options_with_sizes, map_display_names_to_paths,
    extract_filename, extract_filenames, is_empty_collection,
    safe_concat_dataframes, show_status_message, validate_file_list,
    create_validation_info_dict,
    create_validation_info_dict
)

try:
    import streamlit_bridge as st_tools
    from streamlit_bridge import get_available_models
except ImportError:
    st_tools = None
    get_available_models = lambda: {}
import mmap
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import shared file selection components
from file_selection import (
    render_directory_browser,
    render_file_filter_and_selection,
    get_files_with_substring
)

# Import shared data input UI (single source of truth for file selection interface)
from data_input import (
    render_data_input_section,
    render_unified_file_browser,
    render_server_directory_selector,
    get_available_servers
)

# Import refactoring utilities for code deduplication
try:
    from workspace_utils import (
        process_items_with_progress, display_file_list, display_categorized_files,
        validate_directory, create_experiment_map, create_spectrum_lookups,
        create_fuzzy_match_helpers, find_matching_experiment,
        show_count_message, show_file_stats, filter_psms_by_spectra,
        normalize_peptide_sequence
    )
except ImportError:
    from utils.workspace_utils import (
        process_items_with_progress, display_file_list, display_categorized_files,
        validate_directory, create_experiment_map, create_spectrum_lookups,
        create_fuzzy_match_helpers, find_matching_experiment,
        show_count_message, show_file_stats, filter_psms_by_spectra,
        normalize_peptide_sequence
    )

# Import optional st_bokeh and tabulate if available
if BOKEH_AVAILABLE:
    from shared_imports import st_bokeh
if TABULATE_AVAILABLE:
    from shared_imports import tabulate

# Import spectrum engine functions for matplotlib-based plots
try:
    from spectrum_engine import (create_peptide_overview_plot, export_overview_plot,
                                 create_matplotlib_mirror_plot, export_matplotlib_mirror_plot)
except ImportError:
    from core.spectrum_engine import (create_peptide_overview_plot, export_overview_plot,
                                      create_matplotlib_mirror_plot, export_matplotlib_mirror_plot)

# Logger setup (already configured in shared_imports, just get reference)
logger = logging.getLogger(__name__)

# ================================================================
# UTILITY FUNCTIONS (CONSOLIDATED)
# ================================================================

def find_matching_mzml(idxml_stem: str, available_mzml_names: List[str]) -> Optional[str]:
    """
    Smart matching of idXML files to mzML files.
    
    Tries multiple strategies:
    1. Exact match (idxml_stem.mzML)
    2. Common prefix match (e.g., ProT013_1 matches ProT013_1_something.mzML)
    3. Strip suffixes like _pout_filtered, _ms2rescore and try again
    
    Args:
        idxml_stem: Stem of idXML filename (without extension)
        available_mzml_names: List of available mzML filenames
    
    Returns:
        Best matching mzML filename or None
    """
    # Strategy 1: Exact match
    exact_match = f"{idxml_stem}.mzML"
    if exact_match in available_mzml_names:
        return exact_match
    
    # Strategy 2: Strip common idXML suffixes and try again
    suffixes_to_strip = ['_pout_filtered', '_ms2rescore', '_filtered', '_rescored', '_percolator']
    stem_cleaned = idxml_stem
    for suffix in suffixes_to_strip:
        if stem_cleaned.endswith(suffix):
            stem_cleaned = stem_cleaned[:-len(suffix)]
            break
    
    # Try exact match with cleaned stem
    if stem_cleaned != idxml_stem:
        cleaned_match = f"{stem_cleaned}.mzML"
        if cleaned_match in available_mzml_names:
            return cleaned_match
            
    # CRITICAL: Strategy 3 (Containment) must be strict to avoid partial matches
    # e.g. "ProC" should not match "ProCa"
    
    # 3a. Check if the full cleaned stem appears in filename
    # Prefer matches where it is delimited by underscore or boundaries
    matches = []
    import re
    
    for f in available_mzml_names:
        # Check if stem_cleaned is in file
        if stem_cleaned in f:
            # Enforce that it's a distinct segment if possible
            # e.g. "ProT002" in "23..._ProT002_..." -> MATCH
            # "ProC" in "23..._ProCa_..." -> NO MATCH (ProCa != ProC)
            # Use regex to check for boundaries: (^|_)PATTERN($|_)
            escaped_stem = re.escape(stem_cleaned)
            if re.search(f"(^|_){escaped_stem}($|_|\\.)", f):
                matches.append(f)
    
    if matches:
        return min(matches, key=len)
    
    # 3b. Fallback: Main part (before first underscore or hyphen)
    # Also enforce strict boundaries
    separators = ['_', '-']
    main_part = stem_cleaned
    for sep in separators:
        if sep in stem_cleaned:
            main_part = stem_cleaned.split(sep)[0]
            break
            
    if main_part and len(main_part) > 2: # Don't search for tiny parts
        matches = []
        for f in available_mzml_names:
            # STRICT check: main_part must be surrounded by delimiters or boundaries
            # e.g. "ProC" must match "..._ProC_..." or "ProC_..."
            escaped_main = re.escape(main_part)
            if re.search(f"(^|_){escaped_main}($|_)", f):
                 matches.append(f)
                 
        if matches:
            return min(matches, key=len)

    return None

def validate_spectrum_index(mapping_df: pd.DataFrame, spectrum_index: int, context: str = "") -> bool:
    """Validate spectrum index."""
    try:
        return 0 <= spectrum_index < len(mapping_df)
    except: return False

def get_spectrum_row(mapping_df: pd.DataFrame, spectrum_index: int):
    """Get spectrum row from mapping dataframe."""
    if not validate_spectrum_index(mapping_df, spectrum_index):
        return None
    if "spectrum_index" in mapping_df.columns:
        matching_rows = mapping_df[mapping_df["spectrum_index"] == spectrum_index]
        return matching_rows.iloc[0] if not matching_rows.empty else None
    return mapping_df.iloc[spectrum_index]

@st.cache_data
def save_uploaded_file(uploaded_file, suffix: Optional[str] = None) -> str:
    """Save uploaded file to temporary directory."""
    temp_dir = tempfile.mkdtemp()
    file_name = getattr(uploaded_file, 'name', f"uploaded_file{suffix or ''}")
    if suffix and not file_name.endswith(suffix):
        file_name = f"{os.path.splitext(file_name)[0]}{suffix}"
    file_path = os.path.join(temp_dir, file_name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getvalue())
    return file_path

# ================================================================
# OPTIMIZED LOADING FUNCTIONS
# ================================================================

def extract_idxml_metadata(idxml_paths: List[str], progress_callback=None) -> Tuple[Dict[str, List[str]], Dict[str, Any], List[str], List[str]]:
    """
    PHASE 1: Extract metadata from idXML files without loading mzML.
    OPTIMIZED: Uses PyOpenMS directly for faster peptide extraction.
    
    Args:
        idxml_paths: List of idXML file paths
        progress_callback: Optional function(current, total, message) for progress updates
    
    Returns:
        (expected_mzmls_by_idxml, all_peptidoforms_data, messages, empty_files)
    """
    messages = []
    empty_files = []  # Track files with no peptide identifications
    expected_mzmls_by_idxml = {}  # {idxml_name: [expected_mzml_files]}
    all_peptidoforms_data = {}  # {peptidoform_str: {data}}
    
    total_files = len(idxml_paths)
    for idx, idxml_path in enumerate(idxml_paths):
        if progress_callback:
            progress_callback(idx, total_files, extract_filename(idxml_path))
        idxml_name = extract_filename(idxml_path)
        
        # Safety check: Skip non-idXML files
        if not idxml_path.lower().endswith(('.idxml', '.mzid', '.mzidentml')):
            messages.append(f"⚠️ Skipping non-idXML file: {idxml_name}")
            logger.warning(f"Skipping non-idXML file: {idxml_path}")
            continue
        
        try:
            # OPTIMIZED: Use PyOpenMS directly for faster extraction
            protein_ids = []
            peptide_ids = []
            oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)
            
            if not peptide_ids:
                raise Exception("contains no PeptideIdentifications")
            
            # OPTIMIZED: Extract unique peptidoforms (only first hit per peptide_id for speed)
            for pep_id in peptide_ids:
                hits = pep_id.getHits()
                if not hits:
                    continue
                try:
                    # OPTIMIZATION: Only process first hit (top-ranked)
                    sequence = hits[0].getSequence()
                    peptidoform_str = sequence.toString()
                    unmodified_seq = sequence.toUnmodifiedString()
                    
                    if peptidoform_str not in all_peptidoforms_data:
                        display_label = f"{unmodified_seq} [{peptidoform_str}]" if peptidoform_str != unmodified_seq else unmodified_seq
                        all_peptidoforms_data[peptidoform_str] = {
                            'peptidoform': peptidoform_str,
                            'sequence': unmodified_seq,
                            'display_label': display_label,
                            'psm_count': 0
                        }
                    
                    all_peptidoforms_data[peptidoform_str]['psm_count'] += 1
                    
                except Exception as e:
                    continue
                    
        except Exception as e:
            error_msg = str(e)
            if "contains no PeptideIdentifications" in error_msg or "Nothing to parse" in error_msg:
                empty_files.append(idxml_name)
            continue
    
    return expected_mzmls_by_idxml, all_peptidoforms_data, messages, empty_files


    return expected_mzmls_by_idxml, all_peptidoforms_data, messages, empty_files


def decode_metadata_value(value):
    """Decode metadata value that might be bytes, list of bytes, or string."""
    # Handle list of bytes: [b'filename']
    if isinstance(value, list) and len(value) > 0:
        value = value[0]
    # Handle bytes: b'filename'
    if isinstance(value, bytes):
        return value.decode("utf-8")
    # Already a string
    return str(value)


def extract_filename(path: str) -> str:
    """Extract filename from path - locally defined to avoid import issues."""
    return Path(path).name


def _fast_search_peptides(filepath: str, peptides: set) -> bool:
    """
    Rapidly check if any of the peptide strings exist in the file using mmap.
    This avoids parsing huge XML files that don't satisfy the filter.
    
    Args:
        filepath: Path to the file
        peptides: Set of peptide strings (plain unmodified sequences) to search for
        
    Returns:
        True if potential match found, False otherwise.
    """
    try:
        # Convert peptides to bytes for mmap search
        # We search for the sequence. In idXML/mzIdentML it usually appears as:
        # sequence="PEPTIDE" or <PeptideSequence>PEPTIDE</PeptideSequence>
        # Searching for the raw string is a safe heuristic (zero false negatives).
        
        # Optimize: Sort by length, maybe? Or just convert all.
        # Filter out very short peptides to avoid false positives (e.g. "K", "R")
        long_peptides = [p.encode('utf-8') for p in peptides if len(p) > 4]
        if not long_peptides:
            return True # If only short peptides, skip optimization to be safe
            
        with open(filepath, 'rb') as f:
            # Handle empty files
            if os.fstat(f.fileno()).st_size == 0:
                return False
                
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as s:
                for pep_bytes in long_peptides:
                    if s.find(pep_bytes) != -1:
                        return True
        return False
    except Exception as e:
        # If optimization fails (e.g. permissions, empty file), fallback to full parse
        # logger.warning(f"Fast search failed for {filepath}: {e}")
        return True


def _filter_single_idxml(args: Tuple) -> Tuple[str, Dict, List[str], bool]:
    """Helper for parallel processing of idXML filtering.
    
    Returns:
        (idxml_name, mzml_spectra, messages, is_empty)
    """
    idxml_path, available_mzml_names, selected_set, unmodified_selected = args
    messages = []
    is_empty = False
    idxml_name = extract_filename(idxml_path)
    mzml_spectra = {}  # {mzml_file: [(spec_id, peptidoform, peptide_idx)]}
    
    try:
        # OPTIMIZATION: Fast string search before XML parsing
        # This speeds up filtering by 10-100x for large datasets
        if not _fast_search_peptides(idxml_path, unmodified_selected):
             # messages.append(f"  ⏭️ Skipped {idxml_name} (pre-scan: no peptide matches)")
             return idxml_name, {}, messages, False

        # OPTIMIZED: Use PyOpenMS directly (same as extraction phase)
        protein_ids = []
        peptide_ids = []
        oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)
        
        if not peptide_ids:
            return idxml_name, {}, messages, True
        
        # FIX: Build a map of identifier -> mzML filename from ALL ProteinIdentification objects
        # Merged idXML files (e.g. from Percolator) contain multiple runs with different identifiers
        run_map = {}  # identifier -> mzML_filename
        default_mzml_file = None
        
        if protein_ids:
            for i, protein_id in enumerate(protein_ids):
                mzml_file = None
                metadata_keys = ["CometAdapter:1:in", "spectra_data", "spectrum_reference", "mzML_file"]
                for key in metadata_keys:
                    if protein_id.metaValueExists(key):
                        value = decode_metadata_value(protein_id.getMetaValue(key))
                        mzml_file = Path(value).name
                        break
                
                identifier = protein_id.getIdentifier()
                if mzml_file:
                    run_map[identifier] = mzml_file
                    # Keep the first one as default/fallback
                    if default_mzml_file is None:
                        default_mzml_file = mzml_file
        
        # If no metadata found at all, use smart matching as fallback for the default
        if not default_mzml_file:
            idxml_stem = Path(idxml_path).stem
            default_mzml_file = find_matching_mzml(idxml_stem, available_mzml_names)
            if not default_mzml_file:
                default_mzml_file = f"{idxml_stem}.mzML"
        
        mzml_files_seen = set()
        
        for peptide_idx, pep_id in enumerate(peptide_ids):
            hits = pep_id.getHits()
            if not hits:
                continue
            
            # OPTIMIZATION: Check if ANY hit matches selected peptides (early exit)
            match_found = False
            for hit in hits:
                try:
                    sequence = hit.getSequence()
                    unmodified_seq = sequence.toUnmodifiedString()
                    
                    # OPTIMIZATION: Check unmodified first (faster)
                    if unmodified_seq not in unmodified_selected:
                        continue
                    
                    # Now check the full peptidoform
                    peptidoform_str = sequence.toString()
                    
                    if peptidoform_str in selected_set or unmodified_seq in selected_set:
                        match_found = True
                        # Get spectrum ID
                        spectrum_id = None
                        if pep_id.metaValueExists("spectrum_reference"):
                            spectrum_id = decode_metadata_value(pep_id.getMetaValue("spectrum_reference"))
                        elif pep_id.metaValueExists("spectrum_id"):
                            spectrum_id = decode_metadata_value(pep_id.getMetaValue("spectrum_id"))
                        else:
                            spectrum_id = str(pep_id.getSpectrumReference())
                        
                        # Determine mzML file
                        # Priority 1: Explicit metadata on peptide
                        mzml_file = None
                        if pep_id.metaValueExists("spectra_data"):
                            spectra_data = decode_metadata_value(pep_id.getMetaValue("spectra_data"))
                            mzml_file = Path(spectra_data).name
                        elif pep_id.metaValueExists("mzML_file"):
                            mzml_file = decode_metadata_value(pep_id.getMetaValue("mzML_file"))
                        elif pep_id.metaValueExists("file_origin"): # Also check file_origin
                            file_origin = decode_metadata_value(pep_id.getMetaValue("file_origin"))
                            mzml_file = Path(file_origin).name
                        
                        # Priority 2: Inferred from Run/ProteinIdentifier
                        if not mzml_file:
                            identifier = pep_id.getIdentifier()
                            # Try exact match first
                            if identifier in run_map:
                                mzml_file = run_map[identifier]
                            else:
                                # Try finding a run that matches the prefix of the peptide identifier
                                # Percolator often appends _number to the run identifier for peptides
                                for run_id, run_file in run_map.items():
                                    if identifier.startswith(run_id):
                                        mzml_file = run_file
                                        break
                        
                        # Priority 3: Default fallback
                        if not mzml_file:
                            mzml_file = default_mzml_file
                        
                        mzml_files_seen.add(mzml_file)
                        
                        if mzml_file not in mzml_spectra:
                            mzml_spectra[mzml_file] = []
                        
                        # CLEANUP: Ensure spectrum_id matches mzML Native ID format
                        # Sometimes spectrum_reference is "file.mzML scan=123" or similar
                        # We want just "scan=123" or whatever the native ID pattern is
                        spec_id_str = None
                        if spectrum_id and spectrum_id.strip() and spectrum_id != "None":
                            raw_id = spectrum_id.strip()
                            clean_id = raw_id
                            
                            # Heuristic 1: If it contains "scan=", "index=", or "spectrum=", assume that starts the ID
                            for tag in ["scan=", "index=", "spectrum=", "controllerType="]:
                                if tag in clean_id:
                                    idx = clean_id.find(tag)
                                    # If the tag is not at the start, check if it's preceded by space or colon
                                    if idx > 0 and clean_id[idx-1] in [' ', ':']:
                                        clean_id = clean_id[idx:]
                                        break
                            
                            # Heuristic 2: If we didn't find a tag but it looks like "file.mzML 123", strip the file
                            # This handles "file.mzML nativeID" where nativeID is just a number
                            if clean_id == raw_id and (".mzML" in clean_id or ".mzXML" in clean_id):
                                 # Split by space
                                 parts = clean_id.split()
                                 if len(parts) > 1:
                                     # Check specific format "file.mzML id"
                                     # If the first part ends in .mzML, assume the rest is ID
                                     if parts[0].endswith(".mzML") or parts[0].endswith(".mzXML") or parts[0].endswith(".mzML\"") or parts[0].endswith(".mzXML\""):
                                         clean_id = " ".join(parts[1:])
                            
                            spec_id_str = clean_id
                        else:
                            spec_id_str = f"peptide_idx_{peptide_idx}"
                        
                        # CRITICAL: Store peptide_idx for O(1) access later
                        # Format: (spec_id, peptidoform, peptide_idx)
                        mzml_spectra[mzml_file].append((spec_id_str, peptidoform_str, peptide_idx))
                        
                        break  # Found a match, move to next peptide ID
                
                except Exception:
                    continue
                
                if match_found:
                    break
        
        # Generate messages
        if mzml_spectra:
            if len(mzml_files_seen) == 1:
                messages.append(f"  ✓ {idxml_name} → {list(mzml_files_seen)[0]}")
            else:
                messages.append(f"  ✓ {idxml_name} → {len(mzml_files_seen)} mzML files: {', '.join(sorted(mzml_files_seen))}")
        
        return idxml_name, mzml_spectra, messages, is_empty

    except Exception as e:
        error_msg = str(e)
        if "contains no PeptideIdentifications" in error_msg or "Nothing to parse" in error_msg:
            return idxml_name, {}, [], True
        return idxml_name, {}, [f"Error processing {idxml_name}: {e}"], False


def filter_idxml_by_selected_peptides(idxml_paths: List[str], 
                                      mzml_paths: List[str],
                                      selected_peptidoforms: List[str],
                                      progress_callback=None,
                                      max_workers: int = 8) -> Tuple[Dict[str, Dict[str, List[Tuple[str, str, int]]]], List[str], List[str]]:
    """
    PHASE 2: Filter idXML files to get only spectrum IDs matching selected peptides.
    
    Args:
        idxml_paths: List of idXML file paths
        mzml_paths: List of available mzML file paths
        selected_peptidoforms: List of selected peptide string
        progress_callback: Progress callback function
        max_workers: Number of parallel workers
    
    Returns:
        (idxml_to_mzml_spectra, messages, empty_files)
        Structure: {idxml_name: {mzml_file: [(spec_id, peptidoform, peptide_idx)]}}
    """
    messages = []
    empty_files = []
    idxml_to_mzml_spectra = {}
    
    available_mzml_names = [extract_filename(p) for p in mzml_paths]
    selected_set = set(selected_peptidoforms)
    
    # Pre-build unmodified set for faster lookup
    unmodified_selected = set(''.join([c for c in pep if c.isalpha()]) for pep in selected_set)
    
    total_files = len(idxml_paths)
    
    # Prepare arguments for parallel execution
    tasks = []
    for idxml_path in idxml_paths:
        tasks.append((idxml_path, available_mzml_names, selected_set, unmodified_selected))
    
    logger.info(f"Filtering {total_files} idXML files with {max_workers} workers...")
    
    processed_count = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {executor.submit(_filter_single_idxml, task): task for task in tasks}
        
        for future in as_completed(future_to_task):
            processed_count += 1
            if progress_callback:
                progress_callback(processed_count, total_files, f"Filtered {processed_count}/{total_files} files")
            
            try:
                idxml_name, mzml_spectra, file_msgs, is_empty = future.result()
                
                if is_empty:
                    empty_files.append(idxml_name)
                
                if mzml_spectra:
                    idxml_to_mzml_spectra[idxml_name] = mzml_spectra
                    messages.extend(file_msgs)
                elif file_msgs: # Error messages
                    messages.extend(file_msgs)
                    
            except Exception as e:
                logger.error(f"Error in future result: {e}")
                messages.append(f"Error: {e}")
    
    if progress_callback:
        progress_callback(total_files, total_files, "✅ Filtered all files")
    
    return idxml_to_mzml_spectra, messages, empty_files


@st.cache_resource(show_spinner=False)
def _load_indexed_mzml(mzml_path: str):
    """Cache indexed mzML file access for selective spectrum loading.
    
    Uses OnDiscMSExperiment for memory-efficient indexed access.
    Only metadata is loaded into memory, spectra loaded on-demand.
    
    Note: Uses @st.cache_resource because OnDiscMSExperiment maintains
    file handles that cannot be pickled/serialized.
    """
    try:
        # Use OnDiscMSExperiment for indexed access without loading all spectra
        indexed_file = oms.OnDiscMSExperiment()
        indexed_file.openFile(mzml_path)
        return indexed_file
    except Exception as e:
        logger.error(f"Failed to create indexed access for {mzml_path}: {e}")
        return None


@st.cache_data(show_spinner=False, ttl=3600)
def _load_idxml_peptides(idxml_path: str) -> Tuple[List, List]:
    """Cache idXML peptide identifications to avoid repeated loading.
    
    MAJOR OPTIMIZATION: idXML files were being loaded once per mzML file.
    Now loaded once per idXML file and cached for 1 hour.
    
    Returns:
        (protein_ids, peptide_ids) tuple from pyOpenMS
    """
    try:
        protein_ids = []
        peptide_ids = []
        oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)
        logger.info(f"📖 Cached idXML {Path(idxml_path).name}: {len(peptide_ids)} peptides")
        return protein_ids, peptide_ids
    except Exception as e:
        logger.error(f"Failed to load idXML {idxml_path}: {e}")
        return [], []


def _process_single_mzml_for_idxml(args: Tuple) -> Tuple[Optional[pd.DataFrame], List[str], str, str, int]:
    """Process a single mzML file for an idXML (worker function for parallel processing).
    
    Args:
        args: Tuple of (idxml_path, idxml_name, mzml_file, spec_id_peptide_pairs, mzml_path, protein_ids, peptide_ids)
              spec_id_peptide_pairs: List of (spec_id, peptide_sequence, peptide_idx) tuples if optimized,
                                    or (spec_id, peptide_sequence) if legacy.
    
    Returns:
        (mapping_df, messages, idxml_name, mzml_file, num_loaded)
    """
    idxml_path, idxml_name, mzml_file, spec_id_peptide_pairs, mzml_path, protein_ids, peptide_ids = args
    messages = []
    
    try:
        # OPTIMIZED: Use indexed access - only loads metadata, not full file
        indexed_mzml = _load_indexed_mzml(mzml_path)
        if indexed_mzml is None:
            # Check if file exists and is accessible
            if not os.path.exists(mzml_path):
                 messages.append(f"    ❌ File not found: {mzml_file}")
            elif os.path.getsize(mzml_path) == 0:
                 messages.append(f"    ❌ File is empty: {mzml_file}")
            else:
                 messages.append(f"    ❌ Corrupted or invalid mzML: {mzml_file} (Parse Error)")
            
            return None, messages, idxml_name, mzml_file, 0
            messages.append(f"    ❌ Failed to open indexed mzML: {mzml_file}")
            return None, messages, idxml_name, mzml_file, 0
        
        mappings = []
        spectrum_counter = 0
        
        # Iterate over REQUESTED spectra instead of scanning ALL peptides in the idXML
        # This reduces complexity from O(N_peptides) to O(N_requested_spectra) which is much smaller
        for item in spec_id_peptide_pairs:
            # Handle both formats (with and without peptide_idx)
            # New format: (spec_id, peptidoform, peptide_idx)
            # Old format: (spec_id, peptidoform)
            peptide_idx = None
            if len(item) >= 3:
                spectrum_id, expected_peptide, peptide_idx = item[0], item[1], item[2]
            else:
                spectrum_id, expected_peptide = item[0], item[1]
                # Fallback: scan needed (should not happen if using updated filter function)
            
            try:
                # If we have the direct index, use it!
                if peptide_idx is not None and 0 <= peptide_idx < len(peptide_ids):
                    pep_id = peptide_ids[peptide_idx]
                else:
                    # Fallback or error: skip this one (or could implement scan)
                    continue

                # Verify matches expected peptide (sanity check)
                hits = pep_id.getHits()
                if not hits:
                    continue
                
                matching_hit = None
                for hit in hits:
                    sequence = hit.getSequence()
                    peptidoform_str = sequence.toString()
                    unmodified_seq = sequence.toUnmodifiedString()
                    # Match either the full peptidoform or unmodified sequence
                    if peptidoform_str == expected_peptide or unmodified_seq == expected_peptide:
                        matching_hit = hit
                        break
                
                if matching_hit:
                    # Create mapping (spectrum_index is local to this filtered set)
                    mapping = {
                        "spectrum_index": spectrum_counter,
                        # "peptide_index": peptide_idx, # Optional
                        "spectrum_id": spectrum_id,
                        "peptide_sequence": matching_hit.getSequence().toString(),
                        "charge": matching_hit.getCharge(),
                        "run": mzml_file.replace('.mzML', ''),
                        "score": matching_hit.getScore(),
                        "mz": pep_id.getMZ(),
                        "rt": pep_id.getRT(),
                        "idxml_source": idxml_name,
                        "mzml_file": mzml_file
                    }
                    mappings.append(mapping)
                    spectrum_counter += 1
                
            except Exception as e:
                # logger.warning(f"Error processing peptide (spec={spectrum_id}): {e}")
                continue
        
        if mappings:
            df = pd.DataFrame(mappings)
            messages.append(f"    ✅ Loaded {len(mappings)} spectra")
            return df, messages, idxml_name, mzml_file, len(mappings)
        else:
            messages.append(f"    ⚠️ No spectra loaded (filtered list had {len(spec_id_peptide_pairs)} IDs)")
            return None, messages, idxml_name, mzml_file, 0
            
    except Exception as e:
        error_msg = str(e)
        if "contains no PeptideIdentifications" in error_msg or "Nothing to parse" in error_msg:
            messages.append(f"    ⚠️ Skipped {mzml_file} (empty idXML)")
        else:
            messages.append(f"    ❌ Error loading {mzml_file}: {error_msg}")
            logger.error(f"Error loading {mzml_path}: {e}")
        return None, messages, idxml_name, mzml_file, 0


def load_only_needed_spectra(idxml_paths: List[str],
                             mzml_paths: List[str],
                             idxml_to_mzml_spectra: Dict[str, Dict[str, List[Tuple[str, str]]]],
                             progress_callback=None,
                             max_workers: int = 4) -> Tuple[pd.DataFrame, Dict, List[str], List[str]]:
    """
    PHASE 3: Load ONLY the needed spectra from mzML files using indexed access.
    
    OPTIMIZED:
    - Uses OnDiscMSExperiment for selective spectrum loading
    - Parallel processing with ThreadPoolExecutor (4 workers by default)
    - Caches idXML loading (was reloading per mzML file, now once per idXML)
    - Only loads requested spectra into memory, not entire mzML files
    
    Args:
        idxml_paths: List of idXML file paths
        mzml_paths: List of mzML file paths
        idxml_to_mzml_spectra: Filtered mapping of which spectra to load
        progress_callback: Optional function(current, total, message) for progress updates
        max_workers: Number of parallel threads (default 4)
    
    Returns:
        (mapping_df, validation_info, messages, empty_files)
    """
    messages = []
    empty_files = []
    all_mapping_dfs = []
    mzml_name_to_path = {extract_filename(p): p for p in mzml_paths}
    
    # Track validation info
    validation_info = {
        "expected_files": [],
        "uploaded_files": list(mzml_name_to_path.keys()),
        "matched_files": [],
        "missing_files": [],
        "extra_files": [],
        "idxml_mappings": {}
    }
    
    # Count total mzML files to process for progress tracking
    total_mzmls = sum(len(mzml_dict) for mzml_dict in idxml_to_mzml_spectra.values())
    mzml_counter = 0
    
    # OPTIMIZED: Pre-load all idXML files once (cached)
    idxml_peptide_cache = {}
    for idxml_path in idxml_paths:
        idxml_name = extract_filename(idxml_path)
        if idxml_name in idxml_to_mzml_spectra:
            protein_ids, peptide_ids = _load_idxml_peptides(idxml_path)
            idxml_peptide_cache[idxml_name] = (protein_ids, peptide_ids, idxml_path)
            messages.append(f"📖 Cached {idxml_name}: {len(peptide_ids)} peptides")
    
    # OPTIMIZED: Build task list for parallel processing
    tasks = []
    for idxml_path in idxml_paths:
        idxml_name = extract_filename(idxml_path)
        
        # Skip if this idXML had no matching spectra during filtering
        if idxml_name not in idxml_to_mzml_spectra:
            messages.append(f"⏭️ Skipping {idxml_name} (no matching spectra)")
            continue
        
        if idxml_name not in idxml_peptide_cache:
            messages.append(f"⚠️ Skipping {idxml_name} (failed to load peptides)")
            continue
        
        protein_ids, peptide_ids, idxml_path = idxml_peptide_cache[idxml_name]
        mzml_spectra = idxml_to_mzml_spectra[idxml_name]
        messages.append(f"📊 Preparing to load spectra for {idxml_name}")
        
        # Track expected files for this idXML
        expected_mzml_files = list(mzml_spectra.keys())
        validation_info["expected_files"].extend(expected_mzml_files)
        validation_info["idxml_mappings"][idxml_name] = expected_mzml_files
        
        # Build tasks for each mzML file
        for mzml_file, spec_id_peptide_pairs in mzml_spectra.items():
            if mzml_file not in mzml_name_to_path:
                messages.append(f"  ⚠️ Missing: {mzml_file}")
                validation_info["missing_files"].append(mzml_file)
                continue
            
            mzml_path = mzml_name_to_path[mzml_file]
            task = (idxml_path, idxml_name, mzml_file, spec_id_peptide_pairs, mzml_path, protein_ids, peptide_ids)
            tasks.append(task)
    
    if not tasks:
        messages.append("⚠️ No tasks to process")
        mapping_df = pd.DataFrame()
        return mapping_df, validation_info, messages, empty_files
    
    # OPTIMIZED: Process tasks in parallel
    messages.append(f"⚡ Starting parallel processing with {max_workers} workers for {len(tasks)} mzML files")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {executor.submit(_process_single_mzml_for_idxml, task): task for task in tasks}
        
        # Process completed tasks
        for future in as_completed(future_to_task):
            mzml_counter += 1
            task = future_to_task[future]
            _, _, mzml_file, spec_id_peptide_pairs, _, _, _ = task
            
            if progress_callback:
                progress_callback(mzml_counter, total_mzmls, f"Loading mzML {mzml_counter}/{total_mzmls}: {mzml_file} ({len(spec_id_peptide_pairs)} spectra)")
            
            try:
                df, task_messages, task_idxml_name, task_mzml_file, num_loaded = future.result()
                
                # Add task messages
                messages.extend(task_messages)
                
                # Add dataframe if successful
                if df is not None and not df.empty:
                    all_mapping_dfs.append(df)
                    validation_info["matched_files"].append(task_mzml_file)
                
                # Check for empty files
                if "empty idXML" in " ".join(task_messages):
                    if task_idxml_name not in empty_files:
                        empty_files.append(task_idxml_name)
                
            except Exception as e:
                error_msg = str(e)
                messages.append(f"    ❌ Error processing task: {error_msg}")
                logger.error(f"Error in parallel task: {e}")
    
    # Combine all mappings
    if all_mapping_dfs:
        mapping_df = pd.concat(all_mapping_dfs, ignore_index=True)
    else:
        mapping_df = pd.DataFrame()
    
    # Calculate extra files
    matched_set = set(validation_info["matched_files"])
    uploaded_set = set(validation_info["uploaded_files"])
    validation_info["extra_files"] = list(uploaded_set - matched_set - set(validation_info["missing_files"]))
    
    return mapping_df, validation_info, messages, empty_files

# DEPRECATED: Old loading function - kept for compatibility but redirects to optimized version
# Remove @st.cache_data decorator to prevent memory accumulation
# @st.cache_data
def load_and_validate_files_multi(mzml_paths: List[str], idxml_path: str) -> Tuple[bool, str, Optional[pd.DataFrame], Optional[Dict]]:
    """Load and validate multiple mzML files with idXML. Returns (success, message, mapping_df, validation_info)."""
    try:
        logger.info(f"Validating files: {len(mzml_paths)} mzML files, {idxml_path}")
        
        # Extract expected mzML files from idXML
        expected_files = st_tools.extract_mzml_files_from_idxml(idxml_path)
        uploaded_files = extract_filenames(mzml_paths)
        
        # Validate file mapping
        matched_files, missing_files, extra_files = st_tools.validate_mzml_files_mapping(uploaded_files, expected_files)
        
        validation_info = create_validation_info_dict(expected_files, uploaded_files, 
                                                      matched_files, missing_files, extra_files)
        
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
        
        mapping_df = st_tools.extract_spectrum_mapping_multi(mzml_paths, idxml_path)
        if mapping_df.empty:
            return False, "No matching spectra found", None, validation_info
        
        success_msg = f"Successfully loaded {len(mapping_df)} spectrum-peptide mappings from {len(mzml_paths)} files"
        return True, success_msg, st_tools.ensure_streamlit_compatible_dataframe(mapping_df), validation_info
        
    except Exception as e:
        return False, f"Error validating files: {e}", None, None

@st.cache_data
def render_server_idxml_browser():
    """Simplified server browser - only for idXML files."""
    st.subheader("🔍 Select idXML Files from Server")
    
    # Get default directory from settings
    default_dir = st.session_state.settings.get("data_directories", {}).get("default", "/mnt/volume/")
    
    # Initialize current directory in session state
    if 'idxml_directory' not in st.session_state:
        st.session_state['idxml_directory'] = default_dir
    
    # Directory navigation
    col1, col2 = st.columns([4, 1])
    
    with col1:
        new_dir = st.text_input(
            "Directory Path", 
            value=st.session_state['idxml_directory'],
            key="idxml_dir_input",
            help="Enter the full path to directory containing idXML files"
        )
    
    with col2:
        if st.button("🔄 Refresh", key="idxml_refresh"):
            st.rerun()
    
    # Update directory if changed
    if new_dir != st.session_state['idxml_directory'] and os.path.exists(new_dir):
        st.session_state['idxml_directory'] = new_dir
        st.rerun()
    
    current_dir = st.session_state['idxml_directory']
    
    # Validate directory
    if not os.path.exists(current_dir) or not os.path.isdir(current_dir):
        st.error(f"❌ Directory not found: {current_dir}")
        return []
    
    st.success(f"📂 Current directory: {current_dir}")
    
    # Find idXML files
    idxml_files = get_files_with_substring(current_dir, "", ['.idxml', '.idXML', '.mzid', '.mzidentml'])
    
    if idxml_files:
        st.info(f"Found {len(idxml_files)} identification files")
        
        # Option to select all
        select_all = st.checkbox(f"Select all {len(idxml_files)} files", key="select_all_idxml_server")
        
        if select_all:
            return idxml_files
        else:
            # Individual file selection
            file_options = create_file_options_with_sizes(idxml_files)
            selected_display = st.multiselect(
                "Select identification files",
                options=file_options,
                key="idxml_server_select"
            )
            return map_display_names_to_paths(selected_display, idxml_files)
    else:
        st.info("No identification files found in this directory")
        return []


def render_idxml_upload_section():
    """Simplified upload section - only for idXML files."""
    st.subheader("🔍 Upload idXML Files")
    
    max_size_gb = st.session_state.settings.get("performance", {}).get("max_file_size_gb", 2)
    
    idxml_files = st.file_uploader(
        "Select identification files", 
        type=["idxml", "idXML", "mzid", "mzidentml"],
        accept_multiple_files=True, 
        help=f"Peptide identification files (up to {max_size_gb} GB each)",
        key="idxml_uploader"
    )
    
    if idxml_files:
        total_size_mb = sum(f.size for f in idxml_files) / 1024 / 1024
        st.success(f"✅ Selected {len(idxml_files)} ID files ({total_size_mb:.1f} MB)")
    
    return idxml_files or []


def process_idxml_files_only(idxml_files):
    """
    Process only idXML files to extract peptide information.
    No mzML loading or validation.
    """
    if st.button("📋 Process idXML Files", type="primary", key="process_idxml"):
        # Check if files are server paths or uploaded files
        is_server_files = idxml_files and isinstance(idxml_files[0], str)
        
        if is_server_files:
            # Server files - use paths directly
            st.info("Processing idXML files from server...")
            idxml_paths = idxml_files
        else:
            # Uploaded files - save to temp directory
            st.info("📤 Saving uploaded idXML files...")
            
            idxml_paths = process_items_with_progress(
                idxml_files,
                save_uploaded_file,
                "🔍 Saving {name}...",
                lambda x: x.name
            )
        
        # Extract peptides from idXML files only
        with st.spinner("📋 Extracting peptides from idXML files..."):
            expected_mzmls, peptidoforms_data, messages, empty_files = extract_idxml_metadata(idxml_paths)
            
            # Store in session state
            st.session_state.update({
                'idxml_processed': True,
                'idxml_paths': idxml_paths,
                'expected_mzmls': expected_mzmls,  # Which mzML files we'll need later
                'peptidoforms_data': peptidoforms_data,
                'peptides_extracted': True
            })
            
            # Show results
            if peptidoforms_data:
                st.success(f"✅ Extracted {len(peptidoforms_data)} unique peptidoforms from {len(idxml_paths)} files")
                
                # Show which mzML files will be needed
                all_expected = set()
                for mzml_list in expected_mzmls.values():
                    all_expected.update(mzml_list)
                
                st.info(f"💡 After peptide filtering, we'll need to locate {len(all_expected)} mzML files")
                
                if empty_files:
                    st.warning(f"⚠️ {len(empty_files)} files had no peptide identifications and were skipped")
                
                st.rerun()
            else:
                st.error("No peptidoforms found in the uploaded files")


def render_data_input_section_local():
    """Wrapper for render_data_input_section with spectrum_analysis specific settings.
    Shows file selection UI and confirmation button.
    """
    mzml_files, idxml_files = render_data_input_section(
        session_key_prefix="spectrum_analysis",
        show_process_button=False,  # We handle processing separately
        process_callback=None
    )
    
    # Show confirmation button only if both file types are selected
    if mzml_files and idxml_files:
        st.markdown("---")
        
        if st.button("✅ Confirm File Selection", type="primary", key="confirm_files", use_container_width=True):
            # Create progress indicators
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # IMPORTANT: Only process idXML files, not mzML files
            # mzML files are large and will be loaded later on-demand
            status_text.info(f"📋 Extracting peptides from {len(idxml_files)} idXML files...")
            
            # Extract peptides without loading mzML files
            try:
                total_files = len(idxml_files)
                
                # Wrap extraction with progress updates
                def progress_callback(current, total, message):
                    progress_bar.progress(current / total)
                    status_text.info(f"📋 Processing file {current}/{total}: {message}")
                
                expected_mzmls_by_idxml, all_peptidoforms_data, messages, empty_files = extract_idxml_metadata(
                    idxml_files, 
                    progress_callback=progress_callback
                )
                
                progress_bar.progress(1.0)
                status_text.empty()
                
                # Store paths and extracted peptides in session state
                st.session_state.update({
                    'files_confirmed': True,
                    'mzml_paths': mzml_files,
                    'idxml_paths': idxml_files,
                    'all_peptidoforms_data': all_peptidoforms_data,
                    'expected_mzmls_by_idxml': expected_mzmls_by_idxml,
                    'empty_idxml_files': empty_files,
                    'files_loaded': False  # Explicitly set to False until peptides are selected
                })
                
                st.success(f"✅ Found {len(all_peptidoforms_data)} unique peptidoforms")
                st.info("👇 **Next: Select peptides below to filter before loading spectra**")
                st.rerun()
                
            except Exception as e:
                st.error(f"❌ Error extracting peptides: {e}")
                logger.error(f"Peptide extraction failed: {e}")
                import traceback
                st.code(traceback.format_exc())
                return
    
    return mzml_files, idxml_files

# NOTE: render_unified_file_browser and render_server_directory_selector
# are now imported from data_input.py (single source of truth)

def render_peptide_filter_for_preloading():
    """
    Render peptide filtering interface BEFORE loading mzML files.
    Adapted from render_peptide_filter_and_spectrum_selector but for pre-loading stage.
    """
    st.header("📿 Peptide Filter")
    
    all_peptidoforms_data = st.session_state.get('all_peptidoforms_data', {})
    
    if not all_peptidoforms_data:
        st.error("No peptides found in idXML files. Please check your data.")
        if st.button("🔄 Start Over"):
            st.session_state.clear()
            st.rerun()
        return
    
    # Summary statistics
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("Unique Peptidoforms", len(all_peptidoforms_data))
    
    with col2:
        total_psms = sum(data['psm_count'] for data in all_peptidoforms_data.values())
        st.metric("Total PSMs", total_psms)
    
    # Filter mode selection
    filter_mode = st.radio(
        "Choose how to select peptides:",
        options=["Type/Paste peptides", "Select from Available Peptides"],
        horizontal=True,
        help="Paste: Copy-paste peptides (one per line)\nSelect: Choose from multiselect",
        key="preload_filter_mode"
    )
    
    selected_peptidoforms = []
    
    if filter_mode == "Type/Paste peptides":
        # Mode 1: Paste peptides
        pasted_text = st.text_area(
            "**Type/Paste peptide sequences** (one per line, without modifications):",
            height=150,
            placeholder="PEPTIDE\nSEQENCE\nANOTHER",
            help="Enter unmodified peptide sequences, one per line",
            key="preload_paste_peptides"
        )
        
        if pasted_text.strip():
            # Parse pasted peptides
            pasted_peptides = [line.strip().upper() for line in pasted_text.strip().split('\n') if line.strip()]
            st.info(f"📝 Pasted {len(pasted_peptides)} peptide(s)")
            
            # Match against available peptidoforms (ignoring modifications)
            matched_peptidoforms = []
            for pep, data in all_peptidoforms_data.items():
                unmodified_seq = data['sequence'].upper()
                # Check if this peptide matches any pasted sequence
                if unmodified_seq in pasted_peptides:
                    matched_peptidoforms.append(pep)
            
            selected_peptidoforms = matched_peptidoforms
            
            if not selected_peptidoforms:
                            st.warning("⚠️ No matches found for pasted peptides")
        
    else:  # Select from Available Peptides
        # Mode 2: Multiselect from available
        col1, col2 = st.columns([3, 1])
        
        with col1:
            # Search/filter box
            search_query = st.text_input(
                "🔍 Search peptides",
                placeholder="Enter sequence or part of sequence...",
                key="preload_peptide_search",
                help="Filter peptides by sequence"
            )
            
            # Filter peptides by search
            if search_query:
                filtered_peptides = [
                    pep for pep in all_peptidoforms_data.keys()
                    if search_query.upper() in pep.upper() or 
                       search_query.upper() in all_peptidoforms_data[pep]['sequence'].upper()
                ]
                st.caption(f"Showing {len(filtered_peptides)} peptides matching '{search_query}'")
            else:
                filtered_peptides = list(all_peptidoforms_data.keys())
        
        with col2:
            select_all = st.checkbox("Select All", value=False, key="preload_select_all_peptides")
        
        # Peptide selection
        if select_all:
            selected_peptidoforms = filtered_peptides
        else:
            # Create display options with PSM counts
            display_options = [
                f"{all_peptidoforms_data[pep]['display_label']} ({all_peptidoforms_data[pep]['psm_count']} PSMs)"
                for pep in filtered_peptides
            ]
            
            selected_displays = st.multiselect(
                "Select peptides to analyze",
                options=display_options,
                default=[],
                key="preload_selected_peptide_displays",
                help="Choose which peptides to include in the analysis"
            )
            
            # Map back to peptidoform strings
            selected_peptidoforms = []
            for display in selected_displays:
                # Extract peptidoform from display string
                for pep, data in all_peptidoforms_data.items():
                    expected_display = f"{data['display_label']} ({data['psm_count']} PSMs)"
                    if display == expected_display:
                        selected_peptidoforms.append(pep)
                        break
    
    # Show selection summary
    if selected_peptidoforms:
        total_psms = sum(all_peptidoforms_data[pep]['psm_count'] for pep in selected_peptidoforms)
        st.success(f"✅ Selected {len(selected_peptidoforms)} peptidoforms covering {total_psms} PSMs")
        
        # NOW process files with the selected peptides
        if st.button("🚀 Load Spectra for Selected Peptides", type="primary", key="load_spectra", use_container_width=True):
            st.session_state['selected_peptidoforms_for_loading'] = selected_peptidoforms
            
            # PHASE 2: Filter idXML by selected peptides (with progress tracking)
            idxml_paths = st.session_state['idxml_paths']
            mzml_paths = st.session_state['mzml_paths']
            st.info(f"🔍 Filtering {len(idxml_paths)} idXML file(s) by {len(selected_peptidoforms)} selected peptide(s)...")
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            idxml_to_mzml_spectra, filter_messages, filter_empty = filter_idxml_by_selected_peptides(
                idxml_paths,
                mzml_paths,
                selected_peptidoforms,
                progress_callback=lambda i, total, msg: (
                    progress_bar.progress(i / total),
                    status_text.text(msg)
                )
            )
            
            progress_bar.empty()
            status_text.empty()
            
            # PHASE 3: Load ONLY needed spectra from mzML files (with progress tracking)
            st.info(f"📊 Loading filtered spectra from mzML files...")
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            mapping_df, validation_info, load_messages, load_empty = load_only_needed_spectra(
                idxml_paths,
                st.session_state['mzml_paths'],
                idxml_to_mzml_spectra,
                progress_callback=lambda i, total, msg: (
                    progress_bar.progress(i / total),
                    status_text.text(msg)
                )
            )
            
            progress_bar.empty()
            status_text.empty()
            
            if not mapping_df.empty:
                st.success(f"✅ Loaded {len(mapping_df)} spectra for {len(selected_peptidoforms)} peptides")
                
                # Show loading details for debugging
                with st.expander("📋 Loading Details", expanded=False):
                    st.markdown("**Filter Phase:**")
                    for msg in filter_messages:
                        st.text(msg)
                    st.markdown("**Load Phase:**")
                    for msg in load_messages:
                        st.text(msg)
                
                # Store in session state
                st.session_state.update({
                    'files_loaded': True,
                    'mapping_df': mapping_df,
                    'validation_info': validation_info
                })
                
                # Show file mapping validation
                if validation_info:
                    st.markdown("---")
                    render_file_mapping_validation(validation_info)
                
                st.rerun()
            else:
                st.error("❌ No spectra loaded for selected peptides")
                # Show debug info
                with st.expander("🔍 Debug Information", expanded=True):
                    st.markdown("**Filter Messages:**")
                    for msg in filter_messages:
                        st.text(msg)
                    st.markdown("**Load Messages:**")
                    for msg in load_messages:
                        st.text(msg)
    else:
        st.info("👆 Select at least one peptide to continue")

def render_peptide_filter_and_spectrum_selector(mapping_df: pd.DataFrame) -> Optional[int]:
    """
    Render peptide filtering interface with two modes:
    1. Paste a list of peptides
    2. Select from available peptides
    Then show ALL matching spectra across all idXMLs/mzMLs.
    """
    st.header("📿 Peptide Filter")
    
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
    
    # Filter mode selection
    filter_mode = st.radio(
        "Choose how to filter peptides:",
        options=["Type peptides manually", "Select from Available Peptides"],
        horizontal=True,
        help="Paste List: Copy-paste peptides (one per line)\nSelect from Available: Choose from multiselect"
    )
    
    selected_peptides = []
    
    if filter_mode == "Type peptides manually":
        # Mode 1: Paste peptides
        pasted_text = st.text_area(
            "**Type/Paste peptide sequences** (one per line, without modifications):",
            height=150,
            placeholder="PEPTIDE\nSEQENCE\nANOTHER",
            help="Enter unmodified peptide sequences, one per line"
        )
        
        if pasted_text.strip():
            # Parse pasted peptides
            pasted_peptides = [line.strip().upper() for line in pasted_text.strip().split('\n') if line.strip()]
            st.info(f"📝 Pasted {len(pasted_peptides)} peptide(s)")
            
            # Match against available peptides (ignoring modifications)
            if "peptide_sequence" in mapping_df.columns:
                # Extract unmodified sequences from mapping_df
                mapping_df['unmodified_seq'] = mapping_df['peptide_sequence'].apply(
                    lambda x: ''.join([c for c in str(x) if c.isalpha()]).upper()
                )
                
                # Find matches
                matched_peptides = mapping_df[mapping_df['unmodified_seq'].isin(pasted_peptides)]['peptide_sequence'].unique()
                selected_peptides = list(matched_peptides)
                
                if selected_peptides:
                    st.success(f"✅ Found {len(selected_peptides)} matching peptide(s) across all modifications/charges")
                else:
                    st.warning("⚠️ No matches found for pasted peptides")
        
    else:  # Select from Available Peptides
        # Mode 2: Multiselect from available
        if "peptide_sequence" in mapping_df.columns:
            # Get unique peptides
            unique_peptides = sorted(mapping_df["peptide_sequence"].unique())
            
            # Convert to display format
            peptide_display = {st_tools.convert_to_mass_shift_format(str(p)): p for p in unique_peptides}
            
            selected_displays = st.multiselect(
                "Select peptides to analyze:",
                options=sorted(peptide_display.keys()),
                help=f"Choose from {len(unique_peptides)} unique peptides (with all modification states)"
            )
            
            if selected_displays:
                selected_peptides = [peptide_display[d] for d in selected_displays]
    
    # If peptides are selected, filter and show results
    spectrum_index = None
    
    if selected_peptides and "peptide_sequence" in mapping_df.columns:
        st.markdown("---")
        st.subheader("📊 Filtered Spectra Results")
        
        # Filter for selected peptides (match any)
        peptide_spectra = mapping_df[mapping_df["peptide_sequence"].isin(selected_peptides)].copy()
        
        if len(peptide_spectra) > 0:
            # Show summary of filtered results
            st.success(f"✅ Found **{len(peptide_spectra)}** spectra across **{len(selected_peptides)}** peptide(s)")
            
            # Show which files contain matches - check multiple possible column names
            mzml_column = None
            if "mzml_file" in peptide_spectra.columns:
                mzml_column = "mzml_file"
            elif "source_file" in peptide_spectra.columns:
                mzml_column = "source_file"
            elif "run" in peptide_spectra.columns:
                mzml_column = "run"
            
            # Create condensed display dataframe showing only essential information
            display_df = peptide_spectra.copy()
            
            # Extract only filenames (not full paths) - check multiple possible column names
            # idXML source: can be "idxml_source", "idxml_file", or extract from full path
            if "idxml_source" in display_df.columns:
                display_df["idXML"] = display_df["idxml_source"].apply(
                    lambda x: Path(x).stem if isinstance(x, str) else str(x)
                )
            elif "idxml_file" in display_df.columns:
                display_df["idXML"] = display_df["idxml_file"].apply(
                    lambda x: Path(x).stem if isinstance(x, str) else str(x)
                )
            else:
                display_df["idXML"] = "N/A"
            
            # mzML source: can be "mzml_file", "source_file", "run", or extract from full path
            if "mzml_file" in display_df.columns:
                display_df["mzML"] = display_df["mzml_file"].apply(
                    lambda x: Path(x).stem if isinstance(x, str) else str(x)
                )
            elif "source_file" in display_df.columns:
                display_df["mzML"] = display_df["source_file"].apply(
                    lambda x: Path(x).stem if isinstance(x, str) else str(x)
                )
            elif "run" in display_df.columns:
                display_df["mzML"] = display_df["run"]
            else:
                display_df["mzML"] = "N/A"
            
            # Show full peptidoform with modifications
            if "peptide_sequence" in display_df.columns:
                display_df["Peptide"] = display_df["peptide_sequence"]
            else:
                display_df["Peptide"] = "N/A"
            
            # Get spectrum ID (native ID from the mzML file)
            if "spectrum_id" in display_df.columns:
                display_df["Spectrum ID"] = display_df["spectrum_id"]
            elif "native_id" in display_df.columns:
                display_df["Spectrum ID"] = display_df["native_id"]
            elif "spectrum_index" in display_df.columns:
                display_df["Spectrum ID"] = display_df["spectrum_index"].astype(str)
            else:
                display_df["Spectrum ID"] = "N/A"
            
            # Store spectrum_index for selection (preserve original index)
            # Check for multiple possible column names
            if "spectrum_index" in display_df.columns:
                display_df["spec_idx"] = display_df["spectrum_index"]
            elif "index" in display_df.columns:
                display_df["spec_idx"] = display_df["index"]
            else:
                display_df["spec_idx"] = range(len(display_df))
            
            # Initialize selection state if not exists
            if 'selected_spectra_indices' not in st.session_state:
                st.session_state['selected_spectra_indices'] = set(display_df["spec_idx"].tolist())
            
            # Add checkbox column - all checked by default
            display_df["Select"] = display_df["spec_idx"].apply(
                lambda idx: idx in st.session_state['selected_spectra_indices']
            )
            
            # Create final condensed table with checkbox and essential columns
            condensed_columns = ["Select", "Peptide", "idXML", "mzML", "Spectrum ID"]
            
            # Show condensed filtered spectra table with editable checkboxes
            st.markdown("**Spectra found across all idXML files:**")
            st.caption("💡 Uncheck spectra to exclude them from analysis. All are selected by default.")
            
            # Controls for bulk selection
            col1, col2, col3 = st.columns([1, 1, 4])
            with col1:
                if st.button("✅ Select All", use_container_width=True):
                    st.session_state['selected_spectra_indices'] = set(display_df["spec_idx"].tolist())
                    st.rerun()
            with col2:
                if st.button("❌ Deselect All", use_container_width=True):
                    st.session_state['selected_spectra_indices'] = set()
                    st.rerun()
            
            # Display editable dataframe (include peptide_sequence for matching, but hide it)
            editor_columns = condensed_columns + ["spec_idx"]
            if "peptide_sequence" in display_df.columns:
                editor_columns.append("peptide_sequence")
            
            edited_df = st.data_editor(
                display_df[editor_columns],
                use_container_width=True,
                height=400,
                hide_index=True,
                disabled=["Peptide", "idXML", "mzML", "Spectrum ID", "spec_idx", "peptide_sequence"],
                column_config={
                    "Select": st.column_config.CheckboxColumn(
                        "Select",
                        help="Check to include spectrum in analysis",
                        default=True,
                    ),
                    "spec_idx": None,  # Hide the index column
                    "peptide_sequence": None  # Hide but preserve the full peptidoform
                },
                key="spectrum_selector"
            )
            
            # Update session state with selections
            selected_indices = edited_df[edited_df["Select"]]["spec_idx"].tolist()
            st.session_state['selected_spectra_indices'] = set(selected_indices)
            
            # Show selection summary
            num_selected = len(selected_indices)
            num_total = len(display_df)
            if num_selected > 0:
                st.success(f"✅ **{num_selected}** of **{num_total}** spectra selected for analysis")
                # Use the first selected spectrum for visualization
                spectrum_index = sorted(selected_indices)[0] if selected_indices else None
            else:
                st.warning("⚠️ No spectra selected. Please select at least one spectrum.")
                spectrum_index = None
        else:
            st.warning("⚠️ No spectra found matching the selected peptides")
            spectrum_index = None
    
    elif not selected_peptides:
        spectrum_index = None
    
    return spectrum_index

def render_spectrum_selector_only(mapping_df: pd.DataFrame):
    """
    Render spectrum selection table only (peptides already filtered during pre-loading).
    This shows the spectra that match the pre-selected peptides.
    """
    st.header("📊 Select Spectra for Analysis")
    
    if mapping_df.empty:
        st.warning("No spectra available. Please check your file loading.")
        return
    
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
        if "mzml_file" in mapping_df.columns:
            unique_files = mapping_df["mzml_file"].nunique()
        elif "source_file" in mapping_df.columns:
            unique_files = mapping_df["source_file"].nunique()
        elif "run" in mapping_df.columns:
            unique_files = mapping_df["run"].nunique()
        else:
            unique_files = 0
        st.metric("Source Files", unique_files)
    
    st.info("💡 These spectra match the peptides you selected earlier. Choose which spectra to analyze.")
    
    # Create condensed display dataframe
    display_df = mapping_df.copy()
    
    # Extract idXML source
    if "idxml_source" in display_df.columns:
        display_df["idXML"] = display_df["idxml_source"].apply(
            lambda x: Path(x).stem if isinstance(x, str) else str(x)
        )
    elif "idxml_file" in display_df.columns:
        display_df["idXML"] = display_df["idxml_file"].apply(
            lambda x: Path(x).stem if isinstance(x, str) else str(x)
        )
    else:
        display_df["idXML"] = "N/A"
    
    # Extract mzML source
    if "mzml_file" in display_df.columns:
        display_df["mzML"] = display_df["mzml_file"].apply(
            lambda x: Path(x).stem if isinstance(x, str) else str(x)
        )
    elif "source_file" in display_df.columns:
        display_df["mzML"] = display_df["source_file"].apply(
            lambda x: Path(x).stem if isinstance(x, str) else str(x)
        )
    elif "run" in display_df.columns:
        display_df["mzML"] = display_df["run"]
    else:
        display_df["mzML"] = "N/A"
    
    # Show full peptidoform with modifications
    if "peptide_sequence" in display_df.columns:
        display_df["Peptide"] = display_df["peptide_sequence"]
    else:
        display_df["Peptide"] = "N/A"
    
    # Get spectrum ID
    if "spectrum_id" in display_df.columns:
        display_df["Spectrum ID"] = display_df["spectrum_id"]
    elif "native_id" in display_df.columns:
        display_df["Spectrum ID"] = display_df["native_id"]
    elif "spectrum_index" in display_df.columns:
        display_df["Spectrum ID"] = display_df["spectrum_index"].astype(str)
    else:
        display_df["Spectrum ID"] = "N/A"
    
    # Store spectrum_index for selection
    if "spectrum_index" in display_df.columns:
        display_df["spec_idx"] = display_df["spectrum_index"]
    elif "index" in display_df.columns:
        display_df["spec_idx"] = display_df["index"]
    else:
        display_df["spec_idx"] = range(len(display_df))
    
    # Initialize selection state if not exists - select all by default
    if 'selected_spectra_indices' not in st.session_state:
        st.session_state['selected_spectra_indices'] = set(display_df["spec_idx"].tolist())
    
    # Add checkbox column
    display_df["Select"] = display_df["spec_idx"].apply(
        lambda idx: idx in st.session_state['selected_spectra_indices']
    )
    
    # Create final table with checkbox and essential columns
    condensed_columns = ["Select", "Peptide", "idXML", "mzML", "Spectrum ID"]
    
    st.markdown("**Select spectra for analysis:**")
    st.caption("💡 Uncheck spectra to exclude them from analysis. All are selected by default.")
    
    # Controls for bulk selection
    col1, col2, col3 = st.columns([1, 1, 4])
    with col1:
        if st.button("✅ Select All", use_container_width=True):
            st.session_state['selected_spectra_indices'] = set(display_df["spec_idx"].tolist())
            st.rerun()
    with col2:
        if st.button("❌ Deselect All", use_container_width=True):
            st.session_state['selected_spectra_indices'] = set()
            st.rerun()
    
    # Display editable dataframe (include peptide_sequence for matching, but hide it)
    editor_columns = condensed_columns + ["spec_idx"]
    if "peptide_sequence" in display_df.columns:
        editor_columns.append("peptide_sequence")
    
    edited_df = st.data_editor(
        display_df[editor_columns],
        use_container_width=True,
        height=400,
        hide_index=True,
        disabled=["Peptide", "idXML", "mzML", "Spectrum ID", "spec_idx", "peptide_sequence"],
        column_config={
            "Select": st.column_config.CheckboxColumn(
                "Select",
                help="Check to include spectrum in analysis",
                default=True,
            ),
            "spec_idx": None,  # Hide the index column
            "peptide_sequence": None  # Hide but preserve the full peptidoform
        },
        key="spectrum_selector_only"
    )
    
    # Update session state with selections
    selected_indices = edited_df[edited_df["Select"]]["spec_idx"].tolist()
    st.session_state['selected_spectra_indices'] = set(selected_indices)
    
    # Show selection summary
    num_selected = len(selected_indices)
    num_total = len(display_df)
    if num_selected > 0:
        st.success(f"✅ **{num_selected}** of **{num_total}** spectra selected for analysis")
    else:
        st.warning("⚠️ No spectra selected. Please select at least one spectrum.")

def display_plot_with_fallbacks(plot, debug_mode: bool) -> bool:
    """Display plot with multiple fallback methods."""
    plot_displayed = False
    
    # Method 1: Streamlit's built-in bokeh_chart
    try:
        if BOKEH_AVAILABLE:
            st.bokeh_chart(plot, use_container_width=True)
            plot_displayed = True
            if debug_mode:
                st.success("✅ Using st.bokeh_chart")
        else:
            raise ImportError("Bokeh not available")
    except ImportError:
        if debug_mode:
            st.info("Bokeh not available, trying fallback method")
    except Exception as e:
        if debug_mode:
            st.warning(f"st.bokeh_chart failed: {e}")
    
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
            if BOKEH_AVAILABLE:
                html = file_html(plot, CDN, "Spectrum Visualization")
                st.components.v1.html(html, height=750, scrolling=True)
                plot_displayed = True
                if debug_mode:
                    st.success("✅ Using HTML export")
            else:
                raise ImportError("Bokeh not available")
        except Exception as e:
            if debug_mode:
                st.error(f"HTML export failed: {e}")
    
    return plot_displayed

def render_file_mapping_validation(validation_info: Dict):
    """Render file mapping validation information for multiple idXML files."""
    st.header("🔗 File Mapping Validation")
    
    if not validation_info:
        st.info("No validation information available")
        return
    
    expected_files = validation_info.get("expected_files", [])
    uploaded_files = validation_info.get("uploaded_files", [])
    matched_files = validation_info.get("matched_files", [])
    missing_files = validation_info.get("missing_files", [])
    extra_files = validation_info.get("extra_files", [])
    idxml_mappings = validation_info.get("idxml_mappings", {})
    
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
    
    # Detailed breakdown for multiple idXML files
    if idxml_mappings:
        st.subheader("📋 idXML File Mappings")
        
        for idxml_name, expected_mzml_files in idxml_mappings.items():
            with st.expander(f"📄 {idxml_name}"):
                st.write(f"**Expected mzML files:** {len(expected_mzml_files)}")
                
                # Check which files are matched, missing, or extra for this idXML
                matched_for_this_idxml = [f for f in expected_mzml_files if f in uploaded_files]
                missing_for_this_idxml = [f for f in expected_mzml_files if f not in uploaded_files]
                
                if matched_for_this_idxml:
                    display_file_list(
                        matched_for_this_idxml,
                        "success",
                        f"✅ **{len(matched_for_this_idxml)} files matched:**"
                    )
                
                if missing_for_this_idxml:
                    display_file_list(
                        missing_for_this_idxml,
                        "error",
                        f"❌ **{len(missing_for_this_idxml)} files missing:**"
                    )
                
                if not expected_mzml_files:
                    st.warning("⚠️ No mzML file references found in this idXML")
    
    # Overall status
    display_categorized_files(
        matched_files,
        missing_files,
        extra_files,
        "files matched successfully overall:",
        "files are missing overall:",
        "unexpected files uploaded:"
    )
    
    if missing_files:
        with st.expander("💡 Tips for Missing Files"):
            st.markdown("""
            **Common solutions:**
            - Check file names match exactly (case-sensitive)
            - Ensure all files from the same experiment are uploaded
            - Verify each idXML file contains correct file references
            - Try uploading files in the correct order
            - For multiple idXML files, ensure they reference the correct mzML subsets
            """)
    
    if extra_files:
        st.info("These files may not be referenced in any identification file.")
    
    # Show expected files if validation fails
    if expected_files and (missing_files or not matched_files):
        with st.expander("📋 All Expected Files Across idXML Files"):
            st.markdown("**All identification files collectively expect these mzML files:**")
            for file in sorted(set(expected_files)):
                status = "✅" if file in matched_files else "❌"
                st.write(f"{status} {file}")
            
            if len(set(expected_files)) != len(expected_files):
                duplicates = len(expected_files) - len(set(expected_files))
                st.info(f"ℹ️ Note: {duplicates} duplicate file references found across idXML files")

# ================================================================
# MAIN WORKSPACE
# ================================================================

def main():
    """Main workspace function with optimized analysis pipeline.
    
    OPTIMIZED WORKFLOW (prevents unnecessary mzML loading):
    1. Upload/select mzML and idXML files → Confirm selection
    2. Extract peptides from idXML (fast, no mzML loading)
    3. User selects peptides to analyze → WORKFLOW HALTS HERE
    4. Filter idXML by selected peptides
    5. Load ONLY needed spectra from mzML using indexed access
        st.info("👆 Select files or folders and click 'Confirm File Selection' to continue.")
    """
    st.title("🔬 Spectrum Analysis Workspace")
    
    # STEP 1 & 2: File Input (only show if not yet confirmed)
    if not st.session_state.get('files_confirmed', False):
        render_data_input_section_local()
        st.info("👆 Select files or folders and click 'Confirm File Selection & Extract Peptides' to continue.")
        return
    
    # STEP 3 & 4: Show peptide selection after files are confirmed
    if not st.session_state.get('files_loaded', False):
        render_peptide_filter_for_preloading()
        return
    
    # STEP 5+: Only proceed if files are fully loaded with selected peptides
    if not st.session_state.get('files_loaded', False):
        st.info("👆 Please select peptides and load spectra to continue.")
        return
    
    # Step 2: Spectrum Selection (peptides already filtered during pre-loading)
    if 'mapping_df' in st.session_state and not st.session_state['mapping_df'].empty:
        st.markdown("---")
        render_spectrum_selector_only(st.session_state['mapping_df'])
        
        # Only proceed to analysis if spectra are selected
        if not st.session_state.get('selected_spectra_indices'):
            st.info("👆 Select spectra above to proceed with analysis.")
            return
    else:
        st.info("👆 Load your data files to see available spectra for selection.")
        return
    
    # Step 3: Comprehensive Analysis Settings (based on SpecViz.ipynb parameters)
    st.markdown("---")
    render_simplified_analysis_settings()


def shorten_filenames(filenames: List[str]) -> Dict[str, str]:
    """Create shortened, unique versions of filenames by removing common/repeated parts.
    
    This function identifies common prefixes, suffixes, and repeated substrings across
    all filenames and removes them to create concise, unique identifiers.
    
    Args:
        filenames: List of full filenames (e.g., mzML or idXML filenames)
        
    Returns:
        Dictionary mapping original filename to shortened version
    """
    if not filenames or len(filenames) == 0:
        return {}
    
    if len(filenames) == 1:
        fname = filenames[0]
        stem = Path(fname).stem if '.' in fname else fname
        return {fname: stem}
    
    # Remove file extensions
    stems = [Path(fname).stem if '.' in fname else fname for fname in filenames]
    
    # Find common prefix
    def find_common_prefix(strings):
        if not strings:
            return ""
        prefix = strings[0]
        for s in strings[1:]:
            while not s.startswith(prefix):
                prefix = prefix[:-1]
                if not prefix:
                    return ""
        return prefix
    
    # Find common suffix
    def find_common_suffix(strings):
        if not strings:
            return ""
        reversed_strings = [s[::-1] for s in strings]
        common_prefix = find_common_prefix(reversed_strings)
        return common_prefix[::-1]
    
    # Remove common prefix
    common_prefix = find_common_prefix(stems)
    if len(common_prefix) > 5:
        stems = [s[len(common_prefix):] for s in stems]
    
    # Remove common suffix
    common_suffix = find_common_suffix(stems)
    if len(common_suffix) > 5:
        stems = [s[:-len(common_suffix)] for s in stems]
    
    # Clean up underscores
    stems = [s.strip('_').replace('__', '_') for s in stems]
    
    # Check uniqueness
    if len(set(stems)) != len(stems):
        stems = [Path(fname).stem if '.' in fname else fname for fname in filenames]
    
    # Limit length if too long
    max_length = 50
    final_stems = []
    for stem in stems:
        if len(stem) > max_length:
            shortened = stem[:20] + "..." + stem[-25:]
            final_stems.append(shortened)
        else:
            final_stems.append(stem)
    
    return {fname: short for fname, short in zip(filenames, final_stems)}
    
    # Step 4: Automatic Analysis and Visualization
    if st.session_state.get('analysis_ready', False):
        render_automated_analysis()

def render_simplified_analysis_settings():
    """Simplified analysis settings based on SpecViz.ipynb parameters."""
    st.header("⚙️ Analysis Settings")
    
    # Early exit if analysis is already ready - no need to recompute everything
    if st.session_state.get('analysis_ready', False):
        st.success("✅ Analysis configured and ready. See results below.")
        if st.button("🔄 Reconfigure Analysis", use_container_width=True):
            st.session_state['analysis_ready'] = False
            st.rerun()
        return
    
    # Get pre-selected spectra from peptide filter step
    selected_spectra_indices = st.session_state.get('selected_spectra_indices', set())
    if not selected_spectra_indices:
        st.error("No spectra selected. Please select spectra in the previous step.")
        return
    
    mapping_df = st.session_state.get('mapping_df', pd.DataFrame())
    if mapping_df.empty:
        st.error("No PSM data available. Please check your files.")
        return
    
    # Extract peptides from the selected spectra (no double filtering)
    try:
        # Get the selected spectra rows from mapping_df
        if "spectrum_index" in mapping_df.columns:
            selected_rows = mapping_df[mapping_df["spectrum_index"].isin(selected_spectra_indices)]
        elif "index" in mapping_df.columns:
            selected_rows = mapping_df[mapping_df["index"].isin(selected_spectra_indices)]
        else:
            st.error("Could not find spectrum index column in mapping data.")
            return
        
        if selected_rows.empty:
            st.error("No matching spectra found for the selected indices.")
            return
        
        # Extract unique peptide sequences from selected spectra
        if "peptide_sequence" not in selected_rows.columns:
            st.error("No peptide sequence information available.")
            return
        
        selected_peptides = selected_rows["peptide_sequence"].unique().tolist()
        
        # Display info about pre-selected peptides and spectra
        st.info(f"📊 Using **{len(selected_spectra_indices)}** selected spectra covering **{len(selected_peptides)}** unique peptide(s) from the previous step.")
        
        # Build a set of selected spectrum IDs from the mapping_df for fast lookup
        # This ensures we only analyze PSMs that correspond to user-selected spectra
        selected_spectrum_ids = set()
        if "spectrum_id" in selected_rows.columns:
            selected_spectrum_ids = set(selected_rows["spectrum_id"].dropna().unique())
        
        # Also track the selected spectrum_index values for additional filtering
        selected_spec_indices_list = list(selected_spectra_indices)
        
        logger.info(f"Filtering PSMs: {len(selected_spectrum_ids)} unique spectrum IDs from {len(selected_spectra_indices)} selected spectra")
        
        # Check if PSMs are already cached in session state to avoid expensive re-loading
        cache_key = f"psm_cache_{hash(tuple(sorted(st.session_state.get('idxml_paths', []))))}"
        
        if cache_key in st.session_state:
            # Use cached PSM lists
            all_psm_lists = st.session_state[cache_key]
        else:
            # Load PSM data from idXML files (expensive operation - cache it!)
            idxml_paths = st.session_state.get('idxml_paths')
            if not idxml_paths:
                st.error("No idXML file paths found in session state.")
                return
            
            all_psm_lists = []
            with st.spinner(f"Loading PSM data from {len(idxml_paths)} idXML files..."):
                for idxml_path in idxml_paths:
                    try:
                        psm_list = read_file(idxml_path, file_format="idXML")
                        all_psm_lists.append((psm_list, Path(idxml_path).name))
                    except Exception as e:
                        logger.error(f"Error loading {idxml_path}: {e}")
                        st.warning(f"⚠️ Could not load {Path(idxml_path).name}")
            
            # Cache the loaded PSMs
            st.session_state[cache_key] = all_psm_lists
        
        # Create cache key for filtered PSMs based on selected spectra and peptides
        filtered_cache_key = f"filtered_psm_cache_{hash((tuple(sorted(selected_spectrum_ids)), tuple(sorted(selected_peptides)), cache_key))}"
        
        # Check if filtered results are already cached
        if filtered_cache_key in st.session_state:
            # Use cached filtered PSMs and peptidoforms data
            all_psm_lists_filtered = st.session_state[filtered_cache_key]['filtered_psms']
            unique_peptidoforms_data = st.session_state[filtered_cache_key]['peptidoforms_data']
        else:
            # Need to filter PSMs - do it once and cache the results
            all_psm_lists_filtered = []  # Store filtered PSMs only
            unique_peptidoforms_data = {}
            
            # OPTIMIZATION: Pre-build lookup sets for fast matching
            selected_peptides_set = set(selected_peptides)
            selected_peptides_unmod = set(''.join([c for c in str(p) if c.isalpha()]).upper() for p in selected_peptides)
            
            for psm_list, idxml_name in all_psm_lists:
                try:
                    # OPTIMIZED: Filter PSMs using fast set lookups
                    filtered_psm_list = []
                    
                    for psm in psm_list:
                        try:
                            peptidoform_str = str(psm.peptidoform)
                            unmodified_seq = Peptidoform(peptidoform_str).sequence
                            
                            # OPTIMIZED: Quick spectrum ID check using set lookup
                            if hasattr(psm, 'spectrum_id') and psm.spectrum_id and psm.spectrum_id not in selected_spectrum_ids:
                                continue
                            
                            # OPTIMIZED: Quick peptide check using set lookup
                            if peptidoform_str not in selected_peptides_set and unmodified_seq not in selected_peptides_unmod:
                                continue
                            
                            # Both filters passed - include this PSM
                            filtered_psm_list.append(psm)
                            
                            # Initialize or update peptidoform data
                            if peptidoform_str not in unique_peptidoforms_data:
                                # Create display label showing modifications if present
                                if peptidoform_str != unmodified_seq:
                                    display_label = f"{unmodified_seq} [{peptidoform_str}]"
                                else:
                                    display_label = unmodified_seq
                                
                                unique_peptidoforms_data[peptidoform_str] = {
                                    'peptidoform': peptidoform_str,
                                    'sequence': unmodified_seq,
                                    'display_label': display_label,
                                    'psm_count': 0
                                }
                            
                            # Increment PSM count for this peptidoform
                            unique_peptidoforms_data[peptidoform_str]['psm_count'] += 1
                            
                        except Exception as e:
                            logger.warning(f"Could not extract peptidoform from {psm.peptidoform}: {e}")
                            continue
                    
                    # Store the filtered PSM list
                    all_psm_lists_filtered.append((filtered_psm_list, idxml_name))
                except Exception as e:
                    logger.error(f"Error loading {idxml_name}: {e}")
                    st.error(f"Error loading {idxml_name}: {e}")
                    continue
            
            # Cache the filtered results to avoid recomputation
            st.session_state[filtered_cache_key] = {
                'filtered_psms': all_psm_lists_filtered,
                'peptidoforms_data': unique_peptidoforms_data
            }
        
        # Sort peptidoforms by PSM count (descending), then by unmodified sequence
        unique_peptidoforms_list = sorted(list(unique_peptidoforms_data.values()), 
                                        key=lambda x: (-x['psm_count'], x['sequence'], x['peptidoform']))
        
        if not unique_peptidoforms_list:
            st.error("No valid peptidoforms found in the selected spectra.")
            return
        
        # Extract peptidoforms and sequences for use in analysis
        selected_peptidoforms = [item['peptidoform'] for item in unique_peptidoforms_list]
        selected_sequences = [item['sequence'] for item in unique_peptidoforms_list]
        
        # Analysis summary
        total_filtered_psms = sum(len(filtered_list) for filtered_list, _ in all_psm_lists_filtered)
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Selected Spectra", len(selected_spectra_indices))
        with col2:
            st.metric("Unique Peptidoforms", len(unique_peptidoforms_list))
        with col3:
            st.metric("Total PSMs", total_filtered_psms)
        with col4:
            st.metric("Selected Peptides", len(selected_peptides))
        
        # ==========================================
        # SECTION 1: COMMON SETTINGS IN ONE ROW
        # ==========================================
        st.caption("**Settings that apply to all analysis modes:**")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # R² Metric Selection
            r2_metric_options = {
                "Both": "both",
                "Matching Ions Only": "matching",
                "Binned Intensities Only": "binned"
            }
            r2_metric_display = st.selectbox(
                "📊 R² Metric Display",
                options=list(r2_metric_options.keys()),
                index=0,  # Default to "Both"
                help="Choose which R² correlation metric to display in spectrum comparisons:\n"
                     "• Both: Show both matching ions and binned R²\n"
                     "• Matching Ions: Only show R² calculated from annotated matching ions\n"
                     "• Binned Intensities: Only show R² from binned spectrum correlation",
                key="analysis_r2_metric_select"
            )
            r2_metric = r2_metric_options[r2_metric_display]
            
            # Store R² metric preference in session state
            st.session_state['r2_metric'] = r2_metric
        
        with col2:
            # Fragment mass tolerance (matches provided code convention)
            fragment_mass_tolerance = st.number_input(
                "📏 Fragment Mass Tolerance (Da)",
                min_value=0.001,
                max_value=1.0,
                value=st.session_state.settings.get("spectrum_analysis", {}).get("default_fragment_mass_tolerance", 0.02),  
                step=0.001,
                help="The fragment mass tolerance of theoretical and experimental spectrum matching",
                format="%.3f",
                key="analysis_fragment_tolerance"
            )
        
        with col3:
            # Maximum fragment charge
            max_frag_charge = st.number_input(
                "⚡ Max Fragment Charge",
                min_value=1,
                max_value=5,
                value=st.session_state.settings.get("spectrum_analysis", {}).get("default_max_fragment_charge", 2),
                help="Maximum charge state for fragment ions",
                key="analysis_max_frag_charge"
            )
        
        # Ion type selection - split into two rows for all ion types
        st.write("🧪 **Ion Types for Annotation:**")
        st.caption("Select fragment ion types to include in the theoretical spectrum")
        
        # First row: Standard fragment ions
        col_a, col_b, col_c, col_x, col_y, col_z = st.columns(6)
        
        with col_a:
            include_a_ions = st.checkbox("a", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("a", False), help="a-type fragment ions", key="analysis_ion_a")
        with col_b:
            include_b_ions = st.checkbox("b", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("b", True), help="b-type fragment ions", key="analysis_ion_b")
        with col_c:
            include_c_ions = st.checkbox("c", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("c", False), help="c-type fragment ions", key="analysis_ion_c")
        with col_x:
            include_x_ions = st.checkbox("x", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("x", False), help="x-type fragment ions", key="analysis_ion_x")
        with col_y:
            include_y_ions = st.checkbox("y", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("y", True), help="y-type fragment ions", key="analysis_ion_y")
        with col_z:
            include_z_ions = st.checkbox("z", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("z", False), help="z-type fragment ions", key="analysis_ion_z")
        
        # Second row: Additional ion types
        col_loss, col_precursor, col_prefix, col_immonium = st.columns(4)
        
        with col_loss:
            include_loss_ions = st.checkbox("loss", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("loss", False), help="Neutral loss ions (e.g., b-H2O, y-NH3)", key="analysis_ion_loss")
        with col_precursor:
            include_precursor_ions = st.checkbox("precursor", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("precursor", False), help="Precursor ions", key="analysis_ion_precursor")
        with col_prefix:
            include_prefix_ions = st.checkbox("prefix", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("prefix", False), help="First prefix ions", key="analysis_ion_prefix")
        with col_immonium:
            include_immonium_ions = st.checkbox("immonium", value=st.session_state.settings.get("spectrum_analysis", {}).get("default_ion_types", {}).get("immonium", False), help="Abundant immonium ions", key="analysis_ion_immonium")
        
        # Create list of included ions
        included_ions = []
        ion_selections = {
            'a': include_a_ions,
            'b': include_b_ions,
            'c': include_c_ions,
            'x': include_x_ions,
            'y': include_y_ions,
            'z': include_z_ions,
            'loss': include_loss_ions,
            'precursor': include_precursor_ions,
            'prefix': include_prefix_ions,
            'immonium': include_immonium_ions
        }
        
        for ion_type, is_selected in ion_selections.items():
            if is_selected:
                included_ions.append(ion_type)
        
        if not included_ions:
            st.warning("⚠️ Please select at least one ion type for analysis.")
            return
        
        st.info(f"Selected ion types: {', '.join(included_ions)}")
        
        st.markdown("---")  # Divider
        
        # ==========================================
        # SECTION 2: ANALYSIS MODE SELECTION
        # ==========================================
        st.subheader("🔬 Analysis Mode Selection")
        st.caption("Choose one or more analysis modes to generate different types of mirror plots")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            enable_predicted_mode = st.checkbox(
                "✨ Predicted Spectrum",
                value=False,
                help="Compare experimental spectra against AI-predicted spectra (Koina models)\n"
                     "Generates: Exp vs Predicted mirror plots for each PSM",
                key="analysis_mode_predicted"
            )
        
        with col2:
            enable_all_vs_all_mode = st.checkbox(
                "🔄 All-vs-All Replicates", 
                value=False,
                help="Compare all replicates of the same peptide against each other\n"
                     "Generates: All pairwise comparisons within each peptidoform",
                key="analysis_mode_all_vs_all"
            )
        
        with col3:
            enable_one_vs_all_mode = st.checkbox(
                "📍 One-vs-All Replicates",
                value=True,
                help="Select one reference spectrum and compare all other replicates against it\n"
                     "Generates: Reference vs each replicate comparison",
                key="analysis_mode_one_vs_all"
            )
        
        # Check if at least one mode is selected
        if not (enable_predicted_mode or enable_all_vs_all_mode or enable_one_vs_all_mode):
            st.warning("⚠️ Please select at least one analysis mode.")
            return
        
        # Mode-specific settings
        selected_model = None
        collision_energy = None
        instrument_type = None
        reference_runs = []
        
        # Settings for Predicted Spectrum mode
        if enable_predicted_mode:
            st.write("**⚙️ Predicted Spectrum Mode Settings:**")
            col_pred1, col_pred2, col_pred3 = st.columns(3)
            
            with col_pred1:
                # Model selection
                if st_tools and hasattr(st_tools, 'get_available_models'):
                    try:
                        available_models = st_tools.get_available_models()
                    except Exception as e:
                        logger.error(f"Analysis settings error: {e}")
                        st.warning(f"⚠️ Could not load prediction models: {e}")
                        available_models = {"Default": "default"}
                else:
                    available_models = {"Default": "default"}
                
                selected_model = st.selectbox(
                    "🤖 Prediction Model",
                    list(available_models.keys()),
                    index=0,
                    help="AI model for fragment ion prediction",
                    key="analysis_pred_model_select"
                )
                selected_model = available_models[selected_model]  # Get model ID
            
            with col_pred2:
                # Collision Energy
                collision_energy = st.number_input(
                    "💥 Collision Energy",
                    min_value=1,
                    max_value=100,
                    value=st.session_state.settings.get("spectrum_analysis", {}).get("default_collision_energy", 35),
                    help="Collision energy for fragmentation",
                    key="analysis_pred_collision_energy"
                )
            
            with col_pred3:
                # Instrument Type
                instrument_options = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]
                default_instrument = st.session_state.settings.get("spectrum_analysis", {}).get("default_instrument_type", "TIMSTOF")
                default_index = instrument_options.index(default_instrument) if default_instrument in instrument_options else 0
                
                instrument_type = st.selectbox(
                    "🔬 Instrument Type",
                    options=instrument_options,
                    index=default_index,
                    help="Mass spectrometer type",
                    key="analysis_pred_instrument_type"
                )
        
        # Settings for All-vs-All mode
        group_by_unmodified = False
        if enable_all_vs_all_mode:
            st.write("**⚙️ All-vs-All Mode Settings:**")
            group_by_unmodified = st.checkbox(
                "🔗 Group different modifications of the same peptide together",
                value=False,
                help="When enabled, peptides with different modifications (e.g., oxidation, phosphorylation) "
                     "will be compared against each other as if they were the same peptide.\n"
                     "When disabled, only PSMs with identical modifications will be compared.",
                key="analysis_all_vs_all_group_mods"
            )
            
            if group_by_unmodified:
                st.caption("✅ Will compare all PSMs with the same unmodified sequence (ignoring modifications)")
            else:
                st.caption("ℹ️ Will only compare PSMs with identical peptidoform (same modifications)")
        
        # Settings for One-vs-All mode
        if enable_one_vs_all_mode:
            st.write("**⚙️ One-vs-All Mode Settings:**")
            
            # Get available mzML files/runs from the uploaded data
            mzml_paths = st.session_state.get('mzml_paths', [])
            if mzml_paths:
                # Extract run names (stems) from mzML files
                available_runs = [Path(mzml).stem for mzml in mzml_paths]
                
                reference_runs = st.multiselect(
                    "📍 Select Reference Spectra (Standards)",
                    options=available_runs,
                    help="Choose one or more reference mzML files. All other replicates will be compared against these standards.",
                    key="analysis_one_vs_all_reference"
                )
                
                if reference_runs:
                    if len(reference_runs) == 1:
                        st.info(f"Reference: **{reference_runs[0]}** (will appear as lower/mirrored spectrum in plots)")
                    else:
                        st.info(f"References: **{', '.join(reference_runs)}** (will appear as lower/mirrored spectra in plots)")
                else:
                    st.warning("Please select at least one reference spectrum to enable One-vs-All mode.")
                
                # Output directory for exported plots
                col_dir1, col_dir2 = st.columns([3, 1])
                with col_dir1:
                    plot_output_dir = st.text_input(
                        "📁 Plot Export Directory",
                        value=st.session_state.get('plot_output_dir', './output_plots'),
                        help="Directory where PNG and SVG files will be saved",
                        key="analysis_plot_output_dir_input"
                    )
                    st.session_state['plot_output_dir'] = plot_output_dir
                
                with col_dir2:
                    st.write("")  # Spacing
                    st.write("")  # Spacing
                    if st.button("📂 Browse", help="Open file browser"):
                        st.info(f"Current: {plot_output_dir}")
                
                st.caption(f"Plots will be saved as PNG and SVG to: `{plot_output_dir}`")
            else:
                st.error("No mzML files found. Please load data first.")
                return
        
        st.markdown("---")  # Divider
        
        # Show peptide statistics using the counts already calculated
        total_peptides = len(unique_peptidoforms_list)
        total_psms = sum(item['psm_count'] for item in unique_peptidoforms_list)
        st.info(f"📊 Analyzing **{total_peptides} peptidoforms** with **{total_psms} total PSMs** from selected spectra across {len(all_psm_lists)} idXML files")
        
        # Show per-peptide counts (already visible in selection, but detailed breakdown)
        with st.expander("📋 Detailed PSM Breakdown by File (Filtered)", expanded=False):
            # Only compute breakdown when expander is opened (lazy loading)
            # Build a lookup dict once for performance
            peptidoform_to_display = {item['peptidoform']: item['display_label'] 
                                     for item in unique_peptidoforms_list}
            
            for filtered_psm_list, idxml_name in all_psm_lists_filtered:
                file_total = len(filtered_psm_list)
                
                if file_total > 0:
                    # Quick count without iterating through display labels multiple times
                    file_peptide_matches = {}
                    for psm in filtered_psm_list:
                        try:
                            psm_peptidoform = str(psm.peptidoform)
                            file_peptide_matches[psm_peptidoform] = file_peptide_matches.get(psm_peptidoform, 0) + 1
                        except:
                            continue
                    
                    st.write(f"**{idxml_name}**: {file_total} PSMs (from selected spectra)")
                    for peptidoform, count in file_peptide_matches.items():
                        display_label = peptidoform_to_display.get(peptidoform, peptidoform)
                        st.write(f"    • {display_label}: {count} PSMs")
                else:
                    st.write(f"**{idxml_name}**: No matching PSMs")
        
        # Store FILTERED PSM lists for later use (only PSMs matching selected spectra and peptides)
        st.session_state['all_psm_lists'] = all_psm_lists_filtered
        
        # Build PSM to idXML mapping for export functionality
        psm_to_idxml_map = {}
        for psm_list, idxml_name in all_psm_lists_filtered:
            for psm in psm_list:
                psm_key = f"{psm.spectrum_id}_{psm.run}"
                psm_to_idxml_map[psm_key] = idxml_name
        st.session_state['psm_to_idxml_map'] = psm_to_idxml_map
        
        # Start analysis button
        if st.button("🚀 Start Analysis", type="primary", use_container_width=True):
            # Store settings directly (avoid intermediate dict creation)
            st.session_state['selected_peptides'] = selected_sequences
            st.session_state['selected_peptidoforms'] = selected_peptidoforms
            st.session_state['selected_peptides_display'] = [item['display_label'] for item in unique_peptidoforms_list]
            st.session_state['selected_peptidoforms_data'] = unique_peptidoforms_list
            st.session_state['selected_model'] = selected_model if enable_predicted_mode else None
            st.session_state['max_frag_charge'] = max_frag_charge
            st.session_state['fragment_mass_tolerance'] = fragment_mass_tolerance
            st.session_state['collision_energy'] = collision_energy if enable_predicted_mode else 35
            st.session_state['instrument_type'] = instrument_type if enable_predicted_mode else "TIMSTOF"
            st.session_state['included_ions'] = included_ions
            st.session_state['enable_predicted_mode'] = enable_predicted_mode
            st.session_state['enable_all_vs_all_mode'] = enable_all_vs_all_mode
            st.session_state['enable_one_vs_all_mode'] = enable_one_vs_all_mode
            st.session_state['reference_runs'] = reference_runs if enable_one_vs_all_mode else []
            st.session_state['group_by_unmodified'] = group_by_unmodified
            st.session_state['analysis_ready'] = True
            st.session_state['selected_spectrum_ids'] = selected_spectrum_ids
            st.session_state['num_selected_spectra'] = len(selected_spectra_indices)
            st.rerun()
    
    except Exception as e:
        st.error(f"Error in analysis settings: {e}")
        logger.error(f"Analysis settings error: {e}", exc_info=True)

def render_automated_analysis():
    """Automated analysis following SpecViz.ipynb workflow."""
    st.header("📈 Analysis Results")
    
    try:
        # Get parameters from session state - now supporting multiple peptides
        selected_peptides = st.session_state['selected_peptides']  # Unmodified sequences for compatibility
        selected_peptidoforms = st.session_state['selected_peptidoforms']  # Full peptidoforms with modifications
        selected_peptides_display = st.session_state['selected_peptides_display']  # Display labels
        selected_peptidoforms_data = st.session_state['selected_peptidoforms_data']  # Complete data
        selected_model = st.session_state['selected_model']
        max_frag_charge = st.session_state['max_frag_charge']
        fragment_mass_tolerance = st.session_state['fragment_mass_tolerance']  # Updated variable name
        collision_energy = st.session_state['collision_energy']  # New parameter
        instrument_type = st.session_state['instrument_type']  # New parameter
        included_ions = st.session_state['included_ions']  # New parameter for ion types
        
        # Create cache key for analysis results (excludes export directory to prevent recomputation)
        analysis_cache_key = (
            tuple(sorted(selected_peptidoforms)),
            selected_model,
            max_frag_charge,
            fragment_mass_tolerance,
            collision_energy,
            instrument_type,
            tuple(sorted(included_ions)),
            st.session_state.get('enable_predicted_mode', False),
            st.session_state.get('enable_all_vs_all_mode', False),
            st.session_state.get('enable_one_vs_all_mode', False),
            tuple(sorted(st.session_state.get('reference_runs', []))),
            st.session_state.get('num_selected_spectra', 0)
        )
        
        # Check if analysis results are already cached
        if 'analysis_cache_key' in st.session_state and st.session_state['analysis_cache_key'] == analysis_cache_key:
            # Use cached results
            if 'cached_analysis_results' in st.session_state:
                logger.info("Using cached analysis results")
                cached_results = st.session_state['cached_analysis_results']
                # Restore cached variables
                aligned_df = cached_results['aligned_df']
                final_psms = cached_results['final_psms']
                psm_to_idxml_map = cached_results['psm_to_idxml_map']
                exp_pred_plots = cached_results.get('exp_pred_plots', [])
                exp_exp_plots = cached_results.get('exp_exp_plots', [])
                one_vs_all_plots = cached_results.get('one_vs_all_plots', [])
                exp_pred_export_data = cached_results.get('exp_pred_export_data', [])
                exp_exp_export_data = cached_results.get('exp_exp_export_data', [])
                one_vs_all_export_data = cached_results.get('one_vs_all_export_data', [])
                
                # Skip to display section
                use_cached = True
            else:
                use_cached = False
        else:
            use_cached = False
        
        if not use_cached:
            # Get file paths - now supporting multiple idXML files
            if 'mzml_paths' in st.session_state:
                mzml_paths = st.session_state['mzml_paths']
            else:
                mzml_paths = [st.session_state['mzml_path']]
            idxml_paths = st.session_state['idxml_paths']
            
            # Verification: Show what PSMs are being analyzed
            all_psm_lists = st.session_state.get('all_psm_lists', [])
            num_selected_spectra = st.session_state.get('num_selected_spectra', 0)
            
            # Count total PSMs across all files
            total_psms_to_analyze = sum(len(psm_list) for psm_list, _ in all_psm_lists)
            
            st.success(f"✅ Analysis ready: **{total_psms_to_analyze} PSMs** from **{num_selected_spectra} selected spectra** across **{len(selected_peptidoforms)} peptidoforms**")
            
            # Show detailed breakdown
            with st.expander("🔍 Verification: PSMs included in analysis"):
                st.caption("These PSMs passed both filters: selected spectra AND selected peptides")
                for psm_list, idxml_name in all_psm_lists:
                    st.write(f"**{idxml_name}**: {len(psm_list)} PSMs")
                    # Group by peptidoform
                    peptide_counts = {}
                    for psm in psm_list:
                        peptidoform_str = str(psm.peptidoform)
                        peptide_counts[peptidoform_str] = peptide_counts.get(peptidoform_str, 0) + 1
                    
                    for peptidoform, count in sorted(peptide_counts.items()):
                        display_label = next((item['display_label'] for item in selected_peptidoforms_data 
                                            if item['peptidoform'] == peptidoform), peptidoform)
                        st.write(f"    • {display_label}: {count} PSMs")
        else:
            # Using cached results - skip heavy computation
            st.info(f"⚡ Using cached analysis results for {len(selected_peptidoforms)} peptidoforms")
        
        if not use_cached:
            with st.spinner("🔄 Running analysis workflow..."):
                # Use the pre-filtered PSM lists from analysis settings
                # These have already been filtered by selected spectra AND selected peptides
                psms_with_spectra = []
                for psm_list, idxml_name in all_psm_lists:
                    psms_with_spectra.extend(psm_list)
                
                if len(psms_with_spectra) == 0:
                    st.error("No PSMs available for analysis. Please check your selection.")
                    logger.error("PIPELINE STOPPED: No PSMs available")
                    return
                
                # STEP 1: Load experiments and create lookups (optimized)
                experiments = create_experiment_map(mzml_paths, show_ui=False)
                lookup_map = create_spectrum_lookups(experiments)
                
                # STEP 2: Attach spectra to PSMs (optimized)
                exp_stem_to_name, exp_base_to_name = create_fuzzy_match_helpers(experiments)
                
                # OPTIMIZATION: Only process PSMs that don't have spectra yet, filter inline
                final_psms = []
                for psm in psms_with_spectra:
                    # Skip if already has spectrum
                    if hasattr(psm, "spectrum") and psm.spectrum is not None:
                        final_psms.append(psm)
                        continue
                    
                    # Try to attach spectrum
                    exp_name, lookup_key = find_matching_experiment(
                        psm.run, experiments, exp_stem_to_name, exp_base_to_name
                    )
                    
                    if exp_name and lookup_key and lookup_key in lookup_map:
                        try:
                            idx = lookup_map[lookup_key].findByNativeID(psm.spectrum_id)
                            psm.spectrum = experiments[exp_name][idx]
                            final_psms.append(psm)
                        except:
                            pass
                
                if len(final_psms) == 0:
                    st.error("No experimental spectra could be loaded for the selected PSMs.")
                    return
                
                # STEP 3: Generate alignments (optimized - reuse generator)
                if not final_psms:
                    st.error("No experimental spectra could be loaded for the selected PSMs.")
                    return
                
                aligned_rows = []
                
                # OPTIMIZATION: Create generator once with settings, reuse for all PSMs
                gen = oms.TheoreticalSpectrumGenerator()
                params = gen.getParameters()
                for ion in ["a", "b", "y", "c", "x", "z"]: 
                    params.setValue(f"add_{ion}_ions", "true" if ion in included_ions else "false")
                params.setValue("add_metainfo", "true")
                params.setValue("add_losses", "true" if "loss" in included_ions else "false")
                params.setValue("add_precursor_peaks", "true" if "precursor" in included_ions else "false")
                params.setValue("add_first_prefix_ion", "true" if "prefix" in included_ions else "false")
                params.setValue("add_abundant_immonium_ions", "true" if "immonium" in included_ions else "false")
                gen.setParameters(params)
                
                # OPTIMIZATION: Create aligner once, reuse for all PSMs
                aligner = oms.SpectrumAlignment()
                prm = aligner.getParameters()
                prm.setValue("tolerance", fragment_mass_tolerance)
                prm.setValue("is_relative_tolerance", "false")
                aligner.setParameters(prm)
                
                for psm in final_psms:
                    try:
                        # Generate theoretical spectrum (reuse generator)
                        theo_spec = oms.MSSpectrum()
                        gen.getSpectrum(theo_spec, oms.AASequence(str(psm.peptidoform).split('/')[0].replace('[','(').replace(']',')')), 
                                        1, max_frag_charge)
                        theo_spec.sortByPosition()
                        
                        sda, ida = theo_spec.getStringDataArrays(), theo_spec.getIntegerDataArrays()
                        names = [n.decode() for n in sda[0]] if sda and sda[0].size() > 0 else [""] * theo_spec.size()
                        charges = list(ida[0]) if ida and ida[0].size() > 0 else [None] * theo_spec.size()
                        
                        # Align spectra (reuse aligner)
                        psm.spectrum.sortByPosition()
                        pairs = []
                        aligner.getSpectrumAlignment(pairs, theo_spec, psm.spectrum)
                        
                        # Create aligned dataframe
                        rows = []
                        for theo_idx, exp_idx in pairs:
                            rows.append({
                                'exp_idx': exp_idx,
                                'theo_idx': theo_idx,
                                'mz': theo_spec[theo_idx].getMZ(),
                                'intensity': theo_spec[theo_idx].getIntensity(),
                                    'annotation': names[theo_idx] if theo_idx < len(names) else "",
                                    'charge': charges[theo_idx] if theo_idx < len(charges) else None,
                                    'exp_mz': psm.spectrum[exp_idx].getMZ(),
                                    'exp_intensity': psm.spectrum[exp_idx].getIntensity()
                                })
                        
                        if rows:
                            aligned_df_psm = pd.DataFrame(rows)
                            aligned_df_psm["psm_id"] = psm.spectrum_id
                            aligned_rows.append(aligned_df_psm.assign(run=psm.run))
                    except Exception as e:
                        logger.error(f"Error processing PSM {psm.spectrum_id}: {e}")
                        continue
                
                if not aligned_rows:
                    st.error("No successful alignments could be generated.")
                    return
                
                aligned_df = pd.concat(aligned_rows, ignore_index=True, sort=False)
                st.success(f"✅ Generated {len(aligned_df)} peak alignments from {len(aligned_rows)} PSMs")
        
        # Get analysis mode flags
        enable_predicted_mode = st.session_state.get('enable_predicted_mode', True)
        enable_all_vs_all_mode = st.session_state.get('enable_all_vs_all_mode', True)
        enable_one_vs_all_mode = st.session_state.get('enable_one_vs_all_mode', False)
        reference_runs = st.session_state.get('reference_runs', [])
        
        # Get the PSM to idXML mapping
        psm_to_idxml_map = st.session_state.get('psm_to_idxml_map', {})
        
        # Initialize data collections
        exp_pred_plots = []
        exp_exp_plots = []
        one_vs_all_plots = []
        
        # Initialize export data collections
        exp_pred_export_data = []
        exp_exp_export_data = []
        one_vs_all_export_data = []
        
        # MODE 1: Experimental vs Predicted
        if enable_predicted_mode:
            st.subheader("✨ Experimental vs Predicted Spectra")
            for psm in final_psms:
                try:
                    psm_key = f"{psm.spectrum_id}_{psm.run}"
                    source_idxml = psm_to_idxml_map.get(psm_key, "Unknown")
                    plot, r2_matching, r2_binned = create_exp_vs_pred_plot_streamlit(psm, aligned_df, selected_model, collision_energy, instrument_type, source_idxml)
                    if plot:
                        exp_pred_plots.append(plot)
                        # Store export data (including plot object for export)
                        exp_pred_export_data.append({
                            'plot': plot,  # Include the plot object
                            'psm': psm,
                            'aligned_df': aligned_df,
                            'model': selected_model,
                            'collision_energy': collision_energy,
                            'instrument_type': instrument_type,
                            'source_idxml': source_idxml,
                            'peptidoform': str(psm.peptidoform),
                            'fragment_mass_tolerance': fragment_mass_tolerance,
                            'r2_metric': st.session_state.get('r2_metric', 'both'),
                            'r2_matching': r2_matching,  # Store actual R² value
                            'r2_binned': r2_binned  # Store actual R² value
                        })
                except Exception as e:
                    logger.error(f"Error creating exp vs pred plot for {psm.spectrum_id}: {e}")
            
            if exp_pred_plots:
                for plot in exp_pred_plots:
                    st.bokeh_chart(plot, use_container_width=True)
                st.success(f"✅ Generated {len(exp_pred_plots)} experimental vs predicted plots")
            
            # Store export data in session state
            st.session_state['exp_pred_export_data'] = exp_pred_export_data
        
        # MODE 2: All-vs-All
        if enable_all_vs_all_mode and len(final_psms) >= 2:
            st.subheader("🔄 All-vs-All Replicate Comparisons")
            
            # Get grouping preference
            group_by_unmodified = st.session_state.get('group_by_unmodified', False)
            
            # Group PSMs based on user preference
            if group_by_unmodified:
                # Group by unmodified sequence (ignoring modifications)
                st.info("🔗 Grouping by unmodified sequence - comparing PSMs with different modifications")
                psms_by_group = {}
                for psm in final_psms:
                    # Extract unmodified sequence from peptidoform
                    try:
                        unmodified_seq = Peptidoform(str(psm.peptidoform)).sequence
                    except:
                        # Fallback: strip out modification characters
                        unmodified_seq = ''.join([c for c in str(psm.peptidoform) if c.isalpha()]).upper()
                    
                    if unmodified_seq not in psms_by_group:
                        psms_by_group[unmodified_seq] = []
                    psms_by_group[unmodified_seq].append(psm)
                
                group_label = "unmodified sequence"
            else:
                # Group by exact peptidoform (same modifications)
                st.info("ℹ️ Grouping by peptidoform - only comparing PSMs with identical modifications")
                psms_by_group = {}
                for psm in final_psms:
                    peptidoform_str = str(psm.peptidoform)
                    if peptidoform_str not in psms_by_group:
                        psms_by_group[peptidoform_str] = []
                    psms_by_group[peptidoform_str].append(psm)
                
                group_label = "peptidoform"
            
            total_comparisons = 0
            groups_with_comparisons = 0
            
            for group_key, psm_group in psms_by_group.items():
                if len(psm_group) >= 2:
                    num_comparisons = len(list(itertools.combinations(psm_group, 2)))
                    groups_with_comparisons += 1
                    
                    # If grouped by unmodified, show all peptidoforms in this group
                    if group_by_unmodified:
                        unique_peptidoforms = set([str(psm.peptidoform) for psm in psm_group])
                        if len(unique_peptidoforms) > 1:
                            st.write(f"**{group_key}** ({len(unique_peptidoforms)} peptidoforms): {num_comparisons} comparison(s)")
                            with st.expander(f"Show peptidoforms for {group_key}"):
                                for pf in sorted(unique_peptidoforms):
                                    count = sum(1 for psm in psm_group if str(psm.peptidoform) == pf)
                                    st.write(f"  • {pf} ({count} PSMs)")
                        else:
                            st.write(f"**{group_key}**: {num_comparisons} comparison(s)")
                    else:
                        st.write(f"**{group_key}**: {num_comparisons} comparison(s)")
                    
                    for psm1, psm2 in itertools.combinations(psm_group, 2):
                        try:
                            psm1_key = f"{psm1.spectrum_id}_{psm1.run}"
                            psm2_key = f"{psm2.spectrum_id}_{psm2.run}"
                            source_idxml1 = psm_to_idxml_map.get(psm1_key, "Unknown")
                            source_idxml2 = psm_to_idxml_map.get(psm2_key, "Unknown")
                            plot, r2_matching, r2_binned, comparison_info = create_exp_vs_exp_plot_streamlit(psm1, psm2, aligned_df, source_idxml1, source_idxml2)
                            if plot:
                                exp_exp_plots.append(plot)
                                total_comparisons += 1
                                # Store export data (including plot object for export)
                                exp_exp_export_data.append({
                                    'plot': plot,  # Include the plot object
                                    'psm1': psm1,
                                    'psm2': psm2,
                                    'aligned_df': aligned_df,
                                    'source_idxml1': source_idxml1,
                                    'source_idxml2': source_idxml2,
                                    'peptidoform': group_key,
                                    'fragment_mass_tolerance': fragment_mass_tolerance,
                                    'r2_metric': st.session_state.get('r2_metric', 'both'),
                                    'r2_matching': r2_matching,  # Store actual R² value
                                    'r2_binned': r2_binned  # Store actual R² value
                                })
                        except Exception as e:
                            logger.error(f"Error creating exp vs exp plot: {e}")
            
            if exp_exp_plots:
                for plot in exp_exp_plots:
                    st.bokeh_chart(plot, use_container_width=True)
                st.success(f"✅ Generated {total_comparisons} replicate comparison plots from {groups_with_comparisons} {group_label} group(s)")
            else:
                st.warning(f"No replicate comparison plots could be generated. Each {group_label} needs at least 2 PSMs for comparison.")
            
            # Store export data in session state
            st.session_state['exp_exp_export_data'] = exp_exp_export_data
        
        # MODE 3: One-vs-All
        if enable_one_vs_all_mode and reference_runs and len(final_psms) >= 2:
            st.subheader("📍 One-vs-All Replicate Comparisons")
            if len(reference_runs) == 1:
                st.info(f"Reference: **{reference_runs[0]}** (lower/mirrored spectrum)")
            else:
                st.info(f"References: **{', '.join(reference_runs)}** (lower/mirrored spectra)")
            
            psms_by_peptidoform = {}
            for psm in final_psms:
                peptidoform_str = str(psm.peptidoform)
                if peptidoform_str not in psms_by_peptidoform:
                    psms_by_peptidoform[peptidoform_str] = []
                psms_by_peptidoform[peptidoform_str].append(psm)
            
            for peptidoform_str, psm_group in psms_by_peptidoform.items():
                # Collect PSMs from all reference runs
                reference_psms = [psm for psm in psm_group if psm.run in reference_runs]
                # Collect PSMs from all non-reference runs
                other_psms = [psm for psm in psm_group if psm.run not in reference_runs]
                
                if reference_psms and other_psms:
                    st.write(f"**{peptidoform_str}**: {len(reference_psms) * len(other_psms)} comparison(s)")
                    for ref_psm in reference_psms:
                        for other_psm in other_psms:
                            try:
                                ref_key = f"{ref_psm.spectrum_id}_{ref_psm.run}"
                                other_key = f"{other_psm.spectrum_id}_{other_psm.run}"
                                source_idxml_ref = psm_to_idxml_map.get(ref_key, "Unknown")
                                source_idxml_other = psm_to_idxml_map.get(other_key, "Unknown")
                                plot, r2_matching, r2_binned, comparison_info = create_exp_vs_exp_plot_streamlit(other_psm, ref_psm, aligned_df, source_idxml_other, source_idxml_ref)
                                if plot:
                                    one_vs_all_plots.append(plot)
                                    # Store export data (including plot object for export)
                                    one_vs_all_export_data.append({
                                        'plot': plot,  # Include the plot object
                                        'other_psm': other_psm,
                                        'ref_psm': ref_psm,
                                        'aligned_df': aligned_df,
                                        'source_idxml_other': source_idxml_other,
                                        'source_idxml_ref': source_idxml_ref,
                                        'peptidoform': peptidoform_str,
                                        'fragment_mass_tolerance': fragment_mass_tolerance,
                                        'r2_metric': st.session_state.get('r2_metric', 'both'),
                                        'r2_matching': r2_matching,  # Store actual R² value
                                        'r2_binned': r2_binned  # Store actual R² value
                                    })
                            except Exception as e:
                                logger.error(f"Error creating one-vs-all plot: {e}")
            
            if one_vs_all_plots:
                for plot in one_vs_all_plots:
                    st.bokeh_chart(plot, use_container_width=True)
                st.success(f"✅ Generated {len(one_vs_all_plots)} one-vs-all comparison plots")
            
            # Store export data in session state
            st.session_state['one_vs_all_export_data'] = one_vs_all_export_data
        
        # Cache analysis results if this was a fresh computation (not using cached data)
        if not use_cached:
            st.session_state['analysis_cache_key'] = analysis_cache_key
            st.session_state['cached_analysis_results'] = {
                'aligned_df': aligned_df,
                'final_psms': final_psms,
                'psm_to_idxml_map': psm_to_idxml_map,
                'exp_pred_plots': exp_pred_plots,
                'exp_exp_plots': exp_exp_plots,
                'one_vs_all_plots': one_vs_all_plots,
                'exp_pred_export_data': exp_pred_export_data,
                'exp_exp_export_data': exp_exp_export_data,
                'one_vs_all_export_data': one_vs_all_export_data
            }
            logger.info("Cached analysis results for future use")
        
        # ==========================================
        # PLOT EXPORT SECTION (All Modes)
        # ==========================================
        # Collect all export data
        all_export_data = []
        export_mode_labels = []
        
        if exp_pred_export_data:
            all_export_data.extend([{'data': d, 'mode': 'exp_vs_pred'} for d in exp_pred_export_data])
            export_mode_labels.append(f"Exp vs Pred ({len(exp_pred_export_data)})")
        
        if exp_exp_export_data:
            all_export_data.extend([{'data': d, 'mode': 'exp_vs_exp'} for d in exp_exp_export_data])
            export_mode_labels.append(f"All-vs-All ({len(exp_exp_export_data)})")
        
        if one_vs_all_export_data:
            all_export_data.extend([{'data': d, 'mode': 'one_vs_all'} for d in one_vs_all_export_data])
            export_mode_labels.append(f"One-vs-All ({len(one_vs_all_export_data)})")
        
        if all_export_data:
            st.markdown("---")
            st.subheader("💾 Export All Plots")
            st.info(f"Export publication-ready plots: summary plots (peptide counts & R² correlations) + individual mirror plots for all comparisons (PNG/SVG)")
            st.caption(f"📊 Available data: {', '.join(export_mode_labels)} = **{len(all_export_data)} comparisons**")
            
            # Use form to prevent rerunning analysis when export button is clicked
            with st.form("export_form"):
                # Configuration
                export_dir = st.text_input(
                    "📁 Export Directory",
                    value=st.session_state.get('plot_output_dir', './output_plots'),
                    help="Directory where plot files will be saved"
                )
                
                # Export button (submit button in form)
                export_clicked = st.form_submit_button("📤 Export All Plots", type="primary", use_container_width=True)
            
            # Handle export outside form to avoid nested logic issues
            if export_clicked:
                # Store in session state for next time
                st.session_state['plot_output_dir'] = export_dir
                
                # Create export directory if it doesn't exist
                export_path = Path(export_dir)
                try:
                    export_path.mkdir(parents=True, exist_ok=True)
                    st.info(f"📁 Export directory: {export_path.absolute()}")
                except Exception as e:
                    st.error(f"❌ Failed to create export directory: {e}")
                    logger.error(f"Directory creation failed: {e}")
                    st.stop()
                
                with st.spinner("Creating summary plots with matplotlib..."):
                    exported_files = []
                    
                    # Prepare R² data for overview plots from all export data
                    r2_data_rows = []
                    
                    def safe_float(value, default=None):
                        """Convert value to float, handling None, strings, and invalid values."""
                        if value is None:
                            return default
                        try:
                            return float(value)
                        except (ValueError, TypeError):
                            return default
                    
                    for item in all_export_data:
                        plot_data = item['data']
                        mode = item['mode']
                        peptidoform_str = plot_data.get('peptidoform', 'Unknown')
                        
                        # Extract R² values - now using stored numeric values
                        r2_matching = safe_float(plot_data.get('r2_matching'))
                        r2_binned = safe_float(plot_data.get('r2_binned'))
                        
                        # Only add rows with at least one valid R² value
                        if r2_matching is not None or r2_binned is not None:
                            r2_data_rows.append({
                                'peptidoform': peptidoform_str,
                                'r2_matching': r2_matching,
                                'r2_binned': r2_binned,
                                'mode': mode
                            })
                    
                    if r2_data_rows:
                        r2_df = pd.DataFrame(r2_data_rows)
                        
                        # Create overview plots for each metric
                        for metric in ['r2_matching', 'r2_binned']:
                            metric_label = 'matching' if metric == 'r2_matching' else 'binned'
                            
                            try:
                                # Create matplotlib figure
                                fig, summary_df = create_peptide_overview_plot(
                                    r2_data=r2_df,
                                    r2_metric=metric,
                                    top_n=None  # Include all peptides
                                )
                                
                                if fig is not None:
                                    # Export as PNG
                                    png_path = export_overview_plot(fig, str(export_path), f"r2_{metric_label}")
                                    if png_path:
                                        exported_files.append(png_path)
                                    
                                    # Also export as SVG
                                    try:
                                        svg_filename = f"peptide_count_r2_{metric_label}.svg"
                                        svg_path = export_path / svg_filename
                                        fig.savefig(svg_path, dpi=300, bbox_inches="tight", format='svg')
                                        exported_files.append(str(svg_path))
                                        logger.info(f"Saved SVG: {svg_path}")
                                    except Exception as e:
                                        logger.warning(f"SVG export failed: {e}")
                                    
                                    # Close figure
                                    import matplotlib.pyplot as plt
                                    plt.close(fig)
                                    
                                    st.success(f"✅ Created {metric_label} R² summary plot")
                                else:
                                    st.warning(f"⚠️ Could not create {metric_label} overview plot (no valid R² data)")
                            
                            except Exception as e:
                                st.error(f"❌ Failed to create {metric_label} plot: {e}")
                                logger.error(f"Overview plot creation failed for {metric}: {e}", exc_info=True)
                        
                        if exported_files:
                            st.success(f"✅ Successfully exported {len(exported_files)} summary plot file(s) to `{export_path.absolute()}`")
                            with st.expander("📂 View exported files"):
                                for filepath in exported_files:
                                    st.code(filepath)
                        else:
                            st.warning("No plots were exported")
                    else:
                        st.warning("No R² data available for summary plots")
                
                # Export individual mirror plots
                with st.spinner("Creating individual mirror plots with matplotlib..."):
                    individual_plot_files = []
                    
                    for idx, item in enumerate(all_export_data):
                        try:
                            plot_data = item['data']
                            mode = item['mode']
                            peptidoform_str = plot_data.get('peptidoform', 'Unknown')
                            r2_matching = plot_data.get('r2_matching')
                            r2_binned = plot_data.get('r2_binned')
                            
                            # Create appropriate mirror plot based on mode
                            if mode == 'exp_vs_pred':
                                # Experimental vs Predicted
                                psm = plot_data.get('psm')
                                aligned_df = plot_data.get('aligned_df')
                                
                                if psm and aligned_df is not None:
                                    # Extract spectrum data properly
                                    psm_aligned = aligned_df[aligned_df.psm_id == psm.spectrum_id] if not aligned_df.empty else pd.DataFrame()
                                    
                                    # Get experimental data using proper extraction
                                    exp_df = st_tools.extract_annotations(psm.spectrum, psm_aligned)
                                    
                                    # Get predicted data - need to predict it
                                    try:
                                        model = plot_data.get('model')
                                        collision_energy = plot_data.get('collision_energy', 35)
                                        instrument_type = plot_data.get('instrument_type', 'TIMSTOF')
                                        pred_df = st_tools.predict_koina_ions(
                                            str(psm.peptidoform), 
                                            psm.peptidoform.precursor_charge,
                                            model,
                                            collision_energy,
                                            instrument_type
                                        )
                                    except Exception as pred_error:
                                        logger.warning(f"Failed to predict spectrum for export: {pred_error}")
                                        pred_df = pd.DataFrame({'mz': [], 'intensity': [], 'annotation': []})
                                    
                                    if not exp_df.empty or not pred_df.empty:
                                        title = f"Exp vs Pred: {peptidoform_str}"
                                        source_idxml = plot_data.get('source_idxml', 'Unknown')
                                        if source_idxml != 'Unknown':
                                            title = f"{title}\n{source_idxml}"
                                        
                                        # Create matplotlib mirror plot
                                        from spectrum_engine import EXP_COLORS, PRED_COLORS
                                        fig = create_matplotlib_mirror_plot(
                                            exp_df, pred_df,
                                            upper_colors=EXP_COLORS,
                                            lower_colors=PRED_COLORS,
                                            title=title,
                                            r2_matching=r2_matching,
                                            r2_binned=r2_binned,
                                            r2_metric=plot_data.get('r2_metric', 'both')
                                        )
                                        
                                        if fig:
                                            # Generate safe filename
                                            safe_peptide = peptidoform_str.replace('/', '_').replace('[', '').replace(']', '')
                                            filename_base = f"exp_vs_pred_{safe_peptide}_{idx}"
                                            
                                            # Export as PNG and SVG
                                            files = export_matplotlib_mirror_plot(fig, str(export_path), filename_base)
                                            individual_plot_files.extend(files)
                            
                            elif mode in ['exp_vs_exp', 'one_vs_all']:
                                # Experimental vs Experimental comparisons
                                if mode == 'exp_vs_exp':
                                    psm1 = plot_data.get('psm1')
                                    psm2 = plot_data.get('psm2')
                                    source_idxml1 = plot_data.get('source_idxml1', 'Unknown')
                                    source_idxml2 = plot_data.get('source_idxml2', 'Unknown')
                                else:  # one_vs_all
                                    psm1 = plot_data.get('other_psm')
                                    psm2 = plot_data.get('ref_psm')
                                    source_idxml1 = plot_data.get('source_idxml_other', 'Unknown')
                                    source_idxml2 = plot_data.get('source_idxml_ref', 'Unknown')
                                
                                aligned_df = plot_data.get('aligned_df')
                                
                                if psm1 and psm2 and aligned_df is not None:
                                    # Extract spectrum data for both PSMs using the proper extraction function
                                    psm1_aligned = aligned_df[aligned_df.psm_id == psm1.spectrum_id] if not aligned_df.empty else pd.DataFrame()
                                    psm2_aligned = aligned_df[aligned_df.psm_id == psm2.spectrum_id] if not aligned_df.empty else pd.DataFrame()
                                    
                                    # Use extract_annotations to properly format the spectrum data
                                    exp1_df = st_tools.extract_annotations(psm1.spectrum, psm1_aligned)
                                    exp2_df = st_tools.extract_annotations(psm2.spectrum, psm2_aligned)
                                    
                                    if not exp1_df.empty or not exp2_df.empty:
                                        mode_label = "One-vs-All" if mode == 'one_vs_all' else "All-vs-All"
                                        title = f"{mode_label}: {peptidoform_str}\n{source_idxml1} vs {source_idxml2}"
                                        
                                        # Create matplotlib mirror plot
                                        from spectrum_engine import EXP_COLORS
                                        fig = create_matplotlib_mirror_plot(
                                            exp1_df, exp2_df,
                                            upper_colors=EXP_COLORS,
                                            lower_colors=EXP_COLORS,
                                            title=title,
                                            r2_matching=r2_matching,
                                            r2_binned=r2_binned,
                                            r2_metric=plot_data.get('r2_metric', 'both')
                                        )
                                        
                                        if fig:
                                            # Generate safe filename
                                            safe_peptide = peptidoform_str.replace('/', '_').replace('[', '').replace(']', '')
                                            filename_base = f"{mode}_{safe_peptide}_{idx}"
                                            
                                            # Export as PNG and SVG
                                            files = export_matplotlib_mirror_plot(fig, str(export_path), filename_base)
                                            individual_plot_files.extend(files)
                        
                        except Exception as e:
                            logger.error(f"Failed to create mirror plot {idx}: {e}", exc_info=True)
                    
                    if individual_plot_files:
                        st.success(f"✅ Successfully exported {len(individual_plot_files)} individual mirror plot file(s) to `{export_path.absolute()}`")
                        st.info(f"📊 Created {len(individual_plot_files) // 2} mirror plots (PNG + SVG for each)")
                    else:
                        st.warning("No individual mirror plots were created")
                
                # Export CSV tables
                with st.spinner("Exporting CSV tables..."):
                    csv_files = []
                    
                    # First pass: collect all unique idXML filenames from plot data
                    all_idxml_files = set()
                    for item in all_export_data:
                        plot_data = item['data']
                        mode = item['mode']
                        if mode == 'exp_vs_pred':
                            source = plot_data.get('source_idxml', 'Unknown')
                            if source != 'Unknown':
                                all_idxml_files.add(source)
                        elif mode == 'exp_vs_exp':
                            source1 = plot_data.get('source_idxml1', 'Unknown')
                            source2 = plot_data.get('source_idxml2', 'Unknown')
                            if source1 != 'Unknown':
                                all_idxml_files.add(source1)
                            if source2 != 'Unknown':
                                all_idxml_files.add(source2)
                        elif mode == 'one_vs_all':
                            source1 = plot_data.get('source_idxml_other', 'Unknown')
                            source2 = plot_data.get('source_idxml_ref', 'Unknown')
                            if source1 != 'Unknown':
                                all_idxml_files.add(source1)
                            if source2 != 'Unknown':
                                all_idxml_files.add(source2)
                    
                    # Create shortened filename mapping
                    short_name_map = shorten_filenames(list(all_idxml_files))
                    
                    # 1. Export R² correlation matrix
                    try:
                        # Build spectrum identifier for each comparison
                        spectrum_identifiers = {}
                        r2_matrix_data = []
                        
                        for idx, item in enumerate(all_export_data):
                            plot_data = item['data']
                            mode = item['mode']
                            r2_matching = plot_data.get('r2_matching')
                            r2_binned = plot_data.get('r2_binned')
                            peptidoform_str = plot_data.get('peptidoform', 'Unknown')
                            
                            # Create spectrum identifiers based on mode (using shortened names)
                            if mode == 'exp_vs_pred':
                                psm = plot_data.get('psm')
                                source_idxml = plot_data.get('source_idxml', 'Unknown')
                                short_name = short_name_map.get(source_idxml, source_idxml)
                                if psm:
                                    spec1_id = f"{short_name}|{peptidoform_str}|{psm.spectrum_id}"
                                    spec2_id = f"Predicted|{peptidoform_str}|AI"
                                    r2_matrix_data.append({
                                        'Spectrum_1': spec1_id,
                                        'Spectrum_2': spec2_id,
                                        'Peptidoform': peptidoform_str,
                                        'R2_Matching': r2_matching if r2_matching is not None else 'N/A',
                                        'R2_Binned': r2_binned if r2_binned is not None else 'N/A',
                                        'Comparison_Type': 'Exp_vs_Pred'
                                    })
                            
                            elif mode == 'exp_vs_exp':
                                psm1 = plot_data.get('psm1')
                                psm2 = plot_data.get('psm2')
                                source_idxml1 = plot_data.get('source_idxml1', 'Unknown')
                                source_idxml2 = plot_data.get('source_idxml2', 'Unknown')
                                short_name1 = short_name_map.get(source_idxml1, source_idxml1)
                                short_name2 = short_name_map.get(source_idxml2, source_idxml2)
                                if psm1 and psm2:
                                    spec1_id = f"{short_name1}|{peptidoform_str}|{psm1.spectrum_id}"
                                    spec2_id = f"{short_name2}|{peptidoform_str}|{psm2.spectrum_id}"
                                    r2_matrix_data.append({
                                        'Spectrum_1': spec1_id,
                                        'Spectrum_2': spec2_id,
                                        'Peptidoform': peptidoform_str,
                                        'R2_Matching': r2_matching if r2_matching is not None else 'N/A',
                                        'R2_Binned': r2_binned if r2_binned is not None else 'N/A',
                                        'Comparison_Type': 'Replicate_AllvsAll'
                                    })
                            
                            elif mode == 'one_vs_all':
                                psm1 = plot_data.get('other_psm')
                                psm2 = plot_data.get('ref_psm')
                                source_idxml1 = plot_data.get('source_idxml_other', 'Unknown')
                                source_idxml2 = plot_data.get('source_idxml_ref', 'Unknown')
                                short_name1 = short_name_map.get(source_idxml1, source_idxml1)
                                short_name2 = short_name_map.get(source_idxml2, source_idxml2)
                                if psm1 and psm2:
                                    spec1_id = f"{short_name1}|{peptidoform_str}|{psm1.spectrum_id}"
                                    spec2_id = f"{short_name2}|{peptidoform_str}|{psm2.spectrum_id}"
                                    r2_matrix_data.append({
                                        'Spectrum_1': spec1_id,
                                        'Spectrum_2': spec2_id,
                                        'Peptidoform': peptidoform_str,
                                        'R2_Matching': r2_matching if r2_matching is not None else 'N/A',
                                        'R2_Binned': r2_binned if r2_binned is not None else 'N/A',
                                        'Comparison_Type': 'OnevsAll'
                                    })
                        
                        if r2_matrix_data:
                            r2_df = pd.DataFrame(r2_matrix_data)
                            r2_csv_path = export_path / "r2_correlation_matrix.csv"
                            r2_df.to_csv(r2_csv_path, index=False)
                            csv_files.append(str(r2_csv_path))
                            logger.info(f"Exported R² correlation matrix: {r2_csv_path}")
                    
                    except Exception as e:
                        logger.error(f"Failed to export R² correlation matrix: {e}", exc_info=True)
                        st.warning(f"⚠️ Could not export R² correlation matrix: {e}")
                    
                    # 2. Export spectrum selection table
                    try:
                        # Get the selected spectra data
                        selected_peptidoforms_data = st.session_state.get('selected_peptidoforms_data', [])
                        
                        if selected_peptidoforms_data:
                            # Build selection table
                            selection_data = []
                            for pep_data in selected_peptidoforms_data:
                                peptidoform = pep_data.get('peptidoform', 'Unknown')
                                display_label = pep_data.get('display_label', peptidoform)
                                num_spectra = pep_data.get('num_spectra', 0)
                                source_files = pep_data.get('source_files', [])
                                
                                selection_data.append({
                                    'Peptidoform': peptidoform,
                                    'Display_Label': display_label,
                                    'Num_Spectra': num_spectra,
                                    'Source_Files': '; '.join(source_files) if source_files else 'N/A'
                                })
                            
                            selection_df = pd.DataFrame(selection_data)
                            selection_csv_path = export_path / "selected_spectra.csv"
                            selection_df.to_csv(selection_csv_path, index=False)
                            csv_files.append(str(selection_csv_path))
                            logger.info(f"Exported spectrum selection table: {selection_csv_path}")
                        
                        # Also export detailed PSM information if available
                        if final_psms:
                            # Collect unique source files from psm_to_idxml_map for shortening
                            psm_source_files = set()
                            for psm in final_psms:
                                psm_key = f"{psm.spectrum_id}_{psm.run}"
                                source_idxml = psm_to_idxml_map.get(psm_key, "Unknown")
                                if source_idxml != "Unknown":
                                    psm_source_files.add(source_idxml)
                            
                            # Create shortened names for PSM sources
                            psm_short_name_map = shorten_filenames(list(psm_source_files))
                            
                            psm_data = []
                            for psm in final_psms:
                                psm_key = f"{psm.spectrum_id}_{psm.run}"
                                source_idxml = psm_to_idxml_map.get(psm_key, "Unknown")
                                short_source = psm_short_name_map.get(source_idxml, source_idxml)
                                
                                psm_data.append({
                                    'Spectrum_ID': psm.spectrum_id,
                                    'Run': psm.run,
                                    'Peptidoform': str(psm.peptidoform),
                                    'Sequence': psm.peptidoform.sequence,
                                    'Charge': psm.peptidoform.precursor_charge,
                                    'Source_File': short_source,
                                    'Source_File_Full': source_idxml,
                                    'Precursor_MZ': getattr(psm.peptidoform, 'precursor_mz', 'N/A')
                                })
                            
                            psm_df = pd.DataFrame(psm_data)
                            psm_csv_path = export_path / "selected_psms_detailed.csv"
                            psm_df.to_csv(psm_csv_path, index=False)
                            csv_files.append(str(psm_csv_path))
                            logger.info(f"Exported detailed PSM information: {psm_csv_path}")
                    
                    except Exception as e:
                        logger.error(f"Failed to export spectrum selection tables: {e}", exc_info=True)
                        st.warning(f"⚠️ Could not export spectrum selection tables: {e}")
                    
                    if csv_files:
                        st.success(f"✅ Successfully exported {len(csv_files)} CSV table(s)")
                        with st.expander("📄 View exported CSV files"):
                            for filepath in csv_files:
                                st.code(filepath)
        
        # Summary
        st.subheader("📋 Analysis Summary")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Selected Peptidoforms", len(selected_peptidoforms))
        with col2:
            st.metric("PSMs with Spectra", len(final_psms))
        with col3:
            st.metric("Peak Alignments", len(aligned_df))
        with col4:
            total_plots = len(exp_pred_plots) + len(exp_exp_plots) + len(one_vs_all_plots)
            st.metric("Total Plots", total_plots)
            
    except Exception as e:
        st.error(f"Error in automated analysis: {e}")
        logger.error(f"Automated analysis error: {e}", exc_info=True)

def create_exp_vs_pred_plot_streamlit(psm, aligned_df, model, collision_energy=35, instrument_type="TIMSTOF", source_idxml="Unknown"):
    """Create experimental vs predicted plot for Streamlit (based on notebook).
    
    Returns:
        tuple: (plot, r2_matching, r2_binned) - plot object and R² values
    """
    try:
        # Extract experimental spectrum data
        psm_aligned = aligned_df[aligned_df.psm_id == psm.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        
        # Always extract experimental data, even if no alignments
        if psm.spectrum is not None and psm.spectrum.size() > 0:
            # Create experimental dataframe from spectrum
            exp_mz = [psm.spectrum[i].getMZ() for i in range(psm.spectrum.size())]
            exp_intensity = [psm.spectrum[i].getIntensity() for i in range(psm.spectrum.size())]
            
            if not psm_aligned.empty:
                # If we have alignments, use the annotated version
                exp_df = st_tools.extract_annotations(psm.spectrum, psm_aligned)
            else:
                # If no alignments, create basic experimental dataframe
                exp_df = pd.DataFrame({
                    'mz': exp_mz,
                    'intensity': exp_intensity,
                    'annotation': ['unassigned'] * len(exp_mz)
                })
                logger.info(f"No alignments found for {psm.spectrum_id}, showing unannotated spectrum")
        else:
            # Empty experimental spectrum
            exp_df = pd.DataFrame({'mz': [], 'intensity': [], 'annotation': []})
            logger.warning(f"No experimental spectrum data for {psm.spectrum_id}")
        
        # Get predicted data
        try:
            pred_df = st_tools.predict_koina_ions(str(psm.peptidoform), psm.peptidoform.precursor_charge, model,
                                                 collision_energy, instrument_type)
        except Exception as pred_error:
            logger.error(f"Failed to predict spectrum for {psm.spectrum_id}: {pred_error}")
            # Create empty predicted dataframe
            pred_df = pd.DataFrame({'mz': [], 'intensity': [], 'annotation': []})
        
        # Create mirror plot with source idXML information
        title = f"Exp vs Pred: {psm.run} ({psm.peptidoform})\nSource: {source_idxml}"
        
        # Add alignment statistics to title
        if not psm_aligned.empty:
            title += f"\n{len(psm_aligned)} aligned peaks"
        else:
            title += f"\nNo aligned peaks (tolerance: {st.session_state.get('fragment_mass_tolerance', 'N/A')} Da)"
        
        # Get R² metric preference from session state
        r2_metric = st.session_state.get('r2_metric', 'both')
        
        plot = st_tools.create_mirror_plot(
            exp_df, pred_df, 
            title=title,
            upper_colors=st_tools.EXP_COLORS,
            lower_colors=st_tools.PRED_COLORS,
            mass_tolerance=st.session_state.get('fragment_mass_tolerance', 0.02),
            r2_metric=r2_metric,
            width=800, height=600
        )
        
        # Extract R² values from plot tags
        r2_matching = None
        r2_binned = None
        if plot and hasattr(plot, 'tags'):
            for tag in plot.tags:
                if isinstance(tag, tuple) and len(tag) == 2:
                    if tag[0] == 'r2_matching':
                        r2_matching = tag[1]
                    elif tag[0] == 'r2_binned':
                        r2_binned = tag[1]
        
        return plot, r2_matching, r2_binned
        
    except Exception as e:
        logger.error(f"Error creating exp vs pred plot: {e}")
        # Return a basic plot showing the error
        try:
            if BOKEH_AVAILABLE:
                error_plot = figure(width=800, height=600, title=f"Error creating plot for {psm.spectrum_id}: {str(e)}")
                error_plot.text([400], [300], text=[f"Failed to create plot: {str(e)}"], text_align="center")
                return error_plot, None, None
            else:
                raise ImportError("Bokeh not available")
        except:
            return None, None, None

def create_exp_vs_exp_plot_streamlit(psm1, psm2, aligned_df, source_idxml1="Unknown", source_idxml2="Unknown"):
    """Create experimental vs experimental plot for Streamlit (based on notebook).
    
    Returns:
        tuple: (plot, r2_matching, r2_binned, comparison_info) - plot object, R² values, and comparison info dict
    """
    try:
        # Extract data for both PSMs
        psm1_aligned = aligned_df[aligned_df.psm_id == psm1.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        exp_df1 = st_tools.extract_annotations(psm1.spectrum, psm1_aligned)
        
        psm2_aligned = aligned_df[aligned_df.psm_id == psm2.spectrum_id] if not aligned_df.empty else pd.DataFrame()
        exp_df2 = st_tools.extract_annotations(psm2.spectrum, psm2_aligned)
        
        # Create mirror plot with source idXML information
        title = f"Exp vs Exp: {psm1.run} ({psm1.peptidoform}) vs {psm2.run} ({psm2.peptidoform})\nSources: {source_idxml1} vs {source_idxml2}"
        
        # Get R² metric preference from session state
        r2_metric = st.session_state.get('r2_metric', 'both')
        
        plot = st_tools.create_mirror_plot(
            exp_df1, exp_df2,
            title=title,
            upper_colors=st_tools.EXP_COLORS,
            lower_colors=st_tools.EXP_COLORS,
            mass_tolerance=st.session_state.get('fragment_mass_tolerance', 0.02),
            r2_metric=r2_metric,
            width=800, height=600
        )
        
        # Extract R² values from plot tags
        r2_matching = None
        r2_binned = None
        if plot and hasattr(plot, 'tags'):
            for tag in plot.tags:
                if isinstance(tag, tuple) and len(tag) == 2:
                    if tag[0] == 'r2_matching':
                        r2_matching = tag[1]
                    elif tag[0] == 'r2_binned':
                        r2_binned = tag[1]
        
        # Create comparison info dictionary for Excel export
        comparison_info = {
            'spectrum_1': psm1.spectrum_id,
            'run_1': psm1.run,
            'peptidoform_1': str(psm1.peptidoform),
            'source_1': source_idxml1,
            'spectrum_2': psm2.spectrum_id,
            'run_2': psm2.run,
            'peptidoform_2': str(psm2.peptidoform),
            'source_2': source_idxml2,
            'r2_matching_ions': r2_matching if r2_matching is not None else 0,
            'r2_binned': r2_binned if r2_binned is not None else 0
        }
        
        return plot, r2_matching, r2_binned, comparison_info
        
    except Exception as e:
        logger.error(f"Error creating exp vs exp plot: {e}")
        return None, None, None, None

if __name__ == "__main__":
    main()