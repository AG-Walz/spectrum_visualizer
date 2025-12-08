"""Core spectrum processing - Backend algorithms for MS/MS analysis.

Provides all fundamental functions: load/parse files, generate theoretical spectra,
align peaks, integrate Koina AI predictions, create Bokeh plots, calculate metrics.
Pure backend logic (no Streamlit UI). Used by spectrum_tools and workspace files.
"""
from __future__ import annotations
import sys; from pathlib import Path as P; sys.path.insert(0, str(P(__file__).parent.parent / "utils"))
from shared_imports import (logging, oms, pd, np, Path, read_file, Peptidoform, PSM, PSMList,
                            Optional, List, Tuple, Dict, Iterable, Any, BOKEH_AVAILABLE)
if BOKEH_AVAILABLE:
    from shared_imports import figure
    from bokeh.models import HoverTool, ColumnDataSource, Legend
    from bokeh.layouts import column
else: figure = HoverTool = ColumnDataSource = Legend = column = None

# Streamlit for caching
import streamlit as st

# Conditional imports
try:
    from koinapy import Koina
    KOINA_AVAILABLE = True
except ImportError:
    KOINA_AVAILABLE = False

try:
    from .correlation import spectrum_correlation, new_spectrum_correlation
except ImportError:
    from correlation import spectrum_correlation, new_spectrum_correlation

# Configure logging
logger = logging.getLogger(__name__)

# Constants for visualization
MODEL = "ms2pip_timsTOF2024"
MAX_FRAG_CHARGE = 2
MZ_TOLERANCE = 0.02
INCLUDED_IONS = ["b", "y", "a", "c", "x", "z", "loss", "precursor", "prefix", "immonium"]
EXP_COLORS = {"b": "#0066CC", "y": "#CC0000", "a": "#FFCC00", "c": "#00CC66", 
              "x": "#C64ADF", "z": "#035F26", "loss": "#FF8C00", 
              "precursor": "#8B008B", "prefix": "#2E8B57", "immonium": "#80031C",
              "unannotated": "#AAAAAA", "other": "#666666"}
PRED_COLORS = {"b": "#66B3FF", "y": "#FF6666", "other": "#0C0909"}


# Helper Functions for loading & filtering data
@st.cache_resource(show_spinner="Opening indexed mzML files...")
def load_mzml_experiments_indexed(file_paths: Iterable[str]) -> dict:
    """Load mzML files using indexed access (memory-efficient).
    
    OPTIMIZED: Uses OnDiscMSExperiment for selective spectrum loading.
    Only metadata is loaded; spectra are loaded on-demand when accessed.
    
    CACHED: Indexed file handles are cached and reused across the session.
    
    Args:
        file_paths (Iterable[str]): List of mzML file paths.
    Returns:
        dict: Dictionary mapping file names to OnDiscMSExperiment objects.
    """
    experiments = {}
    # Convert to tuple for hashability (required for caching)
    file_paths = tuple(file_paths) if not isinstance(file_paths, tuple) else file_paths
    
    for mzml_file in file_paths:
        try:
            indexed_exp = oms.OnDiscMSExperiment()
            indexed_exp.openFile(str(mzml_file))
            experiments[Path(mzml_file).name] = indexed_exp
            logger.info(f"Opened indexed mzML: {Path(mzml_file).name} ({indexed_exp.getNrSpectra()} spectra)")
        except Exception as e:
            logger.error(f"Failed to open indexed mzML {mzml_file}: {e}")
            # Fallback to full loading if indexed fails
            exp = oms.MSExperiment()
            oms.MzMLFile().load(str(mzml_file), exp)
            experiments[Path(mzml_file).name] = exp
            logger.warning(f"Fell back to full loading for {Path(mzml_file).name}")
    
    return experiments


@st.cache_resource(show_spinner="Loading mzML files...")
def load_mzml_experiments(file_paths: Iterable[str]) -> dict:
    """Load mzML files into pyopenms MSExperiment objects (full loading).
    
    NOTE: Consider using load_mzml_experiments_indexed() for better memory efficiency.
    
    CACHED: Files are loaded once and reused across the session.
    
    Args:
        file_paths (Iterable[str]): List of mzML file paths.
    Returns:
        dict: Dictionary mapping file names to MSExperiment objects."""
    experiments = {}
    # Convert to tuple for hashability (required for caching)
    file_paths = tuple(file_paths) if not isinstance(file_paths, tuple) else file_paths
    
    for mzml_file in file_paths:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_file, exp)
        experiments[Path(mzml_file).name] = exp
    return experiments


@st.cache_resource(show_spinner="Building spectrum indexes...")
def create_spectrum_lookup(experiments: dict) -> dict:
    """Create SpectrumLookup maps for loaded experiments.
    
    Works with both MSExperiment and OnDiscMSExperiment objects.
    
    CACHED: Lookups are built once and reused.
    
    Args:
        experiments (dict): Dictionary mapping file names to experiment objects.
    Returns:
        dict: Dictionary mapping experiment names to SpectrumLookup objects."""
    lookup_map = {}
    for exp_name, exp in experiments.items():
        try:
            spectrum_lookup = oms.SpectrumLookup()
            spectrum_lookup.readSpectra(exp, "scan=(?<SCAN>\\d+)")
            lookup_map[Path(exp_name).stem] = spectrum_lookup
        except Exception as e:
            logger.error(f"Failed to create lookup for {exp_name}: {e}")
    return lookup_map


def filter_psms_by_peptide(psm_list: PSMList, target_peptide: str) -> PSMList:
    """Filter PSMs by target peptide sequence.
    Args:
        psm_list (PSMList): List of PSMs to filter by.
        target_peptide (str): Peptide sequence to filter by.
    Returns:
        PSMList: Filtered list of PSMs matching the target peptide.
    """
    logger.info(f"DEBUG: filter_psms_by_peptide called with target_peptide: '{target_peptide}'")
    logger.info(f"DEBUG: Input PSM list has {len(psm_list)} PSMs")
    
    filtered_psms = PSMList(psm_list=[])
    match_count = 0
    
    for i, psm in enumerate(psm_list):
        try:
            unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
            if unmodified_peptide == target_peptide:
                filtered_psms.append(psm)
                match_count += 1
                if match_count <= 3:  # Log first few matches
                    logger.info(f"DEBUG: Match {match_count} - PSM spectrum_id: {psm.spectrum_id}, peptidoform: {psm.peptidoform}, unmodified: {unmodified_peptide}")
        except Exception as e:
            if i < 3:  # Log first few errors
                logger.warning(f"DEBUG: Error processing PSM {i}: {e}")
            continue
    
    logger.info(f"DEBUG: Found {len(filtered_psms)} PSMs matching '{target_peptide}' out of {len(psm_list)} total PSMs")
    
    if len(filtered_psms) == 0:
        # Debug: Show available peptides for troubleshooting
        unique_peptides = set()
        for i, psm in enumerate(psm_list[:10]):  # Sample first 10
            try:
                unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
                unique_peptides.add(unmodified_peptide)
            except:
                continue
        logger.info(f"DEBUG: No matches found. Sample available peptides: {list(unique_peptides)}")
    
    return filtered_psms


def attach_spectrum_to_psm(psm: PSM, experiments: dict, lookups: dict) -> PSM:
    """Attach experimental spectrum to PSM object if available, return the PSM.
    
    Args:
        psm: PSM object to attach spectrum to
        experiments: Dictionary mapping file names to MSExperiment objects
        lookups: Dictionary mapping experiment names to SpectrumLookup objects
    
    Returns:
        PSM: The input PSM object with spectrum attached (if found) or None spectrum
    """
    logger.info(f"DEBUG: attach_spectrum_to_psm called for PSM {psm.spectrum_id} (run: {psm.run})")
    
    if not (hasattr(psm, "spectrum") and psm.spectrum is not None):
        exp_name = f"{psm.run}.mzML"
        logger.info(f"DEBUG: Looking for experiment '{exp_name}' and lookup '{psm.run}'")
        logger.info(f"DEBUG: Available experiments: {list(experiments.keys())}")
        logger.info(f"DEBUG: Available lookups: {list(lookups.keys())}")
        
        if exp_name in experiments and psm.run in lookups:
            try:
                idx = lookups[psm.run].findByNativeID(psm.spectrum_id)
                psm.spectrum = experiments[exp_name][idx]
                logger.info(f"DEBUG: Successfully attached spectrum for PSM {psm.spectrum_id}")
            except Exception as e:
                logger.warning(f"DEBUG: Failed to attach spectrum for PSM {psm.spectrum_id}: {e}")
                psm.spectrum = None
        else:
            logger.warning(f"DEBUG: No matching experiment or lookup found for PSM {psm.spectrum_id}")
            psm.spectrum = None
    else:
        logger.info(f"DEBUG: PSM {psm.spectrum_id} already has spectrum attached")
    
    return psm


# Annotation functions    
def generate_theoretical_spectrum(pf: Peptidoform, max_charge: int = 2, included_ions: List[str] = None) -> oms.MSSpectrum:
    """
    Generate a theoretical spectrum for a given peptide form.
    Args:
        pf (PeptideForm): The peptide form to generate the spectrum for.
        max_charge (int): The maximum charge state to consider.
        included_ions (List[str]): List of ion types to include in generation.
    Returns:
        MSSpectrum: The generated theoretical spectrum.
    """
    if included_ions is None:
        included_ions = INCLUDED_IONS
    
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()

    for ion in ("a","b","y","c","x","z"): 
        params.setValue(f"add_{ion}_ions", "true" if ion in included_ions else "false")

    params.setValue("add_metainfo", "true")
    params.setValue("add_losses", "true" if "loss" in included_ions else "false")
    params.setValue("add_precursor_peaks", "true" if "precursor" in included_ions else "false")
    params.setValue("add_first_prefix_ion", "true" if "prefix" in included_ions else "false")
    params.setValue("add_abundant_immonium_ions", "true" if "immonium" in included_ions else "false")

    tsg.setParameters(params)
    spec = oms.MSSpectrum()
    tsg.getSpectrum(spec, oms.AASequence(str(pf).split('/')[0].replace('[','(').replace(']',')')), 
                    1, max_charge)
    if spec.size() > 1: 
        spec.sortByPosition()
    return spec


def convert_spectrum_to_df(spec: oms.MSSpectrum) -> pd.DataFrame:
    """
    Convert a spectrum object to a DataFrame.
    Args:
        spec (MSSpectrum): The spectrum object to convert.
    Returns:
        pd.DataFrame: The resulting DataFrame.
    """
    if spec is None or spec.size() == 0: 
        return pd.DataFrame(columns=["mz","intensity","annotation","charge"])
    sda, ida = spec.getStringDataArrays(), spec.getIntegerDataArrays()
    names = [n.decode() for n in sda[0]] if sda else [""] * spec.size()
    charges = list(ida[0]) if ida else [None] * spec.size()
    return pd.DataFrame({"mz":[p.getMZ() for p in spec], "intensity":[p.getIntensity() for p in spec], 
                         "annotation":names, "charge":charges})


def align_spectra(exp_spec: oms.MSSpectrum, theo_spec: oms.MSSpectrum, tol: float) -> list:
    """
    Align experimental and theoretical spectra.
    Args:
        exp_spec (MSSpectrum): The experimental spectrum.
        theo_spec (MSSpectrum): The theoretical spectrum.
        tol (float): The tolerance for alignment.
    Returns:
        list: A list of aligned pairs.
    """
    if exp_spec is None or theo_spec is None or exp_spec.size() == 0 or theo_spec.size() == 0: 
        return []
    if exp_spec.size() > 1: 
        exp_spec.sortByPosition()
    if theo_spec.size() > 1: 
        theo_spec.sortByPosition()
    try:
        aligner = oms.SpectrumAlignment()
        prm = aligner.getParameters()
        prm.setValue("tolerance", tol)
        prm.setValue("is_relative_tolerance","false")
        aligner.setParameters(prm)
        pairs = []
        aligner.getSpectrumAlignment(pairs, theo_spec, exp_spec)
        return pairs
    except Exception as e:
        logger.warning(f"Spectrum alignment failed: {e}")
        return []

def annotate_peaks(exp_spec: oms.MSSpectrum, theo_df: pd.DataFrame, pairs: list) -> pd.DataFrame:
    """
    Annotate experimental spectra with theoretical information.
    Args:
        exp_spec (MSSpectrum): The experimental spectrum.
        theo_df (pd.DataFrame): The theoretical spectrum DataFrame.
        pairs (list): The list of aligned pairs.
    Returns:
        pd.DataFrame: The annotated DataFrame.
    """
    if exp_spec is None or not pairs: 
        return pd.DataFrame(columns=["exp_idx","theo_idx","mz","intensity",
                                   "annotation","exp_mz","exp_intensity"])
    
    results = []
    for theo_idx, exp_idx in pairs:
        # Check bounds
        if exp_idx >= exp_spec.size() or theo_idx >= len(theo_df):
            continue
            
        # Get experimental peak info
        exp_peak = exp_spec[exp_idx]
        exp_mz = exp_peak.getMZ()
        exp_intensity = exp_peak.getIntensity()
        
        # Get theoretical peak info - use iloc to access by position in the (potentially filtered) dataframe
        theo_row = theo_df.iloc[theo_idx]
        annotation = theo_row["annotation"]
        theo_mz = theo_row["mz"]
        theo_intensity = theo_row["intensity"]
            
        results.append({
            "exp_idx": exp_idx,
            "theo_idx": theo_idx,
            "mz": theo_mz,
            "intensity": theo_intensity,
            "annotation": annotation,
            "exp_mz": exp_mz,
            "exp_intensity": exp_intensity
        })
    
    return pd.DataFrame(results)



def assign_ion_color(ion_name, color_map: dict) -> str:
    """Map ion name to color based on ion type"""
    if ion_name == "unassigned":
        return color_map.get("unannotated", "#666666")
    
    # Handle specific patterns for new ion types
    if "loss" in ion_name.lower() or "-" in ion_name:
        return color_map.get("loss", color_map.get("other", "#AAAAAA"))
    if "precursor" in ion_name.lower() or ion_name.startswith("M"):
        return color_map.get("precursor", color_map.get("other", "#AAAAAA"))
    if "immonium" in ion_name.lower() or ion_name.startswith("I"):
        return color_map.get("immonium", color_map.get("other", "#AAAAAA"))
    if "prefix" in ion_name.lower():
        return color_map.get("prefix", color_map.get("other", "#AAAAAA"))
    
    # Handle standard ions
    for ion_type in ["b", "y", "a", "c", "x", "z"]:
        if ion_name.startswith(ion_type):
            return color_map.get(ion_type, color_map.get("other", "#AAAAAA"))
    
    return color_map.get("other", "#AAAAAA")

def get_ion_family(annotation: str) -> str:
    """Determine ion family from annotation string."""
    if annotation == "unassigned" or not annotation:
        return "unassigned"
    
    annotation_lower = annotation.lower()
    
    # Check for specific patterns for new ion types
    if "loss" in annotation_lower or "-" in annotation:
        return "loss"
    if "precursor" in annotation_lower or annotation.startswith("M"):
        return "precursor"
    if "immonium" in annotation_lower or annotation.startswith("I"):
        return "immonium"
    if "prefix" in annotation_lower:
        return "prefix"
    
    # Check standard ions (b, y, a, c, x, z)
    for ion_type in ["b", "y", "a", "c", "x", "z"]:
        if annotation.startswith(ion_type):
            return ion_type
    
    return "other"

def render_spectrum_peaks(mirror: bool, peak_df: pd.DataFrame, plot: figure, color_map: dict, legend_name: str) -> List[Tuple[str, List]]:
    """
    Render ions on the plot and return legend items.
    Args:
        mirror (bool): Whether the plot is mirrored.
        peak_df (pd.DataFrame): DataFrame containing peak information.
        plot (figure): Bokeh figure object to render the ions on.
        color_map (dict): Mapping of ion types to colors.
        legend_name (str): Name to use in the legend for these ions.
    Returns:
        List[Tuple[str, List]]: Legend items for the rendered ions.
    """
    if peak_df.empty:
        return []

    # Normalize intensities and flip sign if mirror = True
    peak_df = peak_df.copy()
    max_intensity = peak_df["intensity"].max()
    if max_intensity > 0:
        peak_df["intensity"] = peak_df["intensity"] / max_intensity * (-1 if mirror else 1)
    
    # Add metadata for hover and group by ion family using the new function
    peak_df["spectrum_type"] = legend_name
    peak_df["ion_family"] = peak_df["annotation"].apply(get_ion_family)

    # Render each ion family
    legend_items = []
    for ion_type in peak_df["ion_family"].unique():
        sub = peak_df[peak_df["ion_family"] == ion_type]
        if not sub.empty:
            renderer = plot.segment(x0="mz", y0=0, x1="mz", y1="intensity", 
                                  line_color=assign_ion_color(ion_type, color_map), 
                                  line_width=2, alpha=0.8, source=ColumnDataSource(sub))
            legend_items.append((f"{legend_name} {ion_type}", [renderer]))
    
    return legend_items


def setup_bokeh_plot(title: str, width: int = 800, height: int = 500) -> figure:
    """
    Create standardized Bokeh plot with common configuration.
    Args:
        title (str): Title of the plot.
        width (int): Width of the plot.
        height (int): Height of the plot.
    Returns:
        figure: Bokeh figure object with standardized configuration.
    """
    plot = figure(width=width, height=height, title=title,
               tools="pan,wheel_zoom,box_zoom,reset,save,hover",
               x_axis_label="m/z", y_axis_label="Relative Intensity")

    plot.select_one(HoverTool).tooltips = [
        ("m/z", "@mz{0.000}"), ("Intensity", "@intensity{0.000}"),
        ("Annotation", "@annotation"), ("Spectrum", "@spectrum_type")
    ]
    return plot


def create_mirror_plot(upper_spectrum_data, lower_spectrum_data, upper_colors=None, lower_colors=None, 
                      title="Mirror Plot", upper_label="", lower_label="", mass_tolerance=0.02, r2_metric="both", **kwargs):
    """
    Create a mirror plot with two spectra - one on top (positive) and one on bottom (negative).
    
    Args:
        upper_spectrum_data / lower_spectrum_data (pd.DataFrame or dict): 
            Spectrum data with columns/keys ['mz', 'intensity', 'annotation']
        upper_colors / lower_colors (dict, optional): 
            Color mapping for upper/lower spectrum ion types. Uses EXP_COLORS/PRED_COLORS if None.
        title (str): Plot title
        upper_label / lower_label (str): 
            Label for upper/lower spectrum ions in legend
        mass_tolerance (float): Mass tolerance in Da (unused, kept for compatibility)
        r2_metric (str): Which R² metric to display - 'both', 'matching', or 'binned'
        width (int): Plot width
        height (int): Plot height
    
    Returns:
        figure: Bokeh figure object with mirror plot
    """

    # Calculate R² values between the spectra if both are non-empty
    enhanced_title = title
    r2_matching = None
    r2_binned = None
    
    if not upper_spectrum_data.empty and not lower_spectrum_data.empty:
        try:
            # Method 1: Spectrum correlation based on matching ions
            if r2_metric in ["both", "matching"]:
                # Extract annotated ions for correlation calculation
                upper_ions = []
                lower_ions = []
                
                # Get matching annotated peaks
                for _, peak in upper_spectrum_data.iterrows():
                    if 'annotation' in peak and peak['annotation'] != 'unassigned':
                        upper_ions.append((peak['mz'], peak['intensity'], peak['annotation']))
                
                for _, peak in lower_spectrum_data.iterrows():
                    if 'annotation' in peak and peak['annotation'] != 'unassigned':
                        lower_ions.append((peak['mz'], peak['intensity'], peak['annotation']))
                
                # Calculate R² using matching ions method
                r2_matching, _, _, _ = spectrum_correlation(tuple(upper_ions), tuple(lower_ions))
            
            # Method 2: Binned spectrum correlation
            if r2_metric in ["both", "binned"]:
                r2_binned, _, _, _ = new_spectrum_correlation(upper_spectrum_data, lower_spectrum_data)
            
            # Add R² values to title based on user preference
            if r2_metric == "both":
                if r2_matching is not None and r2_matching > 0 and r2_binned is not None and r2_binned > 0:
                    enhanced_title = f"{title}\nR² (matching ions): {r2_matching:.3f} | R² (binned): {r2_binned:.3f}"
                elif r2_binned is not None and r2_binned > 0:
                    enhanced_title = f"{title}\nR² (binned): {r2_binned:.3f}"
                elif r2_matching is not None and r2_matching > 0:
                    enhanced_title = f"{title}\nR² (matching ions): {r2_matching:.3f}"
            elif r2_metric == "matching":
                if r2_matching is not None and r2_matching > 0:
                    enhanced_title = f"{title}\nR² (matching ions): {r2_matching:.3f}"
            elif r2_metric == "binned":
                if r2_binned is not None and r2_binned > 0:
                    enhanced_title = f"{title}\nR² (binned): {r2_binned:.3f}"
                
        except Exception as e:
            # If R² calculation fails, keep original title
            enhanced_title = title

    # Create plot with enhanced title
    plot = setup_bokeh_plot(enhanced_title, **kwargs)
    
    # Render both spectra and collect legend items
    all_legend_items = []
    if not upper_spectrum_data.empty:
        all_legend_items.extend(render_spectrum_peaks(False, upper_spectrum_data, plot, upper_colors, upper_label))
    if not lower_spectrum_data.empty:
        all_legend_items.extend(render_spectrum_peaks(True, lower_spectrum_data, plot, lower_colors, lower_label))

    # Add deduplicated legend
    if all_legend_items:
        unique_items = [(label, renderers) for label, renderers in all_legend_items 
                       if label not in [item[0] for item in all_legend_items[:all_legend_items.index((label, renderers))]]]
        plot.add_layout(Legend(items=unique_items, orientation="horizontal", 
                           location="top_center", label_text_font_size="9pt"))
    
    # Set axis range and center line
    plot.y_range.start = -1.2 if not lower_spectrum_data.empty else 0
    plot.y_range.end = 1.2

    all_mz = []
    if not (upper_spectrum_data.empty and lower_spectrum_data.empty):
        all_mz = list(upper_spectrum_data["mz"]) + list(lower_spectrum_data["mz"])
        if all_mz:
            plot.line([min(all_mz), max(all_mz)], [0, 0], line_color="gray", line_width=1, alpha=0.5)

    # Store R² values as metadata on the plot for potential retrieval
    if r2_matching is not None:
        plot.tags.append(("r2_matching", r2_matching))
    if r2_binned is not None:
        plot.tags.append(("r2_binned", r2_binned))

    return plot


# Plotting helper functions
def extract_annotations(exp_spec: oms.MSSpectrum, aligned_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare the DataFrame for plotting by combining experimental and aligned data.
    Args:
        exp_spec: Experimental spectrum data.
        aligned_df: Aligned DataFrame.
    Returns:
        pd.DataFrame: DataFrame ready for plotting.
    """
    # Build list of all experimental peaks
    if exp_spec is None:
        return pd.DataFrame(columns=["mz", "intensity", "annotation"])
    
    # Use existing function instead of manual extraction
    exp_df = convert_spectrum_to_df(exp_spec)[["mz", "intensity"]]
    
    if aligned_df.empty:
        exp_df["annotation"] = "unassigned"
        return exp_df
    
    # Efficient merge instead of loops
    aligned_peaks = aligned_df[["exp_mz", "annotation"]].rename(columns={"exp_mz": "mz"})
    result = exp_df.merge(aligned_peaks, on="mz", how="left")
    result["annotation"] = result["annotation"].fillna("unassigned")
    
    return result[["mz", "intensity", "annotation"]]


def extract_psm_spectrum(psm, aligned_df) -> pd.DataFrame:
    """Extract spectrum data for a single PSM with annotations."""
    psm_aligned = aligned_df[aligned_df.psm_id == psm.spectrum_id] if not aligned_df.empty else pd.DataFrame()
    return extract_annotations(psm.spectrum, psm_aligned)


def create_exp_vs_exp_plot(psm1, psm2, aligned_df):
    """
    Create a mirror plot comparing two experimental spectra.
    Args:
        psm1, psm2: PSM objects with spectrum data
        aligned_df: DataFrame with aligned data for both PSMs
    Returns:
        Bokeh figure object
    """
    # Get data for both PSMs
    exp_df1 = extract_psm_spectrum(psm1, aligned_df)
    exp_df2 = extract_psm_spectrum(psm2, aligned_df)

    title = f"Exp vs Exp: \n{psm1.run} ({psm1.peptidoform}) vs {psm2.run} ({psm2.peptidoform})"
    return create_mirror_plot(exp_df1, exp_df2, title=title, upper_colors=EXP_COLORS, lower_colors=EXP_COLORS)


def convert_peptide_to_unimod(peptide_seq: str) -> str:
    """
    Convert peptide sequence to pyopenMS.AASequence object to convert post translational modifications to UniMod format for Koina.
    Args:
        peptide_seq (str): The peptide sequence to convert.
    Returns:
        str: The UniMod formatted peptide sequence.
    """
    aa_seq = oms.AASequence(str(peptide_seq).split('/')[0].replace('[', '(').replace(']', ')'))
    return str(aa_seq.toUniModString()).replace('(', '[').replace(')', ']')


def predict_koina_ions(peptide_seq: str, charge: int, model: str = "ms2pip_timsTOF2024", 
                      collision_energy: int = 35, instrument_type: str = "TIMSTOF") -> pd.DataFrame:
    """Predict fragment ion intensities using Koina.
    Args:
        peptide_seq (str): The peptide sequence in UniMod format.
        charge (int): The precursor charge state.
        model (str): The model to use for prediction.
        collision_energy (int): Collision energy for fragmentation.
        instrument_type (str): Mass spectrometer type.
    Returns:
        pd.DataFrame: DataFrame with columns ['mz', 'intensity', 'annotation']
    """
    unimod_seq = convert_peptide_to_unimod(peptide_seq)
    koina = Koina(model)
    
    # Create input dataframe with all possible parameters
    inputs = pd.DataFrame({
        "peptide_sequences": [unimod_seq],
        "precursor_charges": [charge],
        "collision_energies": [collision_energy],
        "instrument_types": [instrument_type]
    })
    
    # koinapy only takes the input it requires for the current model
    predictions = koina.predict(inputs)
    
    # Remove byte string artifacts from annotations
    predictions["annotation"] = predictions["annotation"].astype(str).str.replace("^b'|'$","", regex=True)
    
    return pd.DataFrame({'mz': predictions['mz'].astype(float),
                         'intensity': predictions['intensities'].astype(float),
                         'annotation': predictions['annotation']
    })


def create_exp_vs_pred_plot(psm, aligned_df, model=MODEL, collision_energy=35, instrument_type="TIMSTOF"):
    """
    Create a mirror plot comparing experimental vs predicted spectrum for a PSM.
    Args:
        psm: PSM object with spectrum data
        aligned_df: DataFrame with aligned experimental data
        model: Koina model to use for predictions
        collision_energy: Collision energy for fragmentation
        instrument_type: Mass spectrometer type
    Returns:
        Bokeh figure object
    """
    # Get experimental data as dataframe
    exp_df = extract_psm_spectrum(psm, aligned_df)

    # Get predicted data using predict_koina_ions
    pred_df = predict_koina_ions(str(psm.peptidoform), psm.peptidoform.precursor_charge, model, 
                                collision_energy, instrument_type)
    
    title=f"Exp vs Pred: \n {psm.run}({psm.peptidoform})"

    return create_mirror_plot(exp_df, pred_df, title=title, upper_colors=EXP_COLORS, lower_colors=PRED_COLORS)


def create_mirror_plot_for_export(upper_spectrum_data, lower_spectrum_data, upper_colors=None, lower_colors=None,
                                   title="Mirror Plot", mass_tolerance=0.02, r2_metric="both", 
                                   width=1200, height=800):
    """
    Create a publication-ready mirror plot without legend, no background, with peak annotations.
    
    Args:
        upper_spectrum_data / lower_spectrum_data: Spectrum data DataFrames
        upper_colors / lower_colors: Color mappings for ion types
        title: Plot title
        mass_tolerance: Mass tolerance in Da
        r2_metric: Which R² to display
        width, height: Plot dimensions
        
    Returns:
        Bokeh figure object optimized for export
    """
    from bokeh.models import Label
    
    # Set default colors
    if upper_colors is None:
        upper_colors = EXP_COLORS
    if lower_colors is None:
        lower_colors = PRED_COLORS
    
    # Calculate R² values
    enhanced_title = title
    r2_matching = None
    r2_binned = None
    
    if not upper_spectrum_data.empty and not lower_spectrum_data.empty:
        try:
            if r2_metric in ["both", "matching"]:
                upper_ions = [(peak['mz'], peak['intensity'], peak['annotation']) 
                             for _, peak in upper_spectrum_data.iterrows() 
                             if 'annotation' in peak and peak['annotation'] != 'unassigned']
                lower_ions = [(peak['mz'], peak['intensity'], peak['annotation']) 
                             for _, peak in lower_spectrum_data.iterrows() 
                             if 'annotation' in peak and peak['annotation'] != 'unassigned']
                r2_matching, _, _, _ = spectrum_correlation(tuple(upper_ions), tuple(lower_ions))
            
            if r2_metric in ["both", "binned"]:
                r2_binned, _, _, _ = new_spectrum_correlation(upper_spectrum_data, lower_spectrum_data)
            
            if r2_metric == "both" and r2_matching is not None and r2_binned is not None:
                enhanced_title = f"{title}\nR² (matching): {r2_matching:.3f} | R² (binned): {r2_binned:.3f}"
            elif r2_metric == "matching" and r2_matching is not None:
                enhanced_title = f"{title}\nR² (matching): {r2_matching:.3f}"
            elif r2_metric == "binned" and r2_binned is not None:
                enhanced_title = f"{title}\nR² (binned): {r2_binned:.3f}"
        except Exception as e:
            logger.error(f"Error calculating R²: {e}")
    
    # Create plot with transparent background
    plot = figure(title=enhanced_title, width=width, height=height,
                 x_axis_label="m/z", y_axis_label="Relative Intensity",
                 toolbar_location="above", tools="pan,wheel_zoom,box_zoom,reset,save",
                 background_fill_color=None, border_fill_color=None,
                 outline_line_color=None)
    
    # Style plot for publication
    plot.title.text_font_size = "14pt"
    plot.title.align = "center"
    plot.xaxis.axis_label_text_font_size = "12pt"
    plot.yaxis.axis_label_text_font_size = "12pt"
    plot.xaxis.major_label_text_font_size = "10pt"
    plot.yaxis.major_label_text_font_size = "10pt"
    
    # Render upper spectrum (positive intensities) with annotations
    if not upper_spectrum_data.empty:
        upper_df = upper_spectrum_data.copy()
        max_int = upper_df["intensity"].max()
        if max_int > 0:
            upper_df["intensity"] = upper_df["intensity"] / max_int
        
        # Group by ion family and plot
        for _, peak in upper_df.iterrows():
            ion_family = get_ion_family(peak['annotation'])
            color = assign_ion_color(peak['annotation'], upper_colors)
            
            # Draw peak stem
            plot.line([peak['mz'], peak['mz']], [0, peak['intensity']], 
                     line_color=color, line_width=2, alpha=0.8)
            
            # Add annotation label at top of peak (only for annotated peaks)
            if peak['annotation'] != 'unassigned' and peak['annotation']:
                label = Label(x=peak['mz'], y=peak['intensity'], 
                            text=peak['annotation'],
                            text_font_size='8pt', text_color=color,
                            text_align='center', text_baseline='bottom',
                            x_offset=-5, y_offset=2)
                plot.add_layout(label)
    
    # Render lower spectrum (negative intensities) with annotations
    if not lower_spectrum_data.empty:
        lower_df = lower_spectrum_data.copy()
        max_int = lower_df["intensity"].max()
        if max_int > 0:
            lower_df["intensity"] = lower_df["intensity"] / max_int * -1  # Mirror
        
        # Group by ion family and plot
        for _, peak in lower_df.iterrows():
            ion_family = get_ion_family(peak['annotation'])
            color = assign_ion_color(peak['annotation'], lower_colors)
            
            # Draw peak stem
            plot.line([peak['mz'], peak['mz']], [0, peak['intensity']], 
                     line_color=color, line_width=2, alpha=0.8)
            
            # Add annotation label at bottom of peak (only for annotated peaks)
            if peak['annotation'] != 'unassigned' and peak['annotation']:
                label = Label(x=peak['mz'], y=peak['intensity'], 
                            text=peak['annotation'],
                            text_font_size='8pt', text_color=color,
                            text_align='center', text_baseline='top',
                            x_offset=-5, y_offset=-2)
                plot.add_layout(label)
    
    # Set axis ranges
    plot.y_range.start = -1.2
    plot.y_range.end = 1.2
    
    # Add center line
    if not (upper_spectrum_data.empty and lower_spectrum_data.empty):
        all_mz = list(upper_spectrum_data["mz"]) + list(lower_spectrum_data["mz"])
        if all_mz:
            plot.line([min(all_mz), max(all_mz)], [0, 0], 
                     line_color="black", line_width=1, alpha=0.5)
    
    # Store R² values as tags
    if r2_matching is not None:
        plot.tags.append(("r2_matching", r2_matching))
    if r2_binned is not None:
        plot.tags.append(("r2_binned", r2_binned))
    
    return plot


def export_plot(plot, filename_base, output_dir, formats=['png', 'svg'], driver=None):
    """
    Export a Bokeh plot to multiple formats.
    
    Args:
        plot: Bokeh figure object
        filename_base: Base filename (without extension)
        output_dir: Directory to save files
        formats: List of formats to export ('png', 'svg', 'html')
        driver: Optional pre-existing webdriver instance (for batch exports)
        
    Returns:
        List of created file paths
    """
    from pathlib import Path
    import os
    
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    created_files = []
    driver_created = False
    
    # Try to import export functions
    try:
        from bokeh.io import export_png, export_svgs
        from selenium import webdriver
        from selenium.webdriver.firefox.options import Options as FirefoxOptions
        export_available = True
    except ImportError as e:
        logger.warning(f"Bokeh export dependencies not available: {e}. Install selenium, pillow, and geckodriver.")
        export_available = False
    
    if not export_available:
        # Fallback: save as HTML
        from bokeh.plotting import save, output_file
        filepath = output_path / f"{filename_base}.html"
        output_file(str(filepath))
        save(plot)
        created_files.append(str(filepath))
        logger.info(f"Exported HTML (fallback): {filepath}")
        return created_files
    
    # Create driver if not provided
    if driver is None and ('png' in formats or 'svg' in formats):
        try:
            firefox_options = FirefoxOptions()
            firefox_options.add_argument('--headless')
            firefox_options.add_argument('--disable-gpu')
            driver = webdriver.Firefox(options=firefox_options)
            driver_created = True
            logger.info("Created Firefox webdriver for export")
        except Exception as driver_error:
            logger.error(f"Failed to create Firefox webdriver: {driver_error}")
            # Fall back to HTML
            from bokeh.plotting import save, output_file
            filepath = output_path / f"{filename_base}.html"
            output_file(str(filepath))
            save(plot)
            created_files.append(str(filepath))
            logger.info(f"Exported HTML (fallback due to driver error): {filepath}")
            return created_files
    
    for fmt in formats:
        filepath = output_path / f"{filename_base}.{fmt}"
        try:
            if fmt == 'png':
                export_png(plot, filename=str(filepath), webdriver=driver)
                created_files.append(str(filepath))
                logger.info(f"Exported PNG: {filepath.name}")
            elif fmt == 'svg':
                plot.output_backend = 'svg'
                export_svgs(plot, filename=str(filepath), webdriver=driver)
                created_files.append(str(filepath))
                logger.info(f"Exported SVG: {filepath.name}")
            elif fmt == 'html':
                from bokeh.plotting import save, output_file
                output_file(str(filepath))
                save(plot)
                created_files.append(str(filepath))
                logger.info(f"Exported HTML: {filepath.name}")
        except Exception as e:
            logger.error(f"Failed to export {fmt} for {filename_base}: {e}")
            # Try HTML fallback for failed PNG/SVG
            if fmt in ['png', 'svg']:
                try:
                    from bokeh.plotting import save, output_file
                    html_path = output_path / f"{filename_base}.html"
                    output_file(str(html_path))
                    save(plot)
                    created_files.append(str(html_path))
                    logger.info(f"Exported HTML (fallback for {fmt}): {html_path.name}")
                except Exception as html_error:
                    logger.error(f"HTML fallback also failed: {html_error}")
    
    # Cleanup: close driver if we created it
    if driver_created and driver is not None:
        try:
            driver.quit()
            logger.info("Closed Firefox webdriver")
        except Exception as e:
            logger.warning(f"Error closing webdriver: {e}")
    
    return created_files


def create_peptide_overview_plot(r2_data: pd.DataFrame, r2_metric: str = 'r2_matching', 
                                 top_n: Optional[int] = None) -> Tuple:
    """
    Create an overview plot showing peptide counts and R² distribution.
    
    Args:
        r2_data: DataFrame with columns ['peptidoform', 'r_squared'] (or r2_matching/r2_binned)
        r2_metric: Which R² metric to use ('r2_matching' or 'r2_binned')
        top_n: Optional limit on number of peptides to display (None = all)
    
    Returns:
        Tuple of (matplotlib.figure.Figure, summary DataFrame)
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend for server
    except ImportError:
        logger.error("matplotlib not available for overview plot")
        return None, None
    
    # Determine R² column name
    r2_column = r2_metric if r2_metric in r2_data.columns else 'r_squared'
    
    if r2_column not in r2_data.columns:
        logger.error(f"R² column '{r2_column}' not found in data. Available: {r2_data.columns.tolist()}")
        return None, None
    
    # Count peptidoforms and collect R² values
    counts = r2_data.groupby("peptidoform").size().rename("peptide_count")
    r2_lists = (
        r2_data.groupby("peptidoform")[r2_column]
              .apply(lambda x: [v for v in x if pd.notna(v) and v != 'N/A'])
              .rename("r2_values")
    )
    
    summary = pd.concat([counts, r2_lists], axis=1)
    
    # Filter out peptides with no valid R² values
    summary = summary[summary['r2_values'].apply(len) > 0].copy()
    
    if summary.empty:
        logger.warning("No valid R² data for overview plot")
        return None, None
    
    # Sort by count (descending), then by peptide name
    summary = summary.sort_values(
        by=["peptide_count", "peptidoform"],
        ascending=[False, True]
    ).reset_index()
    
    # Limit to top N if specified
    if top_n is not None and top_n > 0:
        summary = summary.head(top_n)
    
    # Prepare data for plotting
    peptides = summary["peptidoform"].tolist()
    counts_vals = summary["peptide_count"].to_numpy()
    y_pos = np.arange(len(peptides))
    r2_data_for_boxplot = summary["r2_values"].tolist()
    
    # Create figure
    fig, ax_bar = plt.subplots(figsize=(10, max(6, len(peptides) * 0.3)))
    
    # A) Horizontal bar plot of peptide counts
    bar_container = ax_bar.barh(
        y_pos,
        counts_vals,
        align="center",
        color="#f28e61",
        edgecolor="none",
        label="Peptide Count",
    )
    
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(peptides, fontsize=9)
    ax_bar.set_xlabel("Peptide count")
    ax_bar.invert_yaxis()  # highest count on top
    ax_bar.set_xlim(left=0, right=max(counts_vals) * 1.25)
    
    # B) Twin axis for R² boxplots
    ax_r2 = ax_bar.twiny()
    
    # Boxplot
    box = ax_r2.boxplot(
        r2_data_for_boxplot,
        vert=False,
        positions=y_pos,
        widths=0.6,
        manage_ticks=False,
        patch_artist=True,
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=1),
        medianprops=dict(color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        flierprops=dict(
            marker="o",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=3,
            linestyle="none",
        )
    )
    
    # Overlay individual R² points (jittered vertically)
    all_x = []
    all_y = []
    for idx, r2_list in enumerate(r2_data_for_boxplot):
        all_x.extend(r2_list)
        jitter = (np.random.rand(len(r2_list)) - 0.5) * 0.15
        all_y.extend(idx + jitter)
    
    ax_r2.plot(
        all_x,
        all_y,
        "o",
        color="black",
        markersize=3,
        label=r"Correlation [$R^2$]",
    )
    
    # R² axis formatting
    ax_r2.set_xlim(0, 1.0)
    metric_label = "Matching Ions" if r2_metric == "r2_matching" else "Binned"
    ax_r2.set_xlabel(f"Synthetic peptide correlation [$R^2$ - {metric_label}]")
    ax_r2.set_ylim(ax_bar.get_ylim())
    ax_r2.xaxis.grid(True, linestyle="--", alpha=0.4)
    ax_bar.xaxis.grid(False)
    ax_bar.set_axisbelow(False)
    
    # Legend
    handles = [
        bar_container.patches[0],
        plt.Line2D([0], [0], marker='o', color='black', linestyle='none', markersize=5),
    ]
    labels = ["Peptide Count", r"Correlation [$R^2$]"]
    ax_bar.legend(handles, labels, loc="lower right", frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    return fig, summary


def export_overview_plot(fig, output_dir: str, metric_name: str = "overview") -> Optional[str]:
    """
    Export overview plot to PNG file.
    
    Args:
        fig: Matplotlib figure object
        output_dir: Directory to save the file
        metric_name: Name to include in filename (e.g., 'matching', 'binned')
    
    Returns:
        Path to saved file or None if failed
    """
    if fig is None:
        return None
    
    try:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        filename = f"peptide_count_r2_{metric_name}.png"
        filepath = output_path / filename
        
        fig.savefig(filepath, dpi=300, bbox_inches="tight")
        logger.info(f"Saved overview plot: {filepath}")
        
        # Close figure to free memory
        import matplotlib.pyplot as plt
        plt.close(fig)
        
        return str(filepath)
    
    except Exception as e:
        logger.error(f"Failed to export overview plot: {e}")
        return None