"""
SpecViz module for spectrum visualization and analysis.
Contains functions for generating theoretical spectra, aligning with experimental data,
and creating visualization plots using PyOpenMS and Bokeh.
"""

import pyopenms as oms
import pandas as pd
import itertools
from pathlib import Path
from psm_utils.io import read_file
from psm_utils import Peptidoform, PSM, PSMList
from typing import Iterable, List, Tuple, Dict, Optional
from koinapy import Koina
import logging

from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, Legend
from bokeh.layouts import column

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants for visualization
MODEL = "ms2pip_timsTOF2024"
MAX_FRAG_CHARGE = 2
MZ_TOLERANCE = 0.5
INCLUDED_IONS = ["b", "y", "a", "c", "x", "z"]
EXP_COLORS = {"b": "#0066CC", "y": "#CC0000", "a": "#FFCC00", "c": "#00CC66", 
              "x": "#C64ADF", "z": "#035F26", "unannotated": "#AAAAAA", "other": "#666666"}
PRED_COLORS = {"b": "#66B3FF", "y": "#FF6666", "other": "#0C0909"}


# Helper Functions for loading & filtering data
def load_mzml_experiments(file_paths: Iterable[str]) -> dict:
    """Load mzML files into pyopenms MSExperiment objects.
    Args:
        file_paths (Iterable[str]): List of mzML file paths.
    Returns:
        dict: Dictionary mapping file names to MSExperiment objects."""
    experiments = {}
    for mzml_file in file_paths:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_file, exp)
        experiments[Path(mzml_file).name] = exp
    return experiments


def create_spectrum_lookup(experiments: dict) -> dict:
    """Create SpectrumLookup maps for loaded experiments.
    Args:
        experiments (dict): Dictionary mapping file names to MSExperiment objects.
    Returns:
        dict: Dictionary mapping experiment names to SpectrumLookup objects."""
    lookup_map = {}
    for exp_name, exp in experiments.items():
        spectrum_lookup = oms.SpectrumLookup()
        spectrum_lookup.readSpectra(exp, "scan=(?<SCAN>\\d+)")
        lookup_map[Path(exp_name).stem] = spectrum_lookup
    return lookup_map


def filter_psms_by_peptide(psm_list: PSMList, target_peptide: str) -> PSMList:
    """Filter PSMs by target peptide sequence.
    Args:
        psm_list (PSMList): List of PSMs to filter.
        target_peptide (str): Peptide sequence to filter by.
    Returns:
        PSMList: Filtered list of PSMs matching the target peptide.
    """
    filtered_psms = PSMList(psm_list=[])
    for psm in psm_list:
        unmodified_peptide = Peptidoform(str(psm.peptidoform)).sequence
        if unmodified_peptide == target_peptide:
            filtered_psms.append(psm)
    logger.info(f"Found {len(filtered_psms)} PSMs matching '{target_peptide}'")
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
    if not (hasattr(psm, "spectrum") and psm.spectrum is not None):
        exp_name = f"{psm.run}.mzML"
        if exp_name in experiments and psm.run in lookups:
            idx = lookups[psm.run].findByNativeID(psm.spectrum_id)
            psm.spectrum = experiments[exp_name][idx]
    return psm


# Annotation functions    
def generate_theoretical_spectrum(pf: Peptidoform, max_charge: int = 2) -> oms.MSSpectrum:
    """
    Generate a theoretical spectrum for a given peptide form.
    Args:
        pf (PeptideForm): The peptide form to generate the spectrum for.
        max_charge (int): The maximum charge state to consider.
    Returns:
        MSSpectrum: The generated theoretical spectrum.
    """
    gen = oms.TheoreticalSpectrumGenerator()
    params = gen.getParameters()
    for ion in ("a","b","y","c","x","z"): 
        params.setValue(f"add_{ion}_ions","true")
    params.setValue("add_metainfo","true")
    gen.setParameters(params)
    spec = oms.MSSpectrum()
    gen.getSpectrum(spec, oms.AASequence(str(pf).split('/')[0].replace('[','(').replace(']',')')), 
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
        aligner.getSpectrumAlignment(pairs, exp_spec, theo_spec)
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
    theo_idxed = theo_df.reset_index().rename(columns={"index":"theo_idx"})
    rows = [{"exp_idx":e, "theo_idx":t} for e,t in pairs]
    df = pd.DataFrame(rows).merge(theo_idxed, on="theo_idx", how="left")
    df["exp_mz"] = df.exp_idx.map(lambda i: exp_spec[i].getMZ())
    df["exp_intensity"] = df.exp_idx.map(lambda i: exp_spec[i].getIntensity())
    return df


# Helper functions to create plots
def assign_ion_color(ion_name, color_map: dict) -> str:
    """Map ion name to color based on ion type"""
    if ion_name == "unassigned":
        return color_map.get("unannotated", "#666666")
    for ion_type in INCLUDED_IONS:
        if ion_name.startswith(ion_type):
            return color_map.get(ion_type, color_map.get("other", "#AAAAAA"))
    return color_map.get("other", "#AAAAAA")


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
    
    # Add metadata for hover and group by ion family
    peak_df["spectrum_type"] = legend_name
    peak_df["ion_family"] = peak_df["annotation"].apply(
        lambda x: "unassigned" if x == "unassigned" else (x[0] if x else "other"))

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
                      title="Mirror Plot", upper_label="", lower_label="", **kwargs):
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
        width (int): Plot width
        height (int): Plot height
    
    Returns:
        figure: Bokeh figure object with mirror plot
    """

    # Create plot with hover
    plot = setup_bokeh_plot(title, **kwargs)
    
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

    if not (upper_spectrum_data.empty and lower_spectrum_data.empty):
        all_mz = list(upper_spectrum_data["mz"]) + list(lower_spectrum_data["mz"])
        if all_mz:
            plot.line([min(all_mz), max(all_mz)], [0, 0], line_color="gray", line_width=1, alpha=0.5)

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


def predict_koina_ions(peptide_seq: str, charge: int, model: str = "ms2pip_timsTOF2024") -> pd.DataFrame:
    """Predict fragment ion intensities using Koina.
    Args:
        peptide_seq (str): The peptide sequence in UniMod format.
        charge (int): The precursor charge state.
        model (str): The model to use for prediction.
    Returns:
        pd.DataFrame: DataFrame with columns ['mz', 'intensity', 'annotation']
    """
    unimod_seq = convert_peptide_to_unimod(peptide_seq)
    koina = Koina(model)
    predictions = koina.predict(pd.DataFrame({"peptide_sequences": [unimod_seq],
                                              "precursor_charges": [charge]}))
    # Remove byte string artifacts from annotations
    predictions["annotation"] = predictions["annotation"].astype(str).str.replace("^b'|'$","", regex=True)
    
    return pd.DataFrame({'mz': predictions['mz'].astype(float),
                         'intensity': predictions['intensities'].astype(float),
                         'annotation': predictions['annotation']
    })


def create_exp_vs_pred_plot(psm, aligned_df, model=MODEL):
    """
    Create a mirror plot comparing experimental vs predicted spectrum for a PSM.
    Args:
        psm: PSM object with spectrum data
        aligned_df: DataFrame with aligned experimental data
        model: Koina model to use for predictions
    Returns:
        Bokeh figure object
    """
    # Get experimental data as dataframe
    exp_df = extract_psm_spectrum(psm, aligned_df)

    # Get predicted data using predict_koina_ions
    pred_df = predict_koina_ions(str(psm.peptidoform), psm.peptidoform.precursor_charge, model)
    
    title=f"Exp vs Pred: \n {psm.run}({psm.peptidoform})"

    return create_mirror_plot(exp_df, pred_df, title=title, upper_colors=EXP_COLORS, lower_colors=PRED_COLORS)