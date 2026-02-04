"""
Peak annotation for MS/MS spectra.

Generates theoretical fragment ions and annotates experimental peaks.
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyopenms as oms
from psm_utils import Peptidoform

logger = logging.getLogger(__name__)

# Default ion types for fragmentation
DEFAULT_ION_TYPES = ["b", "y"]
ALL_ION_TYPES = ["a", "b", "c", "x", "y", "z", "precursor", "immonium"]


def peptidoform_to_aasequence(peptidoform_str: str) -> oms.AASequence:
    """
    Convert a peptidoform string to PyOpenMS AASequence.

    Handles bracket notation conversion: [] -> ()
    """
    # Remove charge if present (after /)
    pf_clean = peptidoform_str.split('/')[0]

    # Convert bracket notation for modifications
    # psm_utils uses [] for mods, PyOpenMS uses ()
    pf_openms = pf_clean.replace('[', '(').replace(']', ')')

    try:
        return oms.AASequence.fromString(pf_openms)
    except Exception as e:
        logger.warning(f"Could not parse peptidoform '{peptidoform_str}': {e}")
        # Try with just the sequence
        try:
            pf = Peptidoform(peptidoform_str)
            return oms.AASequence.fromString(pf.sequence)
        except Exception:
            raise ValueError(f"Cannot convert peptidoform: {peptidoform_str}")


def generate_theoretical_fragments(
    peptidoform: str,
    charge: int = 2,
    ion_types: List[str] = None,
    max_fragment_charge: int = 2
) -> pd.DataFrame:
    """
    Generate theoretical fragment ions for a peptide.

    Args:
        peptidoform: Peptidoform string (with or without modifications)
        charge: Precursor charge state
        ion_types: List of ion types to generate (default: ["b", "y"])
        max_fragment_charge: Maximum charge for fragment ions

    Returns:
        DataFrame with columns: mz, intensity (relative), annotation, charge
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    try:
        aa_seq = peptidoform_to_aasequence(peptidoform)
    except ValueError as e:
        logger.error(f"Failed to generate fragments: {e}")
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Configure theoretical spectrum generator
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()

    # Set ion types
    for ion in ["a", "b", "c", "x", "y", "z"]:
        params.setValue(f"add_{ion}_ions", "true" if ion in ion_types else "false")

    params.setValue("add_metainfo", "true")
    params.setValue("add_losses", "true" if "loss" in ion_types else "false")
    params.setValue("add_precursor_peaks", "true" if "precursor" in ion_types else "false")
    params.setValue("add_first_prefix_ion", "true")
    params.setValue("add_abundant_immonium_ions", "true" if "immonium" in ion_types else "false")

    tsg.setParameters(params)

    # Generate spectrum
    theo_spec = oms.MSSpectrum()
    tsg.getSpectrum(theo_spec, aa_seq, 1, max_fragment_charge)

    if theo_spec.size() > 1:
        theo_spec.sortByPosition()

    # Convert to DataFrame
    if theo_spec.size() == 0:
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Extract annotations from StringDataArrays
    string_arrays = theo_spec.getStringDataArrays()
    int_arrays = theo_spec.getIntegerDataArrays()

    annotations = []
    if string_arrays:
        annotations = [n.decode() if isinstance(n, bytes) else str(n) for n in string_arrays[0]]
    else:
        annotations = [""] * theo_spec.size()

    charges = []
    if int_arrays:
        charges = list(int_arrays[0])
    else:
        charges = [1] * theo_spec.size()

    data = {
        "mz": [peak.getMZ() for peak in theo_spec],
        "intensity": [peak.getIntensity() for peak in theo_spec],
        "annotation": annotations,
        "charge": charges
    }

    return pd.DataFrame(data)


def annotate_spectrum(
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    theo_fragments: pd.DataFrame,
    mz_tolerance: float = 0.02
) -> pd.DataFrame:
    """
    Annotate experimental spectrum peaks with theoretical fragment matches.

    Args:
        exp_mz: Experimental m/z values
        exp_intensity: Experimental intensity values
        theo_fragments: DataFrame of theoretical fragments
        mz_tolerance: Tolerance for matching in Da

    Returns:
        DataFrame with columns: mz, intensity, annotation, charge, is_annotated
    """
    if len(exp_mz) == 0:
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge", "is_annotated"])

    # Create experimental spectrum object
    exp_spec = oms.MSSpectrum()
    for mz, intensity in zip(exp_mz, exp_intensity):
        peak = oms.Peak1D()
        peak.setMZ(mz)
        peak.setIntensity(intensity)
        exp_spec.push_back(peak)
    exp_spec.sortByPosition()

    # Create theoretical spectrum object
    theo_spec = oms.MSSpectrum()
    for _, row in theo_fragments.iterrows():
        peak = oms.Peak1D()
        peak.setMZ(row['mz'])
        peak.setIntensity(row['intensity'])
        theo_spec.push_back(peak)
    theo_spec.sortByPosition()

    # Perform alignment
    pairs = []
    if theo_spec.size() > 0 and exp_spec.size() > 0:
        try:
            aligner = oms.SpectrumAlignment()
            params = aligner.getParameters()
            params.setValue("tolerance", mz_tolerance)
            params.setValue("is_relative_tolerance", "false")
            aligner.setParameters(params)
            aligner.getSpectrumAlignment(pairs, theo_spec, exp_spec)
        except Exception as e:
            logger.warning(f"Spectrum alignment failed: {e}")

    # Build result DataFrame
    # Start with all experimental peaks as unannotated
    result_data = {
        "mz": list(exp_mz),
        "intensity": list(exp_intensity),
        "annotation": ["unannotated"] * len(exp_mz),
        "charge": [0] * len(exp_mz),
        "is_annotated": [False] * len(exp_mz)
    }

    # Update annotated peaks
    theo_list = theo_fragments.to_dict('records')
    for theo_idx, exp_idx in pairs:
        if exp_idx < len(result_data["mz"]) and theo_idx < len(theo_list):
            theo_row = theo_list[theo_idx]
            result_data["annotation"][exp_idx] = theo_row.get("annotation", "matched")
            result_data["charge"][exp_idx] = theo_row.get("charge", 0)
            result_data["is_annotated"][exp_idx] = True

    return pd.DataFrame(result_data)


def get_ion_type(annotation: str) -> str:
    """
    Extract ion type from annotation string.

    Examples:
        "b3++" -> "b"
        "y7" -> "y"
        "[M+2H]2+" -> "precursor"
        "unannotated" -> "unannotated"
    """
    if not annotation or annotation == "unannotated":
        return "unannotated"

    annotation_lower = annotation.lower()

    # Check for precursor
    if annotation_lower.startswith("[m") or "precursor" in annotation_lower:
        return "precursor"

    # Check for standard ions
    for ion in ["a", "b", "c", "x", "y", "z"]:
        if annotation_lower.startswith(ion) and len(annotation) > 1:
            next_char = annotation[1]
            if next_char.isdigit():
                return ion

    # Check for immonium
    if "immonium" in annotation_lower or annotation.startswith("I"):
        return "immonium"

    return "unannotated"


def prepare_spectrum_for_plot(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform: str,
    charge: int,
    ion_types: List[str] = None,
    mz_tolerance: float = 0.02,
    normalize: bool = True
) -> pd.DataFrame:
    """
    Prepare a spectrum for plotting with annotations.

    Args:
        mz: m/z array
        intensity: Intensity array
        peptidoform: Peptidoform string
        charge: Precursor charge
        ion_types: Ion types to annotate
        mz_tolerance: Matching tolerance
        normalize: Whether to normalize intensities to [0, 1]

    Returns:
        DataFrame ready for plotting with: mz, intensity, annotation, ion_type, charge
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    # Generate theoretical fragments
    theo_frags = generate_theoretical_fragments(
        peptidoform,
        charge=charge,
        ion_types=ion_types
    )

    # Annotate spectrum
    annotated = annotate_spectrum(
        mz, intensity, theo_frags,
        mz_tolerance=mz_tolerance
    )

    # Add ion type column
    annotated["ion_type"] = annotated["annotation"].apply(get_ion_type)

    # Normalize if requested
    if normalize and len(annotated) > 0:
        max_int = annotated["intensity"].max()
        if max_int > 0:
            annotated["intensity"] = annotated["intensity"] / max_int

    return annotated


def calculate_annotation_stats(annotated_df: pd.DataFrame) -> Dict:
    """
    Calculate annotation statistics for a spectrum.

    Returns dict with: total_peaks, annotated_peaks, coverage, ion_type_counts
    """
    total = len(annotated_df)
    annotated = annotated_df["is_annotated"].sum() if "is_annotated" in annotated_df else 0

    ion_counts = {}
    if "ion_type" in annotated_df:
        counts = annotated_df[annotated_df["ion_type"] != "unannotated"]["ion_type"].value_counts()
        ion_counts = counts.to_dict()

    return {
        "total_peaks": total,
        "annotated_peaks": annotated,
        "coverage": annotated / total if total > 0 else 0,
        "ion_type_counts": ion_counts
    }
