"""
Peak annotation for MS/MS spectra.

Generates theoretical fragment ions and annotates experimental peaks.
"""
import functools
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyopenms as oms
from psm_utils import Peptidoform
from spectrum_utils.spectrum import MsmsSpectrum

from src.common.ion_utils import get_ion_type  # noqa: F401

logger = logging.getLogger(__name__)

# Default ion types for fragmentation
DEFAULT_ION_TYPES = ["b", "y"]
ALL_ION_TYPES = ["a", "b", "c", "x", "y", "z", "precursor", "immonium"]


def peptidoform_to_aasequence(peptidoform_str: str) -> oms.AASequence:
    """
    Convert a ProForma peptidoform string to PyOpenMS AASequence.

    Parses via psm_utils Peptidoform and reconstructs in OpenMS notation,
    correctly handling terminal modifications (e.g. [Acetyl]-PEPTIDE becomes
    .(Acetyl)PEPTIDE rather than the invalid (Acetyl)-PEPTIDE).
    """
    pf_clean = peptidoform_str.split('/')[0]
    pf = Peptidoform(pf_clean)

    # Rebuild in OpenMS parenthesis notation from the parsed representation
    parts = []
    if pf.properties.get("n_term"):
        for m in pf.properties["n_term"]:
            parts.append(f".({m.name})")
    for aa, residue_mods in pf.parsed_sequence:
        parts.append(aa)
        if residue_mods:
            for m in residue_mods:
                parts.append(f"({m.name})")
    if pf.properties.get("c_term"):
        for m in pf.properties["c_term"]:
            parts.append(f".({m.name})")
    oms_str = "".join(parts)

    try:
        return oms.AASequence.fromString(oms_str)
    except Exception as e:
        logger.warning(f"Could not parse peptidoform '{peptidoform_str}': {e}")
        try:
            return oms.AASequence.fromString(pf.sequence)
        except Exception:
            raise ValueError(f"Cannot convert peptidoform: {peptidoform_str}")


@functools.lru_cache(maxsize=1024)
def _generate_fragments_impl(
    peptidoform: str,
    charge: int,
    ion_types: tuple,          # tuple (hashable) — converted from list by caller
    max_fragment_charge: int,
    add_losses: bool,
    add_first_prefix_ion: bool,
) -> pd.DataFrame:
    """
    Cached implementation of theoretical fragment generation.

    Returns the same DataFrame object on repeated calls with identical
    arguments — callers must not mutate the returned DataFrame.
    """
    try:
        aa_seq = peptidoform_to_aasequence(peptidoform)
    except ValueError as e:
        logger.error(f"Failed to generate fragments: {e}")
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Configure theoretical spectrum generator
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()

    ion_types_list = list(ion_types)
    for ion in ["a", "b", "c", "x", "y", "z"]:
        params.setValue(f"add_{ion}_ions", "true" if ion in ion_types_list else "false")

    params.setValue("add_metainfo", "true")
    params.setValue("add_losses", "true" if add_losses else "false")
    params.setValue("add_precursor_peaks", "true")
    params.setValue("add_all_precursor_charges", "true")
    params.setValue("add_first_prefix_ion", "true" if add_first_prefix_ion else "false")
    params.setValue("add_abundant_immonium_ions", "true" if "immonium" in ion_types_list else "false")

    tsg.setParameters(params)

    # Generate spectrum
    theo_spec = oms.MSSpectrum()
    tsg.getSpectrum(theo_spec, aa_seq, 1, max_fragment_charge)

    if theo_spec.size() > 1:
        theo_spec.sortByPosition()

    if theo_spec.size() == 0:
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Use batch API — much faster than iterating Peak1D objects
    mzs, ints = theo_spec.get_peaks()

    # Extract annotations from StringDataArrays
    string_arrays = theo_spec.getStringDataArrays()
    int_arrays = theo_spec.getIntegerDataArrays()

    if string_arrays:
        annotations = [n.decode() if isinstance(n, bytes) else str(n) for n in string_arrays[0]]
    else:
        annotations = [""] * theo_spec.size()

    if int_arrays:
        charges = list(int_arrays[0])
    else:
        charges = [1] * theo_spec.size()

    return pd.DataFrame({
        "mz": mzs,
        "intensity": ints,
        "annotation": annotations,
        "charge": charges,
    })


def generate_theoretical_fragments(
    peptidoform: str,
    charge: int = 2,
    ion_types: List[str] = None,
    max_fragment_charge: int = 2,
    add_losses: bool = False,
    add_first_prefix_ion: bool = True,
) -> pd.DataFrame:
    """
    Generate theoretical fragment ions for a peptide.

    Results are cached by (peptidoform, charge, ion_types, add_losses,
    add_first_prefix_ion) using functools.lru_cache — identical calls within
    the same process (e.g. many query spectra for the same peptidoform) pay
    the PyOpenMS cost only once.

    Args:
        peptidoform: Peptidoform string (with or without modifications)
        charge: Precursor charge state
        ion_types: List of ion types to generate (default: ["b", "y"])
        max_fragment_charge: Maximum charge for fragment ions
        add_losses: Whether to add neutral loss peaks (H2O, NH3)
        add_first_prefix_ion: Whether to add first prefix ion (a1, b1)

    Returns:
        DataFrame with columns: mz, intensity (relative), annotation, charge
        Do NOT modify the returned DataFrame — it is shared across callers.
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    return _generate_fragments_impl(
        peptidoform,
        charge,
        tuple(ion_types),
        max_fragment_charge,
        add_losses,
        add_first_prefix_ion,
    )


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

    # Build experimental MSSpectrum using batch API (faster than push_back loop)
    exp_spec = oms.MSSpectrum()
    exp_spec.set_peaks([list(exp_mz), list(exp_intensity)])
    exp_spec.sortByPosition()

    # Build theoretical MSSpectrum using batch API
    theo_spec = oms.MSSpectrum()
    theo_spec.set_peaks([
        theo_fragments["mz"].tolist(),
        theo_fragments["intensity"].tolist(),
    ])
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


def prepare_spectrum_for_plot(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform: str,
    charge: int,
    ion_types: List[str] = None,
    mz_tolerance: float = 0.02,
    normalize: bool = True,
    add_losses: bool = False,
    add_first_prefix_ion: bool = True,
    remove_precursor_peaks: bool = False
) -> pd.DataFrame:
    """
    Prepare a spectrum for plotting with annotations.

    Args:
        mz: m/z array
        intensity: Intensity array
        peptidoform: Peptidoform string
        charge: Precursor charge
        ion_types: Ion types to annotate (determines coloring)
        mz_tolerance: Matching tolerance
        normalize: Whether to normalize intensities to [0, 1]
        add_losses: Whether to include neutral loss peaks in theoretical spectrum
        add_first_prefix_ion: Whether to include first prefix ion (a1, b1)
        remove_precursor_peaks: Whether to remove precursor peaks from the result

    Returns:
        DataFrame ready for plotting with: mz, intensity, annotation, ion_type, charge
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    # Generate theoretical fragments — lru_cache avoids repeated PyOpenMS calls
    # for the same peptidoform+charge+ion_types combination.
    theo_frags = generate_theoretical_fragments(
        peptidoform,
        charge=charge,
        ion_types=ion_types,
        add_losses=add_losses,
        add_first_prefix_ion=add_first_prefix_ion
    )

    # Annotate spectrum
    annotated = annotate_spectrum(
        mz, intensity, theo_frags,
        mz_tolerance=mz_tolerance
    )

    # Add ion type column
    annotated["ion_type"] = annotated["annotation"].apply(get_ion_type)

    # If "precursor" not in ion_types, mark precursor peaks as unannotated (grey)
    if "precursor" not in ion_types:
        precursor_mask = annotated["ion_type"] == "precursor"
        annotated.loc[precursor_mask, "ion_type"] = "unannotated"
        annotated.loc[precursor_mask, "is_annotated"] = False

    # Remove precursor peaks if requested
    if remove_precursor_peaks:
        def is_precursor_annotation(ann):
            if not ann or ann == "unannotated":
                return False
            ann_lower = ann.lower()
            return ann_lower.startswith("[m") or "precursor" in ann_lower

        precursor_mask = annotated["annotation"].apply(is_precursor_annotation)
        annotated = annotated[~precursor_mask].reset_index(drop=True)

    # Normalize if requested (do this AFTER filtering so max is recalculated)
    if normalize and len(annotated) > 0:
        max_int = annotated["intensity"].max()
        if max_int > 0:
            annotated["intensity"] = annotated["intensity"] / max_int

    return annotated


def _parse_koina_charge(annotation: str) -> int:
    """
    Extract charge from a Koina annotation string.

    Koina formats: 'y3+1' (charge 1), 'b5+2' (charge 2), 'y7++' (charge 2).
    """
    if not annotation:
        return 0
    # Count trailing '+' signs (e.g. 'y3++' → charge 2)
    plus_count = len(annotation) - len(annotation.rstrip("+"))
    if plus_count >= 2:
        return plus_count
    # Try 'y3+1' or 'b5+2' format
    parts = annotation.rsplit("+", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return int(parts[1])
    return 1


def prepare_koina_spectrum_for_plot(
    koina_df: pd.DataFrame,
    ion_types: List[str] = None,
    normalize: bool = True,
    remove_precursor_peaks: bool = False,
    **_kwargs,
) -> pd.DataFrame:
    """
    Prepare a Koina-predicted spectrum for plotting, preserving native annotations.

    Unlike prepare_spectrum_for_plot which re-matches against PyOpenMS theoretical
    fragments, this function uses Koina's own AI-predicted annotations directly.

    Args:
        koina_df: DataFrame from predict_koina_spectrum() with columns [mz, intensity, annotation]
        ion_types: Ion types to highlight (others shown as unannotated)
        normalize: Whether to normalize intensities to [0, 1]
        remove_precursor_peaks: Whether to remove precursor peaks

    Returns:
        DataFrame ready for plotting with: mz, intensity, annotation, ion_type, charge, is_annotated
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    annotated = pd.DataFrame({
        "mz": koina_df["mz"].values,
        "intensity": koina_df["intensity"].values.astype(float),
        "annotation": koina_df["annotation"].values,
        "charge": koina_df["annotation"].apply(_parse_koina_charge).values,
        "is_annotated": True,
    })

    annotated["ion_type"] = annotated["annotation"].apply(get_ion_type)

    # Grey out ions not in the user's selected ion_types
    unwanted_mask = ~annotated["ion_type"].isin(ion_types) & (annotated["ion_type"] != "unannotated")
    annotated.loc[unwanted_mask, "ion_type"] = "unannotated"
    annotated.loc[unwanted_mask, "is_annotated"] = False

    if remove_precursor_peaks:
        precursor_mask = annotated["ion_type"] == "precursor"
        annotated = annotated[~precursor_mask].reset_index(drop=True)

    if normalize and len(annotated) > 0:
        max_int = annotated["intensity"].max()
        if max_int > 0:
            annotated["intensity"] = annotated["intensity"] / max_int

    return annotated


def msms_spectrum_from_loaded_spectrum(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform_str: str,
    precursor_mz: float,
    precursor_charge: int,
) -> "spectrum_utils.spectrum.MsmsSpectrum":
    """
    Build a spectrum_utils MsmsSpectrum from raw arrays.

    Entry point for spectrum_utils-based annotation (e.g. edge-case ProForma mods
    using UNIMOD:N notation that PyOpenMS cannot parse).

    Args:
        mz: m/z array
        intensity: intensity array
        peptidoform_str: ProForma peptidoform string (charge suffix stripped)
        precursor_mz: Precursor m/z
        precursor_charge: Precursor charge state

    Returns:
        spectrum_utils MsmsSpectrum ready for annotate_proforma()
    """
    proforma_str = peptidoform_str.split("/")[0]
    return MsmsSpectrum(
        identifier=proforma_str,
        precursor_mz=float(precursor_mz),
        precursor_charge=int(precursor_charge),
        mz=np.array(mz, dtype=np.float64),
        intensity=np.array(intensity, dtype=np.float32),
    )


def _su_annotation_to_str(ann) -> str:
    """Convert a spectrum_utils FragmentAnnotation to an MSpecViz annotation string."""
    try:
        ion_type = str(ann.ion_type)
        ion_index = getattr(ann, "ion_index", getattr(ann, "fragment_index", "?"))
        charge = int(getattr(ann, "charge", 1))
        nl = getattr(ann, "neutral_loss", None)
        charge_str = "+" * max(1, charge)
        base = f"{ion_type}{ion_index}"
        if nl:
            base += f"-{str(nl).lstrip('-')}"
        return base + charge_str
    except Exception:
        return "matched"


def annotate_with_spectrum_utils(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform_str: str,
    precursor_mz: float,
    precursor_charge: int,
    fragment_tol_mass: float = 0.02,
    fragment_tol_mode: str = "Da",
    ion_types: List[str] = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Annotate a spectrum using spectrum_utils annotate_proforma().

    Secondary annotation path for ProForma 2.0 / UNIMOD:N edge cases that
    PyOpenMS cannot handle.  Falls back to an unannotated DataFrame on error.

    Returns a DataFrame in the same format as the primary annotation pipeline:
    columns: mz, intensity, annotation, charge, is_annotated, ion_type
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    su_ion_types = "".join(t for t in ion_types if t in ("a", "b", "c", "x", "y", "z"))
    if not su_ion_types:
        su_ion_types = "by"

    proforma_str = peptidoform_str.split("/")[0]

    def _unannotated_df():
        int_arr = np.array(intensity, dtype=float)
        if normalize and int_arr.max() > 0:
            int_arr = int_arr / int_arr.max()
        return pd.DataFrame({
            "mz": np.array(mz, dtype=float),
            "intensity": int_arr,
            "annotation": ["unannotated"] * len(mz),
            "charge": [0] * len(mz),
            "is_annotated": [False] * len(mz),
            "ion_type": ["unannotated"] * len(mz),
        })

    try:
        spec = msms_spectrum_from_loaded_spectrum(
            mz, intensity, proforma_str, precursor_mz, precursor_charge
        )
        spec.annotate_proforma(
            proforma_str,
            fragment_tol_mass,
            fragment_tol_mode,
            ion_types=su_ion_types,
        )
    except Exception as e:
        logger.warning(f"spectrum_utils annotation failed for '{peptidoform_str}': {e}")
        return _unannotated_df()

    annotations_out = []
    charges_out = []
    is_annotated_out = []
    ion_types_out = []

    for ann in spec.annotation:
        if ann is None:
            annotations_out.append("unannotated")
            charges_out.append(0)
            is_annotated_out.append(False)
            ion_types_out.append("unannotated")
        else:
            ann_str = _su_annotation_to_str(ann)
            charge = int(getattr(ann, "charge", 1))
            ion_t = str(getattr(ann, "ion_type", ""))[0:1] or "unannotated"
            annotations_out.append(ann_str)
            charges_out.append(charge)
            is_annotated_out.append(True)
            ion_types_out.append(ion_t)

    int_arr = np.array(spec.intensity, dtype=float)
    if normalize and int_arr.max() > 0:
        int_arr = int_arr / int_arr.max()

    return pd.DataFrame({
        "mz": np.array(spec.mz, dtype=float),
        "intensity": int_arr,
        "annotation": annotations_out,
        "charge": charges_out,
        "is_annotated": is_annotated_out,
        "ion_type": ion_types_out,
    })


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
