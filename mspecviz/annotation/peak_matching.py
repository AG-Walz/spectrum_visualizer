"""Match experimental peaks to theoretical fragments via PyOpenMS ``SpectrumAlignment``.

Produces the annotated DataFrames that feed mirror plots, comparison
metrics, and visualization frames; also handles Koina spectra whose
annotations are pre-attached.
"""
import logging
from typing import Dict, List

import numpy as np
import pandas as pd
import pyopenms as oms

from mspecviz.annotation.peak_theory import generate_theoretical_fragments
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.ion_utils import get_ion_type
from mspecviz.common.normalization import normalize_intensities

logger = logging.getLogger(__name__)


def annotate_spectrum(
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    theo_fragments: pd.DataFrame,
    mz_tolerance: float = 0.02
) -> pd.DataFrame:
    """Align experimental peaks to ``theo_fragments`` within ``mz_tolerance`` Da.

    Returns one row per experimental peak with ``mz, intensity, annotation,
    charge, is_annotated`` — unmatched peaks keep ``annotation='unannotated'``.
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
        except RuntimeError as e:
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


def _postprocess_annotated(annotated, ion_types, normalize, remove_precursor_peaks):
    """Filter unwanted ions, remove precursor peaks, and normalize intensities."""
    annotated["ion_type"] = annotated["annotation"].apply(get_ion_type)

    # Mask ions not in selected ion_types as unannotated
    unwanted = ~annotated["ion_type"].isin(ion_types) & (annotated["ion_type"] != "unannotated")
    annotated.loc[unwanted, "ion_type"] = "unannotated"
    annotated.loc[unwanted, "is_annotated"] = False

    if remove_precursor_peaks:
        precursor_mask = annotated["annotation"].apply(
            lambda ann: bool(ann and ann != "unannotated" and
                             (ann.lower().startswith("[m") or "precursor" in ann.lower()))
        )
        annotated = annotated[~precursor_mask].reset_index(drop=True)

    if normalize:
        annotated = normalize_intensities(annotated)

    return annotated


def prepare_spectrum_for_plot(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform: str,
    charge: int,
    ion_types: List[str] = None,
    mz_tolerance: float = 0.02,
    normalize: bool = True,
    add_neutral_losses: bool = False,
    add_first_prefix_ion: bool = True,
    remove_precursor_peaks: bool = False
) -> pd.DataFrame:
    """Annotate a spectrum and return a DataFrame ready for plotting."""
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)

    theo_frags = generate_theoretical_fragments(
        peptidoform, charge=charge, ion_types=ion_types,
        add_neutral_losses=add_neutral_losses, add_first_prefix_ion=add_first_prefix_ion,
    )
    annotated = annotate_spectrum(mz, intensity, theo_frags, mz_tolerance=mz_tolerance)
    return _postprocess_annotated(annotated, ion_types, normalize, remove_precursor_peaks)


def _parse_koina_charge(annotation: str) -> int:
    """Parse Koina annotation charge: ``y3+1`` → 1, ``b5+2`` → 2, ``y7++`` → 2."""
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
    koina_frame: pd.DataFrame,
    ion_types: List[str] = None,
    normalize: bool = True,
    remove_precursor_peaks: bool = False,
    **_kwargs,
) -> pd.DataFrame:
    """Prepare a Koina-predicted spectrum for plotting, preserving native annotations."""
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)

    annotated = pd.DataFrame({
        "mz": koina_frame["mz"].values,
        "intensity": koina_frame["intensity"].values.astype(float),
        "annotation": koina_frame["annotation"].values,
        "charge": koina_frame["annotation"].apply(_parse_koina_charge).values,
        "is_annotated": True,
    })
    return _postprocess_annotated(annotated, ion_types, normalize, remove_precursor_peaks)


def calculate_annotation_stats(annotated_frame: pd.DataFrame) -> Dict:
    """Return ``{total_peaks, annotated_peaks, coverage, ion_type_counts}`` for an annotated frame."""
    total = len(annotated_frame)
    annotated = annotated_frame["is_annotated"].sum() if "is_annotated" in annotated_frame else 0

    ion_counts = {}
    if "ion_type" in annotated_frame:
        counts = annotated_frame[annotated_frame["ion_type"] != "unannotated"]["ion_type"].value_counts()
        ion_counts = counts.to_dict()

    return {
        "total_peaks": total,
        "annotated_peaks": annotated,
        "coverage": annotated / total if total > 0 else 0,
        "ion_type_counts": ion_counts
    }
