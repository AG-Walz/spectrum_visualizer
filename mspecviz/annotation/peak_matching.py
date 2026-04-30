"""Annotate experimental + Koina spectra into the DataFrame schema used downstream.

Experimental annotation is delegated to ``spectrum_utils.MsmsSpectrum.annotate_proforma``
which generates theoretical fragments and aligns them in one call; we then
convert the per-peak ``PeakInterpretation`` objects back into the project's
DataFrame schema. Koina-predicted spectra arrive with native string annotations
that we keep verbatim — no theoretical step is needed.

Output schema for both paths:
    mz (float), intensity (float), annotation (str), ion_type (str),
    charge (int), is_annotated (bool)
"""
import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from spectrum_utils.spectrum import MsmsSpectrum

from mspecviz.common.defaults import DEFAULTS

logger = logging.getLogger(__name__)

_EMPTY_COLUMNS = ["mz", "intensity", "annotation", "ion_type", "charge", "is_annotated"]

# Per spectrum_utils' SUPPORTED_IONS = "?abcxyzIm_prf", the leading character of
# FragmentAnnotation.ion_type identifies the ion class. Koina annotation strings
# follow the same convention ("y3+1", "b5++"). The six fragment letters keep
# their letter so user filters like ['b', 'y'] match the dataframe column;
# uppercase 'I' and lowercase 'p' map to readable bucket names.
_BUCKET_FROM_PREFIX = {
    "a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z",
    "I": "immonium", "p": "precursor",
}


def get_ion_type(annotation: str) -> str:
    """Map a spectrum_utils ion_type or Koina annotation to its colour-bucket key."""
    return _BUCKET_FROM_PREFIX.get(annotation[:1], "unannotated") if annotation else "unannotated"


def _ion_types_to_string(ion_types: List[str]) -> str:
    """Pack ion_types list (e.g. ['b', 'y']) into spectrum_utils' single-string form ('by')."""
    return "".join(t for t in ion_types if t in ("a", "b", "c", "x", "y", "z"))


def _empty_annotated_frame(mz, intensity) -> pd.DataFrame:
    """Return a DataFrame with all peaks marked unannotated — used as fallback."""
    n = len(mz)
    return pd.DataFrame({
        "mz": np.asarray(mz, dtype=float),
        "intensity": np.asarray(intensity, dtype=float),
        "annotation": ["unannotated"] * n,
        "ion_type": ["unannotated"] * n,
        "charge": [0] * n,
        "is_annotated": [False] * n,
    })


def _spec_to_dataframe(spec: MsmsSpectrum) -> pd.DataFrame:
    """Convert an annotated MsmsSpectrum to the project's DataFrame schema.

    spectrum_utils returns ``PeakInterpretation`` per peak, each carrying a list
    of ``FragmentAnnotation`` (multiple matches possible). We keep the first match.
    """
    annotations: List[str] = []
    ion_types: List[str] = []
    charges: List[int] = []
    flags: List[bool] = []

    for peak_ann in spec.annotation:
        fa_list = peak_ann.fragment_annotations if peak_ann else []
        if not fa_list:
            annotations.append("unannotated")
            ion_types.append("unannotated")
            charges.append(0)
            flags.append(False)
            continue
        fa = fa_list[0]
        annotations.append(f"{fa.ion_type}+{fa.charge}")
        ion_types.append(get_ion_type(fa.ion_type))
        charges.append(int(fa.charge))
        flags.append(True)

    return pd.DataFrame({
        "mz": np.asarray(spec.mz, dtype=float),
        "intensity": np.asarray(spec.intensity, dtype=float),
        "annotation": annotations,
        "ion_type": ion_types,
        "charge": charges,
        "is_annotated": flags,
    })


def annotate_spectrum(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform: str,
    charge: int,
    mz_tolerance: float = 0.02,
    ion_types: Optional[List[str]] = None,
    add_neutral_losses: bool = False,
) -> pd.DataFrame:
    """Annotate experimental peaks against theoretical fragments of ``peptidoform``.

    Returns one row per experimental peak in the standard schema.
    """
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)
    if len(mz) == 0:
        return pd.DataFrame(columns=_EMPTY_COLUMNS)

    ion_types_str = _ion_types_to_string(ion_types)
    if not ion_types_str:
        return _empty_annotated_frame(mz, intensity)

    spec = MsmsSpectrum(
        identifier="experimental",
        precursor_mz=0.0,
        precursor_charge=int(charge),
        mz=np.asarray(mz, dtype=float),
        intensity=np.asarray(intensity, dtype=float),
    )
    try:
        spec.annotate_proforma(
            peptidoform.split("/")[0],   # strip ProForma charge suffix if present
            fragment_tol_mass=mz_tolerance,
            fragment_tol_mode="Da",
            ion_types=ion_types_str,
            neutral_losses=add_neutral_losses,
        )
    except (ValueError, KeyError) as e:
        logger.warning(f"annotate_proforma failed for '{peptidoform}': {e}")
        return _empty_annotated_frame(mz, intensity)

    return _spec_to_dataframe(spec)


def _filter_and_normalize(
    annotated: pd.DataFrame, ion_types: List[str], normalize: bool,
) -> pd.DataFrame:
    """Mask peaks whose ion_type isn't in ``ion_types``; optionally scale to max=1."""
    if annotated.empty:
        return annotated

    keep_letters = set(ion_types) | {"unannotated"}
    unwanted = ~annotated["ion_type"].isin(keep_letters)
    if unwanted.any():
        annotated.loc[unwanted, "annotation"] = "unannotated"
        annotated.loc[unwanted, "ion_type"] = "unannotated"
        annotated.loc[unwanted, "is_annotated"] = False
        annotated.loc[unwanted, "charge"] = 0

    if normalize:
        max_intensity = annotated["intensity"].max()
        if max_intensity and max_intensity > 0:
            annotated["intensity"] = annotated["intensity"] / max_intensity

    return annotated


def prepare_spectrum_for_plot(
    mz: np.ndarray,
    intensity: np.ndarray,
    peptidoform: str,
    charge: int,
    ion_types: Optional[List[str]] = None,
    mz_tolerance: float = 0.02,
    normalize: bool = True,
    add_neutral_losses: bool = False,
) -> pd.DataFrame:
    """Annotate a spectrum and return a DataFrame ready for plotting / metrics."""
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)
    annotated = annotate_spectrum(
        mz, intensity, peptidoform, charge,
        mz_tolerance=mz_tolerance,
        ion_types=ion_types,
        add_neutral_losses=add_neutral_losses,
    )
    return _filter_and_normalize(annotated, ion_types, normalize)


def _parse_koina_charge(annotation: str) -> int:
    """Parse Koina annotation charge: ``y3+1`` → 1, ``b5+2`` → 2, ``y7++`` → 2."""
    if not annotation:
        return 0
    plus_count = len(annotation) - len(annotation.rstrip("+"))
    if plus_count >= 2:
        return plus_count
    parts = annotation.rsplit("+", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return int(parts[1])
    return 1


def prepare_koina_spectrum_for_plot(
    koina_frame: pd.DataFrame,
    ion_types: Optional[List[str]] = None,
    normalize: bool = True,
    **_kwargs,
) -> pd.DataFrame:
    """Prepare a Koina-predicted spectrum for plotting; preserves native annotations."""
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)

    annotated = pd.DataFrame({
        "mz": koina_frame["mz"].values,
        "intensity": koina_frame["intensity"].values.astype(float),
        "annotation": koina_frame["annotation"].values,
        "charge": koina_frame["annotation"].apply(_parse_koina_charge).values,
        "is_annotated": True,
    })
    annotated["ion_type"] = annotated["annotation"].apply(get_ion_type)
    return _filter_and_normalize(annotated, ion_types, normalize)


def calculate_annotation_stats(annotated_frame: pd.DataFrame) -> Dict:
    """Return ``{total_peaks, annotated_peaks, coverage, ion_type_counts}``."""
    total = len(annotated_frame)
    annotated = int(annotated_frame["is_annotated"].sum()) if "is_annotated" in annotated_frame else 0

    ion_counts: Dict[str, int] = {}
    if "ion_type" in annotated_frame:
        counts = annotated_frame[annotated_frame["ion_type"] != "unannotated"]["ion_type"].value_counts()
        ion_counts = counts.to_dict()

    return {
        "total_peaks": total,
        "annotated_peaks": annotated,
        "coverage": annotated / total if total > 0 else 0,
        "ion_type_counts": ion_counts,
    }
