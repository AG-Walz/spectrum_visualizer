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

import numpy as np
import pandas as pd
from spectrum_utils.spectrum import MsmsSpectrum

from core.common.defaults import DEFAULTS

logger = logging.getLogger(__name__)

# Per spectrum_utils' SUPPORTED_IONS 
_BUCKET_FROM_PREFIX = {
    "a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z",
    "I": "immonium", "p": "precursor", "m": "internal",
}


def get_ion_type(annotation: str) -> str:
    """Map a spectrum_utils ion_type or Koina annotation to its colour-bucket key."""
    return _BUCKET_FROM_PREFIX.get(annotation[:1], "unannotated")


def _spec_to_dataframe(spec: MsmsSpectrum) -> pd.DataFrame:
    """Convert an annotated MsmsSpectrum to the project's DataFrame schema. """
    # Allocate per-column accumulators sized to the number of peaks
    annotations: list[str] = []
    ion_types: list[str] = []
    charges: list[int] = []
    flags: list[bool] = []

    # Walk peaks: take first FragmentAnnotation, fall back to unannotated row
    for peak_ann in spec.annotation:
        fa_list = peak_ann.fragment_annotations if peak_ann else []
        if not fa_list:
            annotations.append("unannotated")
            ion_types.append("unannotated")
            charges.append(0)
            flags.append(False)
            continue
        fa = fa_list[0]
        loss_suffix = f"-{fa.neutral_loss}" if fa.neutral_loss else ""
        annotations.append(f"{fa.ion_type}{loss_suffix}+{fa.charge}")
        ion_types.append(get_ion_type(fa.ion_type))
        charges.append(int(fa.charge))
        flags.append(True)

    # Pack columns into the standard schema DataFrame
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
    ion_types: list[str] | None = None,
    add_neutral_losses: bool = False,
    precursor_mz: float = 0.0,
) -> pd.DataFrame:
    """Annotate experimental peaks against theoretical fragments of ``peptidoform``.
    Returns one row per experimental peak in the standard schema.

    *precursor_mz* must be set to the experimental precursor m/z for
    ``ion_types`` to include ``"precursor"`` — spectrum_utils needs a real
    value to match precursor peaks in the spectrum. A 0.0 fallback keeps
    backbone-ion annotation working when the caller doesn't have it.
    """
    # Defaults: empty-input + ion-type string + fallback-frame builder
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)
    # Pack UI ion-type names into spectrum_utils' single-letter form:
    # backbone (a/b/c/x/y/z) and the three "special" ions it accepts —
    # m = internal, p = precursor, I = immonium.
    _SU_LETTERS = ("a", "b", "c", "x", "y", "z")
    _SU_SPECIAL = {"internal": "m", "precursor": "p", "immonium": "I"}
    ion_types_string = "".join(
        [t for t in ion_types if t in _SU_LETTERS]
        + [letter for name, letter in _SU_SPECIAL.items() if name in ion_types]
    )
    n = len(mz)

    def _fallback_unannotated_frame() -> pd.DataFrame:
        return pd.DataFrame({
            "mz": np.asarray(mz, dtype=float),
            "intensity": np.asarray(intensity, dtype=float),
            "annotation": ["unannotated"] * n,
            "ion_type": ["unannotated"] * n,
            "charge": [0] * n,
            "is_annotated": [False] * n,
        })

    # Fast exit: empty input or no requested ion letters
    if n == 0 or not ion_types_string:
        return _fallback_unannotated_frame()

    # Build the spectrum_utils MsmsSpectrum carrier
    spec = MsmsSpectrum(
        identifier="experimental",
        precursor_mz=float(precursor_mz),
        precursor_charge=int(charge),
        mz=np.asarray(mz, dtype=float),
        intensity=np.asarray(intensity, dtype=float),
    )

    # Annotate via ProForma; fall back to unannotated frame on parse/match failure
    try:
        spec.annotate_proforma(
            peptidoform.split("/")[0],   # strip ProForma charge suffix if present
            fragment_tol_mass=mz_tolerance,
            fragment_tol_mode="Da",
            ion_types=ion_types_string,
            neutral_losses=add_neutral_losses,
        )
    except (ValueError, KeyError) as error:
        logger.warning(f"annotate_proforma failed for '{peptidoform}': {error}")
        return _fallback_unannotated_frame()

    return _spec_to_dataframe(spec)


def _filter_and_normalize(
    annotated: pd.DataFrame, ion_types: list[str], normalize: bool,
) -> pd.DataFrame:
    """Mask peaks whose ion_type isn't in ``ion_types`` and optionally scale to max=1."""
    # Fast exit: nothing to filter or scale
    if annotated.empty:
        return annotated

    # Mask peaks outside the requested ion_types — reset annotation/type/charge
    keep_letters = set(ion_types) | {"unannotated"}
    unwanted = ~annotated["ion_type"].isin(keep_letters)
    if unwanted.any():
        annotated.loc[unwanted, "annotation"] = "unannotated"
        annotated.loc[unwanted, "ion_type"] = "unannotated"
        annotated.loc[unwanted, "is_annotated"] = False
        annotated.loc[unwanted, "charge"] = 0

    # Optionally rescale intensities so the brightest peak equals 1.0
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
    ion_types: list[str] | None = None,
    mz_tolerance: float = 0.02,
    normalize: bool = True,
    add_neutral_losses: bool = False,
    precursor_mz: float = 0.0,
) -> pd.DataFrame:
    """Annotate a spectrum and return a DataFrame ready for plotting / metrics."""
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)
    annotated = annotate_spectrum(
        mz, intensity, peptidoform, charge,
        mz_tolerance=mz_tolerance,
        ion_types=ion_types,
        add_neutral_losses=add_neutral_losses,
        precursor_mz=precursor_mz,
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
    ion_types: list[str] | None = None,
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
