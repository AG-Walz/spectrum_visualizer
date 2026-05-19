"""Match experimental or Koina peaks to theoretical fragments of a peptidoform (called by ``core.pipeline.annotation_steps``).

Experimental path delegates to ``spectrum_utils.annotate_proforma``; Koina keeps native labels."""
import logging

import numpy as np
import pandas as pd
from spectrum_utils.spectrum import MsmsSpectrum

from core.defaults import DEFAULTS

logger = logging.getLogger(__name__)

# First-char of an annotation -> ion-type bucket used for plot colouring
_BUCKET_FROM_PREFIX = {
    "a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z",
    "I": "immonium", "p": "precursor", "m": "internal",
}


def get_ion_type(annotation: str) -> str:
    """Return the colour-bucket key (``b``, ``y``, ``immonium``, ``precursor``...) for an annotation string."""
    return _BUCKET_FROM_PREFIX.get(annotation[:1], "unannotated")


def _spec_to_dataframe(spec: MsmsSpectrum) -> pd.DataFrame:
    """Convert an annotated ``MsmsSpectrum`` to a DataFrame row per peak in the standard schema."""
    # Accumulators: one entry per peak, filled in order
    annotations: list[str] = []
    ion_types: list[str] = []
    charges: list[int] = []
    flags: list[bool] = []

    # Walk peaks: take the first matched fragment, else mark as unannotated
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

    # Pack: stitch the columns into the project's standard peaks-frame schema
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
    """Match experimental peaks to theoretical fragments of ``peptidoform``; one row per peak.

    Pass real ``precursor_mz`` when ``"precursor"`` is in ``ion_types`` (0.0 only blocks precursor matches)."""
    # Defaults: fill ion_types, build spectrum_utils' letter string, prepare fallback frame
    if ion_types is None:
        ion_types = list(DEFAULTS.ion_types)
    # Translate UI ion-type names into the single-letter alphabet spectrum_utils expects:
    # backbone letters a/b/c/x/y/z pass through; "internal"->m, "precursor"->p, "immonium"->I.
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

    # Fast exit: nothing to annotate or no ion types requested
    if n == 0 or not ion_types_string:
        return _fallback_unannotated_frame()

    # Wrap: build the MsmsSpectrum carrier spectrum_utils operates on
    spec = MsmsSpectrum(
        identifier="experimental",
        precursor_mz=float(precursor_mz),
        precursor_charge=int(charge),
        mz=np.asarray(mz, dtype=float),
        intensity=np.asarray(intensity, dtype=float),
    )

    # Annotate: ProForma fragment matching; any parse/match failure -> unannotated frame
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
    """Drop annotations whose ion type the user didn't ask for and optionally normalise to max=1."""
    # Fast exit: nothing to do on an empty frame
    if annotated.empty:
        return annotated

    # Mask: peaks outside the requested ion_types lose their label but stay in the frame
    keep_letters = set(ion_types) | {"unannotated"}
    unwanted = ~annotated["ion_type"].isin(keep_letters)
    if unwanted.any():
        annotated.loc[unwanted, "annotation"] = "unannotated"
        annotated.loc[unwanted, "ion_type"] = "unannotated"
        annotated.loc[unwanted, "is_annotated"] = False
        annotated.loc[unwanted, "charge"] = 0

    # Normalise: rescale so the tallest peak sits at intensity 1.0 (mirror-plot convention)
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
    """Annotate an experimental spectrum, filter by ion types, optionally normalise — ready for plotting and metrics."""
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
    """Pull the fragment charge out of a Koina label: ``y3+1`` -> 1, ``b5+2`` -> 2, ``y7++`` -> 2."""
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
    """Prepare a Koina-predicted spectrum for plotting; uses the predictor's labels verbatim (no theoretical matching)."""
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
