"""
Koina AI spectrum prediction for MSpecViz.

Generates predicted MS/MS fragment spectra using Koina deep-learning models
(via koinapy). Unlike PyOpenMS' theoretical generator which only produces m/z
positions, Koina models predict both m/z AND relative intensities, making them
useful for quantitative comparison against experimental spectra.

Used as a fallback reference when no experimental reference mzML is available.
"""
import logging
from typing import Optional

import pandas as pd
from psm_utils import Peptidoform

# KOINA_MODELS / INSTRUMENT_TYPES sourced from defaults (single catalogue).
from mspecviz.common.defaults import DEFAULTS, INSTRUMENT_TYPES, KOINA_MODELS

logger = logging.getLogger(__name__)

try:
    from koinapy import Koina
    KOINA_AVAILABLE = True
except ImportError:
    KOINA_AVAILABLE = False


def _warn_none(message: str) -> None:
    """Log *message* at WARNING and return ``None`` — used by the predict-spectrum branches."""
    logger.warning(message)
    return None


def convert_peptide_to_unimod(peptidoform_str: str) -> str:
    """Convert a ProForma peptidoform to Koina's UniMod-bracket notation.

    ``"PEPM[Oxidation]IDEK"`` → ``"PEPM[UNIMOD:35]IDEK"``. Reads the resolved
    UniMod accession from each modification's ``definition`` (psm_utils delegates
    the resolution to pyteomics' UnimodResolver). A modification without a
    UniMod entry — e.g. a bare mass-shift like ``[+15.99]`` or PSI-MOD/XLMOD —
    raises ``KeyError`` here and is caught upstream by ``predict_koina_spectrum``.
    """
    p = Peptidoform(peptidoform_str.split("/")[0])
    fmt = lambda mod: f"[UNIMOD:{mod.definition['id']}]"
    parts = [f"{fmt(m)}-" for m in p.properties.get("n_term") or []]
    for aa, mods in p.parsed_sequence:
        parts.append(aa)
        parts += [fmt(m) for m in (mods or [])]
    parts += [f"-{fmt(m)}" for m in p.properties.get("c_term") or []]
    return "".join(parts)


def _build_koina_input(unimod_seq: str, charge: int, collision_energy: int,
                       instrument_type: str) -> pd.DataFrame:
    """Build the single-row input DataFrame Koina expects."""
    return pd.DataFrame({
        "peptide_sequences": [unimod_seq],
        "precursor_charges": [charge],
        "collision_energies": [collision_energy],
        "instrument_types": [instrument_type],
    })


def _koina_predictions_to_frame(predictions) -> pd.DataFrame:
    """Convert raw Koina output into our standard ``[mz, intensity, annotation]`` frame."""
    result = pd.DataFrame({
        "mz": predictions["mz"].astype(float),
        "intensity": predictions["intensities"].astype(float),
        "annotation": predictions["annotation"].astype(str).str.replace(r"^b'|'$", "", regex=True),
    })
    return result[result["intensity"] > 0].reset_index(drop=True)


def predict_koina_spectrum(
    peptidoform_str: str,
    charge: int,
    model_id: str = DEFAULTS.koina_model,
    collision_energy: int = DEFAULTS.koina_collision_energy,
    instrument_type: str = DEFAULTS.koina_instrument_type,
) -> Optional[pd.DataFrame]:
    """Predict MS/MS spectrum via Koina. Returns DataFrame [mz, intensity, annotation] or None."""
    if not KOINA_AVAILABLE:
        return _warn_none("koinapy not installed -- cannot generate Koina predictions")

    try:
        unimod_seq = convert_peptide_to_unimod(peptidoform_str)
    except Exception as error:
        return _warn_none(f"Koina: failed to convert '{peptidoform_str}' to UniMod: {error}")

    try:
        predictions = Koina(model_id).predict(
            _build_koina_input(unimod_seq, charge, collision_energy, instrument_type)
        )
        result = _koina_predictions_to_frame(predictions)
    except Exception as error:
        return _warn_none(f"Koina prediction failed for '{peptidoform_str}': {error}")

    if result.empty:
        return _warn_none(f"Koina: no non-zero predictions for '{peptidoform_str}' (model={model_id})")

    logger.info(f"Koina: predicted {len(result)} fragments for '{peptidoform_str}' (model={model_id}, CE={collision_energy}, inst={instrument_type})")
    return result


def get_koina_settings(settings: dict) -> dict:
    """Extract Koina config from settings, falling back to defaults."""
    koina_cfg = settings.get("koina", {})
    return {
        "enabled": koina_cfg.get("enabled", True),
        "model": koina_cfg.get("model", DEFAULTS.koina_model),
        "collision_energy": koina_cfg.get("collision_energy", DEFAULTS.koina_collision_energy),
        "instrument_type": koina_cfg.get("instrument_type", DEFAULTS.koina_instrument_type),
    }
