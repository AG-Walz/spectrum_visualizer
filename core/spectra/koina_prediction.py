"""Koina REST API client — runs ML models to predict MS/MS fragmentation spectrum (m/z + intensity) 
from peptide sequences; used as fallback reference when no experimental reference mzML exists 
(called by core.pipeline).
"""
import logging
from typing import Optional
import pandas as pd
from psm_utils import Peptidoform, PSM
from koinapy import Koina

from core.defaults import DEFAULTS, INSTRUMENT_TYPES, KOINA_MODELS

logger = logging.getLogger(__name__)


def convert_peptide_to_unimod(peptidoform: Peptidoform) -> str:
    """Rewrite a peptidoform in Koina's UniMod-bracket notation
    (e.g. Peptidoform("PEPM[Oxidation]IDEK") -> "PEPM[UNIMOD:35]IDEK").
    Only the sequence and its modifications are emitted — any precursor charge is ignored.
    """
    fmt = lambda mod: f"[UNIMOD:{mod.definition['id']}]"
    parts = [f"{fmt(m)}-" for m in peptidoform.properties.get("n_term") or []]
    for aa, mods in peptidoform.parsed_sequence:
        parts.append(aa)
        parts += [fmt(m) for m in (mods or [])]
    parts += [f"-{fmt(m)}" for m in peptidoform.properties.get("c_term") or []]
    return "".join(parts)


def predict_koina_spectrum(
    psm: PSM,
    model_id: str = DEFAULTS.koina_model,
    collision_energy: int = DEFAULTS.koina_collision_energy,
    instrument_type: str = DEFAULTS.koina_instrument_type,
) -> Optional[pd.DataFrame]:
    """Call Koina for a predicted MS/MS spectrum of psm's peptidoform;
    return DataFrame[mz, intensity, annotation] or None on failure.
    The peptidoform (modified sequence + charge) is read from psm; model_id,
    collision_energy, and instrument_type are Koina settings, not PSM properties.
    """
    # The peptidoform carries both the sequence and the charge state.
    peptidoform = psm.peptidoform
    charge = psm.get_precursor_charge() or 0
    sequence = peptidoform.modified_sequence

    # Convert the peptidoform sequence to UniMod notation
    try:
        unimod_seq = convert_peptide_to_unimod(peptidoform)
    except Exception as error:
        logger.warning(f"Koina: failed to convert '{sequence}' to UniMod: {error}")
        return None

    # Build one-row input Dataframe
    try:
        koina_input = pd.DataFrame({
            "peptide_sequences": [unimod_seq],
            "precursor_charges": [charge],
            "collision_energies": [collision_energy],
            "instrument_types": [instrument_type],
        })
        # Predict spectrum with Koina
        predictions = Koina(model_id).predict(koina_input)
        # Clean up annotation strings and convert to numeric types
        result = pd.DataFrame({
            "mz": predictions["mz"].astype(float),
            "intensity": predictions["intensities"].astype(float),
            "annotation": predictions["annotation"].astype(str).str.replace(r"^b'|'$", "", regex=True),
        })
        # Drop zero-intensity predictions (especially for low CE)
        result = result[result["intensity"] > 0].reset_index(drop=True)
    # Returns None on any failure
    except Exception as error:
        logger.warning(f"Koina prediction failed for '{sequence}': {error}")
        return None

    # Final summary
    if result.empty:
        logger.warning(f"Koina: no non-zero predictions for '{sequence}' (model={model_id})")
        return None

    logger.info(f"Koina: predicted {len(result)} fragments for '{sequence}' (model={model_id}, CE={collision_energy}, inst={instrument_type})")
    return result


def get_koina_settings(settings: dict) -> dict:
    """Pull Koina configuration settings from the provided settings, defaulting any missing keys.
    """
    koina_cfg = settings.get("koina", {})
    return {
        "enabled": koina_cfg.get("enabled", True),
        "model": koina_cfg.get("model", DEFAULTS.koina_model),
        "collision_energy": koina_cfg.get("collision_energy", DEFAULTS.koina_collision_energy),
        "instrument_type": koina_cfg.get("instrument_type", DEFAULTS.koina_instrument_type),
    }
