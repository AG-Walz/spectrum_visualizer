"""Koina REST API client — runs ML models to predict MS/MS fragmentation spectrum (m/z + intensity) 
from peptide sequences; used as fallback reference when no experimental reference mzML exists 
(called by ``core.pipeline``)."""

import logging
from typing import Optional
import pandas as pd
from psm_utils import Peptidoform
from koinapy import Koina

from core.defaults import DEFAULTS, INSTRUMENT_TYPES, KOINA_MODELS

logger = logging.getLogger(__name__)


def convert_peptide_to_unimod(peptidoform_str: str) -> str:
    """Rewrite a ProForma peptidoform in Koina's UniMod-bracket notation 
    (e.g. ``"PEPM[Oxidation]IDEK" -> "PEPM[UNIMOD:35]IDEK"``)."""

    p = Peptidoform(peptidoform_str.split("/")[0])
    fmt = lambda mod: f"[UNIMOD:{mod.definition['id']}]"
    parts = [f"{fmt(m)}-" for m in p.properties.get("n_term") or []]
    for aa, mods in p.parsed_sequence:
        parts.append(aa)
        parts += [fmt(m) for m in (mods or [])]
    parts += [f"-{fmt(m)}" for m in p.properties.get("c_term") or []]
    return "".join(parts)


def predict_koina_spectrum(
    peptidoform_str: str,
    charge: int,
    model_id: str = DEFAULTS.koina_model,
    collision_energy: int = DEFAULTS.koina_collision_energy,
    instrument_type: str = DEFAULTS.koina_instrument_type,
) -> Optional[pd.DataFrame]:
    """Call Koina for a predicted MS/MS spectrum; return ``DataFrame[mz, intensity, annotation]`` or ``None`` on failure."""

    # Convert peptide sequence to UniMod notation
    try:
        unimod_seq = convert_peptide_to_unimod(peptidoform_str)
    except Exception as error:
        logger.warning(f"Koina: failed to convert '{peptidoform_str}' to UniMod: {error}")
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
        logger.warning(f"Koina prediction failed for '{peptidoform_str}': {error}")
        return None

    # Final summary
    if result.empty:
        logger.warning(f"Koina: no non-zero predictions for '{peptidoform_str}' (model={model_id})")
        return None

    logger.info(f"Koina: predicted {len(result)} fragments for '{peptidoform_str}' (model={model_id}, CE={collision_energy}, inst={instrument_type})")
    return result


def get_koina_settings(settings: dict) -> dict:
    """Pull Koina configuration settings from the provided settings, defaulting any missing keys."""
    koina_cfg = settings.get("koina", {})
    return {
        "enabled": koina_cfg.get("enabled", True),
        "model": koina_cfg.get("model", DEFAULTS.koina_model),
        "collision_energy": koina_cfg.get("collision_energy", DEFAULTS.koina_collision_energy),
        "instrument_type": koina_cfg.get("instrument_type", DEFAULTS.koina_instrument_type),
    }
