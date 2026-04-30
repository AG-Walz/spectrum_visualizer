"""
Koina AI spectrum prediction for MSpecViz.

Generates predicted MS/MS fragment spectra using Koina deep-learning models
(via koinapy). Unlike PyOpenMS' theoretical generator which only produces m/z
positions, Koina models predict both m/z AND relative intensities, making them
useful for quantitative comparison against experimental spectra.

Used as a fallback reference when no experimental reference mzML is available.
"""
import logging
from typing import Dict, Optional

import pandas as pd

from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.fragment_utils import peptidoform_to_aasequence

logger = logging.getLogger(__name__)

try:
    from koinapy import Koina
    KOINA_AVAILABLE = True
except ImportError:
    KOINA_AVAILABLE = False

KOINA_MODELS: Dict[str, str] = {
    # MS2PIP family
    "MS2PIP TimsTOF 2024": "ms2pip_timsTOF2024",
    "MS2PIP HCD 2021": "ms2pip_HCD2021",
    "MS2PIP Immuno HCD": "ms2pip_Immuno_HCD",
    "MS2PIP TTOF5600": "ms2pip_TTOF5600",
    "MS2PIP CID TMT": "ms2pip_CID_TMT",
    "MS2PIP iTRAQ Phospho": "ms2pip_iTRAQphospho",
    # Prosit family
    "Prosit 2019": "Prosit_2019_intensity",
    "Prosit 2020 HCD": "Prosit_2020_intensity_HCD",
    "Prosit 2020 CID": "Prosit_2020_intensity_CID",
    "Prosit 2020 TMT": "Prosit_2020_intensity_TMT",
    "Prosit 2023 timsTOF": "Prosit_2023_intensity_timsTOF",
    # AlphaPept / PeptDeep
    "AlphaPept Generic": "AlphaPept_ms2_generic",
    # Other
    "UniSpec": "UniSpec",
    "PredFull": "PredFull",
}

INSTRUMENT_TYPES = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]


def convert_peptide_to_unimod(peptidoform_str: str) -> str:
    """Convert a ProForma peptidoform to UniMod bracket notation for Koina input."""
    aa_seq = peptidoform_to_aasequence(peptidoform_str)
    # toUniModString() uses parentheses: PEPTM(UniMod:35)IDE
    # Koina expects brackets: PEPTM[UNIMOD:35]IDE
    unimod_str = str(aa_seq.toUniModString())
    return unimod_str.replace("(", "[").replace(")", "]")


def predict_koina_spectrum(
    peptidoform_str: str,
    charge: int,
    model_id: str = DEFAULTS.koina_model,
    collision_energy: int = DEFAULTS.koina_collision_energy,
    instrument_type: str = DEFAULTS.koina_instrument_type,
) -> Optional[pd.DataFrame]:
    """Predict MS/MS spectrum via Koina. Returns DataFrame [mz, intensity, annotation] or None."""
    if not KOINA_AVAILABLE:
        logger.warning("koinapy not installed -- cannot generate Koina predictions")
        return None

    try:
        unimod_seq = convert_peptide_to_unimod(peptidoform_str)
    except Exception as e:
        logger.warning(f"Koina: failed to convert '{peptidoform_str}' to UniMod: {e}")
        return None

    try:
        koina = Koina(model_id)

        inputs = pd.DataFrame({
            "peptide_sequences": [unimod_seq],
            "precursor_charges": [charge],
            "collision_energies": [collision_energy],
            "instrument_types": [instrument_type],
        })

        predictions = koina.predict(inputs)

        # Clean byte-string artifacts from annotation column (b'y3+1' -> y3+1)
        annotations = (
            predictions["annotation"]
            .astype(str)
            .str.replace(r"^b'|'$", "", regex=True)
        )

        result = pd.DataFrame({
            "mz": predictions["mz"].astype(float),
            "intensity": predictions["intensities"].astype(float),
            "annotation": annotations,
        })

        result = result[result["intensity"] > 0].reset_index(drop=True)

        if result.empty:
            logger.warning(f"Koina: no non-zero predictions for '{peptidoform_str}' (model={model_id})")
            return None

        logger.info(
            f"Koina: predicted {len(result)} fragments for '{peptidoform_str}' "
            f"(model={model_id}, CE={collision_energy}, inst={instrument_type})"
        )
        return result

    except Exception as e:
        logger.warning(f"Koina prediction failed for '{peptidoform_str}': {e}")
        return None


def get_koina_settings(settings: dict) -> dict:
    """Extract Koina config from settings, falling back to defaults."""
    koina_cfg = settings.get("koina", {})
    return {
        "enabled": koina_cfg.get("enabled", True),
        "model": koina_cfg.get("model", DEFAULTS.koina_model),
        "collision_energy": koina_cfg.get("collision_energy", DEFAULTS.koina_collision_energy),
        "instrument_type": koina_cfg.get("instrument_type", DEFAULTS.koina_instrument_type),
    }
