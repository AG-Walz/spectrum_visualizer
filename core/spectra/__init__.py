"""Spectrum loading + Koina prediction.

Re-exports the package's public API so callers can use
``from core.spectra import ...`` without caring which file a name lives in.
"""
from core.spectra.koina_prediction import (
    KOINA_AVAILABLE,
    convert_peptide_to_unimod,
    get_koina_settings,
    predict_koina_spectrum,
)
from core.spectra.loaded_spectrum import LoadedSpectrum, spectrum_lookup_key
from core.spectra.spectrum_streaming import load_spectra_streaming

__all__ = [
    "KOINA_AVAILABLE",
    "LoadedSpectrum",
    "convert_peptide_to_unimod",
    "get_koina_settings",
    "load_spectra_streaming",
    "predict_koina_spectrum",
    "spectrum_lookup_key",
]
