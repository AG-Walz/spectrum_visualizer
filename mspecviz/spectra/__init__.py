"""Spectrum loading + Koina prediction.

Re-exports the public API of the submodules so callers can keep using
``from mspecviz.spectra import ...`` without caring which file a name lives in.
"""
from mspecviz.spectra.koina_prediction import (
    INSTRUMENT_TYPES,
    KOINA_AVAILABLE,
    KOINA_MODELS,
    convert_peptide_to_unimod,
    get_koina_settings,
    predict_koina_spectrum,
)
from mspecviz.spectra.loaded_spectrum import (
    LoadResult,
    LoadedSpectrum,
    StreamingResult,
    spectrum_lookup_key,
)
from mspecviz.spectra.spectrum_streaming import (
    build_mzml_to_matches_map,
    load_spectra_streaming,
)
from mspecviz.spectra.mzml_workers import _load_spectra_from_mzml
from mspecviz.spectra.mzml_lookup import (
    find_spectrum_index,
    open_indexed_mzml_fast,
)
from mspecviz.spectra.mzml_native_ids import (
    ANY_DIGITS_RE,
    INDEX_TAG,
    NATIVE_ID_NUM_RE,
    NATIVE_ID_TAGS,
    SCAN_TAG,
    SPECTRUM_TAG,
    clean_spectrum_id,
)

__all__ = [
    # loaded_spectrum
    "LoadResult",
    "LoadedSpectrum",
    "StreamingResult",
    "spectrum_lookup_key",
    # spectrum_streaming
    "build_mzml_to_matches_map",
    "load_spectra_streaming",
    # mzml_workers
    "_load_spectra_from_mzml",
    # mzml_native_ids
    "ANY_DIGITS_RE",
    "INDEX_TAG",
    "NATIVE_ID_NUM_RE",
    "NATIVE_ID_TAGS",
    "SCAN_TAG",
    "SPECTRUM_TAG",
    "clean_spectrum_id",
    # mzml_lookup
    "find_spectrum_index",
    "open_indexed_mzml_fast",
    # koina_prediction
    "INSTRUMENT_TYPES",
    "KOINA_AVAILABLE",
    "KOINA_MODELS",
    "convert_peptide_to_unimod",
    "get_koina_settings",
    "predict_koina_spectrum",
]
