"""Shared dataclasses for the mzML loader (split out so callers can import the types without pulling in pyopenms).
"""
from dataclasses import dataclass
from typing import List, Optional, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from psm_utils import PSM


@dataclass
class LoadedSpectrum:
    """Single loaded MS2 spectrum object (peak arrays + metadata)
    created in core.spectra.spectrum_streaming._build_loaded_spectrum wrapping
    pyopenms.MSSpectrum and its PSM together, and in
    core.pipeline.query_processing.get_koina_prediction building synthetic
    LoadedSpectra without PSM.
    """
    spectrum_id: str
    mz: np.ndarray
    intensity: np.ndarray
    precursor_mz: float
    precursor_charge: int
    peptidoform: str
    base_sequence: str
    source_file: str
    is_reference: bool = False
    rt: float = 0.0
    match: Optional["PSM"] = None # koina spectra have no PSM match


def source_label(psm) -> str:
    """Return the run name for plot titles / picker tables, read off the PSM.
    Koina references carry run="Koina:<model>"; falls back to "Koina" when unset.
    """
    return psm.run or "Koina"


@dataclass
class StreamingResult:
    """Single query (experimental) spectrum plus its peptide variant key and
    all the matching reference spectra (yielded by load_spectra_streaming
    and consumed by core.pipeline.query_processing.process_query_result).
    """
    query_spectrum: LoadedSpectrum
    reference_spectra: List[LoadedSpectrum]
    variant_key: str


@dataclass
class LoadResult:
    """Per-file output from _load_spectra_from_mzml; 
    loaded is [(variant_key, LoadedSpectrum), ...] so the orchestrator 
    can extract each spectrum by peptide variant without rescanning.
    """
    display_name: str
    is_reference: bool
    loaded: List[Tuple[str, LoadedSpectrum]]
    loaded_count: int
    not_found_count: int
    errors: List[str]


def spectrum_lookup_key(spec: "LoadedSpectrum") -> str:
    """Return "spectrum_id|peptidoform"  — stable per-PSM key for joining peaks/sequence DataFrames.
    Pairs both fields because spectrum IDs alone collide across mzML files and across peptidoforms.
    """
    return f"{spec.spectrum_id}|{spec.peptidoform}"


def psm_lookup_key(psm: "PSM") -> str:
    """PSM-side counterpart to spectrum_lookup_key — same "spectrum_id|peptidoform"
    format so PSMs and LoadedSpectra hash to identical keys when joining tables.
    """
    return f"{psm.spectrum_id}|{psm.peptidoform.modified_sequence}"
