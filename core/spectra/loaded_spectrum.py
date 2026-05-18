"""Shared dataclasses for the mzML loader (split out so callers can import the types without pulling in pyopenms)."""
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from psm_utils import PSM


@dataclass
class LoadedSpectrum:
    """One loaded spectrum (peaks + metadata); fields duplicate ``match.*`` so Koina refs (``match=None``) keep the same shape."""
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
    match: Optional["PSM"] = None


def source_label(spectrum) -> str:
    """Filename stem for plot titles and picker tables; ``"Koina"`` for synthetic refs.

    Module-level (not a method) so cached ``LoadedSpectrum`` instances survive Streamlit reloads."""
    return Path(spectrum.source_file).stem if spectrum.source_file else "Koina"


@dataclass
class StreamingResult:
    """One query spectrum plus its variant key and the matching reference spectra (yielded by ``load_spectra_streaming``)."""
    query_spectrum: LoadedSpectrum
    reference_spectra: List[LoadedSpectrum]
    variant_key: str


@dataclass
class LoadResult:
    """Per-file output from ``_load_spectra_from_mzml``; ``loaded`` is ``[(variant_key, LoadedSpectrum), ...]`` so the orchestrator can bucket by variant without rescanning."""
    display_name: str
    is_reference: bool
    loaded: List[Tuple[str, LoadedSpectrum]]
    loaded_count: int
    not_found_count: int
    errors: List[str]


def spectrum_lookup_key(spec: "LoadedSpectrum") -> str:
    """Return ``"spectrum_id|peptidoform"`` — stable per-PSM key for joining peaks/sequence DataFrames.

    Pairs both fields because spectrum IDs alone collide across mzML files and across peptidoforms."""
    return f"{spec.spectrum_id}|{spec.peptidoform}"
