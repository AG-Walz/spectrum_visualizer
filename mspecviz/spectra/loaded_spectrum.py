"""Dataclasses returned by the mzML spectrum loader.

Kept separate from the loader workers so they can be imported cheaply (no
pyopenms pull-in) by consumers that only need the types.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # avoid circular import at runtime
    from psm_utils import PSM


@dataclass
class LoadedSpectrum:
    """A loaded spectrum; fields duplicating ``match.*`` exist so Koina refs (``match=None``) carry the same shape."""
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
    """Filename stem for plot titles / picker tables; ``"Koina"`` for synthetic refs.

    Defined as a module-level function (not a method on LoadedSpectrum) so that
    older instances cached in Streamlit's session_state across an auto-reload
    keep working — a method/property would only resolve on instances of the
    *currently-loaded* class.
    """
    return Path(spectrum.source_file).stem if spectrum.source_file else "Koina"


@dataclass
class StreamingResult:
    """One query spectrum with its variant key and matching reference spectra."""
    query_spectrum: LoadedSpectrum
    reference_spectra: List[LoadedSpectrum]
    variant_key: str


@dataclass
class LoadResult:
    """Per-file output from ``_load_spectra_from_mzml``.

    Replaces an anonymous 6-tuple. *loaded* is a list of
    ``(variant_key, LoadedSpectrum)`` pairs so the caller can bucket spectra
    by variant without re-scanning.
    """
    display_name: str
    is_reference: bool
    loaded: List[Tuple[str, LoadedSpectrum]]
    loaded_count: int
    not_found_count: int
    errors: List[str]


def spectrum_lookup_key(spec: "LoadedSpectrum") -> str:
    """Return ``spectrum_id|peptidoform`` — stable PSM-unique key for lookup.

    The same physical spectrum may appear with different peptidoforms across
    PSMs, and spectrum IDs alone can collide across mzML files. Pair the key
    with ``build_psm_table()``'s ``psm_lookup`` dict to resolve ``psm_id`` when
    building peaks and sequence DataFrames. Module-level (not a method) for
    the same Streamlit auto-reload reason as :func:`source_label`.
    """
    return f"{spec.spectrum_id}|{spec.peptidoform}"
