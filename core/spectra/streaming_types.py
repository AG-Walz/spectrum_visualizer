"""Streaming carrier types + PSM helpers for the mzML loader (split out so
callers can import the types without pulling in pyopenms). The carrier the
pipeline streams is the psm_utils PSM itself, with its peaks attached to
psm.spectrum; these dataclasses just bundle it with its variant key.
"""
from dataclasses import dataclass
from typing import List, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from psm_utils import PSM


def source_label(psm) -> str:
    """Return the run name for plot titles / picker tables, read off the PSM.
    Koina references carry run="Koina:<model>"; falls back to "Koina" when unset.
    """
    return psm.run or "Koina"


@dataclass
class StreamingResult:
    """Single query spectrum (a carrier PSM) plus its peptide variant key and the
    matching reference PSMs (yielded by load_spectra_streaming and consumed by
    core.pipeline.query_processing.process_query_result).
    """
    query_spectrum: "PSM"
    reference_spectra: List["PSM"]
    variant_key: str


@dataclass
class LoadResult:
    """Per-file output from _load_spectra_from_mzml;
    loaded is [(variant_key, PSM), ...] so the orchestrator can extract each
    spectrum (carrier PSM) by peptide variant without rescanning.
    """
    display_name: str
    is_reference: bool
    loaded: List[Tuple[str, "PSM"]]
    loaded_count: int
    not_found_count: int
    errors: List[str]


def psm_lookup_key(psm: "PSM") -> str:
    """Return "spectrum_id|peptidoform" — the stable per-PSM key for joining the
    peaks/sequence/PSM DataFrames. Pairs both fields because spectrum IDs alone
    collide across mzML files and across peptidoforms.
    """
    return f"{psm.spectrum_id}|{psm.peptidoform.modified_sequence}"
