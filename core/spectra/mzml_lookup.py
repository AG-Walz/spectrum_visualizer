"""mzML spectrum-ID resolution layer — maps native spectrum IDs to array positions for one mzML file.
Builds an ID -> array-position index for a single mzML so callers can fetch a
spectrum by its native ID (scan=N / index=N / spectrum=N) (called by core.spectra.mzml_workers)
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pyopenms as oms

from core.spectra.mzml_native_ids import (
    INDEX_TAG,
    SCAN_TAG,
    SPECTRUM_TAG,
    clean_spectrum_id,
    extract_id_number,
)

logger = logging.getLogger(__name__)


def _detect_id_format(indexed_exp: oms.OnDiscMSExperiment, n_spectra: int) -> Tuple[str, int]:
    """Probe the first spectrum's native ID and classify the file's ID scheme.
    Returns (format, offset) where format is one of "index", "scan", "spectrum" or "unknown". 
    offset is meaningful only for the index format: it is the index=N value of array position 0,
    so 0 for normal files and 1 for the Bruker timsTOF 1-based gotcha. For all
    other formats offset is 0 and unused.
    """
    # Empty file — nothing to probe.
    if n_spectra == 0:
        return "unknown", 0
    # Inspect spectrum 0's native ID; the first matching tag wins.
    first_id = indexed_exp.getSpectrum(0).getNativeID()
    # index=N — the value at array position 0 gives the base offset (0 or 1).
    if INDEX_TAG in first_id:
        first_number = extract_id_number(first_id)
        return "index", (first_number if first_number is not None else 0)
    if SCAN_TAG in first_id:
        return "scan", 0
    if SPECTRUM_TAG in first_id:
        return "spectrum", 0
    # No recognised tag — caller must fall back to a full scan.
    return "unknown", 0


def _try_candidate_indices(
    indexed_exp: oms.OnDiscMSExperiment, clean_id: str, candidates: List[int],
    id_to_index: Dict[str, int],
) -> bool:
    """Probe a few likely array positions for a spectrum whose native ID matches.
    On the first match, records both the cleaned ID and a _num_<n> numeric
    alias in id_to_index and returns True; returns False if no
    candidate matches. Mutates id_to_index in place.
    """
    for candidate_index in candidates:
        # Compare canonicalised IDs so tag/format differences do not block a match.
        check_id = indexed_exp.getSpectrum(candidate_index).getNativeID()
        if clean_id != clean_spectrum_id(check_id):
            continue
        # Hit — record the ID and a numeric alias so later lookups skip the probe.
        id_to_index[clean_id] = candidate_index
        check_number = extract_id_number(check_id)
        if check_number is not None:
            id_to_index[f"_num_{check_number}"] = candidate_index
        return True
    return False


def _binary_search_for_number(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, target_number: int,
    clean_id: str, id_to_index: Dict[str, int],
) -> bool:
    """Binary-search the spectra for one matching scan number.
    Assumes native-ID numbers increase monotonically with array position (the
    usual mzML layout). On a match, records the cleaned ID and a _num_<n>
    alias in id_to_index and returns True; returns False if the
    number is absent or any spectrum read fails. Mutates id_to_index.
    """
    low, high = 0, n_spectra - 1
    while low <= high:
        mid = (low + high) // 2
        # A failed read collapses the assumption of an ordered file — bail out.
        try:
            mid_id = indexed_exp.getSpectrum(mid).getNativeID()
        except Exception:
            return False
        # An unparseable number means the monotonic assumption no longer holds.
        mid_number = extract_id_number(mid_id)
        if mid_number is None:
            return False
        # Match — record the position and its numeric alias.
        if mid_number == target_number:
            id_to_index[clean_id] = mid
            id_to_index[f"_num_{mid_number}"] = mid
            return True
        # Narrow the search window towards the target number.
        if mid_number < target_number:
            low = mid + 1
        else:
            high = mid - 1
    return False


def _resolve_targeted_id(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, spectrum_id: str,
    id_to_index: Dict[str, int],
) -> bool:
    """Locate one needed spectrum ID and record its array position in id_to_index.
    Tries the cheap path first — the scan number itself and its off-by-one
    neighbour as direct array positions — then falls back to a binary search by
    scan number. Returns True if the ID was resolved and recorded, False
    if its number cannot be parsed or no spectrum matches.
    """
    clean_id = clean_spectrum_id(spectrum_id)
    extracted_number = extract_id_number(clean_id)
    # No parseable number — nothing to search on.
    if extracted_number is None:
        return False
    # Off-by-one guess: the number may be the array position directly, or 1-based.
    candidates = [c for c in (extracted_number, extracted_number - 1) if 0 <= c < n_spectra]
    if _try_candidate_indices(indexed_exp, clean_id, candidates, id_to_index):
        return True
    # Guess missed — search the whole file assuming numbers are monotonic.
    return _binary_search_for_number(
        indexed_exp, n_spectra, extracted_number, clean_id, id_to_index,
    )


def _build_full_scan_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, id_to_index: Dict[str, int],
) -> None:
    """Walk every spectrum and index it by native ID — the last-resort full scan.
    Used when the ID format is unrecognised or targeted lookups all failed.
    Records each spectrum's raw native ID plus a _num_<n> numeric alias in
    id_to_index; spectra that fail to read are skipped. Mutates
    id_to_index in place and returns nothing.
    """
    for index in range(n_spectra):
        # Skip unreadable spectra rather than aborting the whole scan.
        try:
            native_id = indexed_exp.getSpectrum(index).getNativeID()
        except Exception:
            continue
        # Index by the raw native ID and by its extracted number.
        id_to_index[native_id] = index
        number = extract_id_number(native_id)
        if number is not None:
            id_to_index[f"_num_{number}"] = index


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """Open one mzML and build its native-ID -> array-position map. Returns (experiment, id_to_index). 
    Three strategies, cheapest first: (1) for the index format the map needs no per-spectrum entries 
    since find_spectrum_index derives positions from the offset directly; 
    (2) when needed_ids is supplied for a scan / spectrum file, only those IDs are resolved; 
    (3) otherwise every spectrum is walked.
    """
    # Open the file on disc and probe its native-ID format and index offset.
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))
    n_spectra = indexed_exp.getNrSpectra()

    # Seed the map with the file-describing sentinels read by find_spectrum_index.
    id_format, index_offset = _detect_id_format(indexed_exp, n_spectra)
    id_to_index: Dict[str, int] = {
        "__format__": id_format,
        "__n_spectra__": n_spectra,
        "__index_offset__": index_offset,
    }

    # Fast path: index=N maps straight to an array position via the offset, so
    # no per-spectrum entries are needed at all.
    if id_format == "index":
        logger.debug(f"Fast path: index format detected (offset={index_offset}), "
                     f"skipping full scan ({n_spectra} spectra)")
        return indexed_exp, id_to_index

    # Targeted lookup: resolve only the scan/spectrum IDs the caller actually needs.
    if needed_ids and id_format in ("scan", "spectrum"):
        found = sum(
            _resolve_targeted_id(indexed_exp, n_spectra, sid, id_to_index)
            for sid in needed_ids
        )
        logger.debug(
            f"Targeted lookup: found {found}/{len(needed_ids)} spectra in "
            f"{Path(mzml_path).name}"
        )
        # Accept the targeted result only if it found something; otherwise the
        # monotonic assumption likely failed and the full scan still runs.
        if found > 0:
            return indexed_exp, id_to_index

    # Fallback: format unknown or targeted lookup found nothing — index everything.
    logger.debug(f"Full scan fallback for {Path(mzml_path).name} ({n_spectra} spectra)")
    _build_full_scan_index(indexed_exp, n_spectra, id_to_index)
    return indexed_exp, id_to_index


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int,
) -> Optional[int]:
    """Resolve one spectrum ID to its array position using a prebuilt map.
    id_to_index is the map returned by open_indexed_mzml_fast and
    n_spectra its spectrum count. Returns the integer array position, or
    None if the ID is empty, has no parseable number, or cannot be located.
    Resolution order: exact cleaned-ID hit, then _num_<n> numeric alias,
    then — for the index format — the position computed from the native
    number minus the auto-detected index offset.
    """
    clean_id = clean_spectrum_id(spectrum_id)
    # Empty ID — nothing to resolve.
    if not clean_id:
        return None
    # Exact match on the cleaned native ID.
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    # No exact entry — fall back to matching on the extracted number.
    extracted = extract_id_number(clean_id)
    if extracted is None:
        return None

    # Numeric alias recorded by the targeted / full-scan index builders.
    sentinel = id_to_index.get(f"_num_{extracted}")
    if sentinel is not None:
        return sentinel

    # index format: derive the position arithmetically, correcting for the
    # 0- vs 1-based offset (Bruker timsTOF emits 1-based index IDs).
    if id_to_index.get("__format__") == "index":
        offset = id_to_index.get("__index_offset__", 0)
        position = extracted - offset
        if 0 <= position < n_spectra:
            return position

    return None
