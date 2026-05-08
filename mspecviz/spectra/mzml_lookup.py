"""mzML spectrum-ID indexing with format-aware fast paths.

CRITICAL: mzML ``index=N`` native IDs are not always 0-based — timsTOF/Bruker
emit 1-based. ``open_indexed_mzml_fast`` auto-detects the offset from the first
spectrum and ``find_spectrum_index`` applies it; removing this offset logic
silently corrupts lookups for Bruker files.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyopenms as oms

from mspecviz.spectra.mzml_native_ids import (
    BARE_DIGITS_PATTERN,
    INDEX_TAG,
    NATIVE_ID_NUMBER_PATTERN,
    SCAN_TAG,
    SPECTRUM_TAG,
    clean_spectrum_id,
    extract_id_number,
)

logger = logging.getLogger(__name__)


def _detect_id_format(indexed_exp: oms.OnDiscMSExperiment, n_spectra: int) -> Tuple[str, int]:
    """Inspect the first spectrum's native ID; return ``(format, offset)``.

    ``offset`` is meaningful only for ``index=N`` (0 for 0-based, 1 for 1-based).
    """
    if n_spectra == 0:
        return "unknown", 0
    first_id = indexed_exp.getSpectrum(0).getNativeID()
    if INDEX_TAG in first_id:
        first_number = extract_id_number(first_id)
        return "index", (first_number if first_number is not None else 0)
    if SCAN_TAG in first_id:
        return "scan", 0
    if SPECTRUM_TAG in first_id:
        return "spectrum", 0
    return "unknown", 0


def _try_candidate_indices(
    indexed_exp: oms.OnDiscMSExperiment, clean_id: str, candidates: List[int],
    id_to_index: Dict[str, int],
) -> bool:
    """Try a small list of candidate array positions; persist the first match."""
    for candidate_index in candidates:
        try:
            check_id = indexed_exp.getSpectrum(candidate_index).getNativeID()
        except Exception:
            continue
        if clean_id != clean_spectrum_id(check_id):
            continue
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
    """Scan IDs are typically monotonic — binary-search by extracted number."""
    low, high = 0, n_spectra - 1
    while low <= high:
        mid = (low + high) // 2
        try:
            mid_id = indexed_exp.getSpectrum(mid).getNativeID()
        except Exception:
            return False
        mid_number = extract_id_number(mid_id)
        if mid_number is None:
            return False
        if mid_number == target_number:
            id_to_index[clean_id] = mid
            id_to_index[f"_num_{mid_number}"] = mid
            return True
        if mid_number < target_number:
            low = mid + 1
        else:
            high = mid - 1
    return False


def _resolve_targeted_id(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, spectrum_id: str,
    id_to_index: Dict[str, int],
) -> bool:
    """Resolve one needed spectrum ID via off-by-one probe then binary search."""
    clean_id = clean_spectrum_id(spectrum_id)
    extracted_number = extract_id_number(clean_id)
    if extracted_number is None:
        return False
    candidates = []
    if 0 <= extracted_number < n_spectra:
        candidates.append(extracted_number)
    if extracted_number > 0 and (extracted_number - 1) < n_spectra:
        candidates.append(extracted_number - 1)
    if _try_candidate_indices(indexed_exp, clean_id, candidates, id_to_index):
        return True
    return _binary_search_for_number(
        indexed_exp, n_spectra, extracted_number, clean_id, id_to_index,
    )


def _build_targeted_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, needed_ids: List[str],
    id_to_index: Dict[str, int],
) -> int:
    """Populate *id_to_index* for each needed ID; returns count found."""
    found = 0
    for spectrum_id in needed_ids:
        if _resolve_targeted_id(indexed_exp, n_spectra, spectrum_id, id_to_index):
            found += 1
    return found


def _build_full_scan_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, id_to_index: Dict[str, int],
) -> None:
    """Fall back to walking every spectrum and recording its native ID."""
    for index in range(n_spectra):
        try:
            native_id = indexed_exp.getSpectrum(index).getNativeID()
        except Exception:
            continue
        id_to_index[native_id] = index
        number_match = BARE_DIGITS_PATTERN.search(native_id)
        if number_match:
            id_to_index[f"_num_{number_match.group(1)}"] = index


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """Open mzML and build a targeted ID→position map (skip full scan when possible).

    ``index=N``: O(1) via detected offset. ``scan=N`` with ``needed_ids``: targeted
    lookup. Otherwise: full scan. The returned dict carries three sentinel keys
    (``__format__``, ``__n_spectra__``, ``__index_offset__``) used by
    :func:`find_spectrum_index`.
    """
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))
    n_spectra = indexed_exp.getNrSpectra()

    id_format, index_offset = _detect_id_format(indexed_exp, n_spectra)
    id_to_index: Dict[str, int] = {
        "__format__": id_format,
        "__n_spectra__": n_spectra,
        "__index_offset__": index_offset,
    }

    if id_format == "index":
        logger.debug(f"Fast path: index format detected (offset={index_offset}), "
                     f"skipping full scan ({n_spectra} spectra)")
        return indexed_exp, id_to_index

    if needed_ids and id_format in ("scan", "spectrum"):
        found = _build_targeted_index(indexed_exp, n_spectra, needed_ids, id_to_index)
        logger.debug(
            f"Targeted lookup: found {found}/{len(needed_ids)} spectra in "
            f"{Path(mzml_path).name}"
        )
        if found > 0:
            return indexed_exp, id_to_index

    logger.debug(f"Full scan fallback for {Path(mzml_path).name} ({n_spectra} spectra)")
    _build_full_scan_index(indexed_exp, n_spectra, id_to_index)
    return indexed_exp, id_to_index


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int,
) -> Optional[int]:
    """Resolve a spectrum ID to its array position via the map from ``open_indexed_mzml_fast``.

    Tries direct lookup, then ``(scan|index|spectrum)=N`` with offset adjustment,
    then a bare-digit fallback.
    """
    clean_id = clean_spectrum_id(spectrum_id)
    if not clean_id:
        return None

    id_format = id_to_index.get("__format__", "unknown")
    index_offset = id_to_index.get("__index_offset__", 0)

    # Strategy 1: direct lookup
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    # Strategy 2: structured ID with (scan|index|spectrum)=N
    scan_match = NATIVE_ID_NUMBER_PATTERN.search(clean_id)
    if scan_match:
        extracted_number = int(scan_match.group(1))

        if id_format == "index":
            array_position = extracted_number - index_offset
            if 0 <= array_position < n_spectra:
                return array_position

        number_key = f"_num_{extracted_number}"
        if number_key in id_to_index:
            return id_to_index[number_key]

    # Strategy 3: bare number fallback
    number_match = BARE_DIGITS_PATTERN.search(clean_id)
    if number_match:
        extracted_number = int(number_match.group(1))
        number_key = f"_num_{extracted_number}"
        if number_key in id_to_index:
            return id_to_index[number_key]

        if id_format == "index":
            array_position = extracted_number - index_offset
            if 0 <= array_position < n_spectra:
                return array_position

    return None
