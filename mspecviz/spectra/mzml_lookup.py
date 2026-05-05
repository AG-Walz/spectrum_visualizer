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
    ANY_DIGITS_RE,
    INDEX_TAG,
    NATIVE_ID_NUM_RE,
    SCAN_TAG,
    SPECTRUM_TAG,
    _extract_id_number,
    clean_spectrum_id,
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
        first_num = _extract_id_number(first_id)
        return "index", (first_num if first_num is not None else 0)
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
    for candidate_idx in candidates:
        try:
            check_id = indexed_exp.getSpectrum(candidate_idx).getNativeID()
        except Exception:
            continue
        if clean_id != clean_spectrum_id(check_id):
            continue
        id_to_index[clean_id] = candidate_idx
        check_num = _extract_id_number(check_id)
        if check_num is not None:
            id_to_index[f"_num_{check_num}"] = candidate_idx
        return True
    return False


def _binary_search_for_num(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, num: int,
    clean_id: str, id_to_index: Dict[str, int],
) -> bool:
    """Scan IDs are typically monotonic — binary-search by extracted number."""
    lo, hi = 0, n_spectra - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        try:
            mid_id = indexed_exp.getSpectrum(mid).getNativeID()
        except Exception:
            return False
        mid_num = _extract_id_number(mid_id)
        if mid_num is None:
            return False
        if mid_num == num:
            id_to_index[clean_id] = mid
            id_to_index[f"_num_{mid_num}"] = mid
            return True
        if mid_num < num:
            lo = mid + 1
        else:
            hi = mid - 1
    return False


def _resolve_targeted_id(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, sid: str,
    id_to_index: Dict[str, int],
) -> bool:
    """Resolve one needed spectrum ID via off-by-one probe then binary search."""
    clean_id = clean_spectrum_id(sid)
    num = _extract_id_number(clean_id)
    if num is None:
        return False
    candidates = []
    if 0 <= num < n_spectra:
        candidates.append(num)
    if num > 0 and (num - 1) < n_spectra:
        candidates.append(num - 1)
    if _try_candidate_indices(indexed_exp, clean_id, candidates, id_to_index):
        return True
    return _binary_search_for_num(indexed_exp, n_spectra, num, clean_id, id_to_index)


def _build_targeted_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, needed_ids: List[str],
    id_to_index: Dict[str, int],
) -> int:
    """Populate *id_to_index* for each needed ID; returns count found."""
    found = 0
    for sid in needed_ids:
        if _resolve_targeted_id(indexed_exp, n_spectra, sid, id_to_index):
            found += 1
    return found


def _build_full_scan_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, id_to_index: Dict[str, int],
) -> None:
    """Fall back to walking every spectrum and recording its native ID."""
    for i in range(n_spectra):
        try:
            native_id = indexed_exp.getSpectrum(i).getNativeID()
        except Exception:
            continue
        id_to_index[native_id] = i
        num_match = ANY_DIGITS_RE.search(native_id)
        if num_match:
            id_to_index[f"_num_{num_match.group(1)}"] = i


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
    scan_match = NATIVE_ID_NUM_RE.search(clean_id)
    if scan_match:
        num = int(scan_match.group(1))

        if id_format == "index":
            arr_pos = num - index_offset
            if 0 <= arr_pos < n_spectra:
                return arr_pos

        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

    # Strategy 3: bare number fallback
    num_match = ANY_DIGITS_RE.search(clean_id)
    if num_match:
        num = int(num_match.group(1))
        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

        if id_format == "index":
            arr_pos = num - index_offset
            if 0 <= arr_pos < n_spectra:
                return arr_pos

    return None
