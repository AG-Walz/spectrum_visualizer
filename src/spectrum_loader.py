"""
Spectrum loading from mzML files and Parquet caches.

Provides memory-efficient indexed access to spectra using PyOpenMS.
Uses pre-resolved mzML mappings from the search phase.
Supports fast Parquet-based loading for pre-converted datasets.
"""
import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable

import numpy as np
import pyopenms as oms

from src.path_utils import disambiguate_paths
from src.peptide_search import PeptideMatch, VariantGroup

logger = logging.getLogger(__name__)


@dataclass
class LoadedSpectrum:
    """A loaded spectrum with metadata."""
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
    match: Optional[PeptideMatch] = None


def build_mzml_to_matches_map(
    variant_groups: Dict[str, VariantGroup],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PeptideMatch, str]]]:
    """
    Build a mapping of mzML path key -> list of (match, variant_key) tuples.

    Uses mzml_resolved_path when available (handles same-name files from
    different directories). Falls back to bare filename matching.

    Args:
        variant_groups: Dict of variant groups with matches
        available_mzml_paths: Dict mapping path key (lowered full path) -> Path object.
            Also supports bare-name keys for backward compatibility.

    Returns:
        Dict mapping path key -> list of (match, variant_key) tuples
    """
    # Build a bare-name fallback lookup: lowercase name -> list of path keys
    name_to_keys: Dict[str, List[str]] = {}
    for key, path in available_mzml_paths.items():
        bare = path.name.lower()
        if bare not in name_to_keys:
            name_to_keys[bare] = []
        name_to_keys[bare].append(key)

    mzml_to_matches: Dict[str, List[Tuple[PeptideMatch, str]]] = {}

    for variant_key, variant_group in variant_groups.items():
        for match in variant_group.matches:
            if not match.mzml_filename:
                continue

            # Prefer resolved path (handles same-name files)
            path_key = None
            if match.mzml_resolved_path:
                resolved_key = str(Path(match.mzml_resolved_path).resolve()).lower()
                if resolved_key in available_mzml_paths:
                    path_key = resolved_key

            # Fallback: bare filename match (only if unambiguous)
            if path_key is None:
                bare_lower = match.mzml_filename.lower()
                candidates = name_to_keys.get(bare_lower, [])
                if len(candidates) == 1:
                    path_key = candidates[0]
                elif len(candidates) > 1:
                    # Ambiguous — pick first and warn
                    path_key = candidates[0]
                    candidate_paths = [str(available_mzml_paths[c]) for c in candidates]
                    logger.warning(
                        f"Ambiguous mzML match for '{match.mzml_filename}': "
                        f"{len(candidates)} files with same name. "
                        f"Using '{candidate_paths[0]}'. "
                        f"All candidates: {candidate_paths}"
                    )

            if path_key:
                if path_key not in mzml_to_matches:
                    mzml_to_matches[path_key] = []
                mzml_to_matches[path_key].append((match, variant_key))

    return mzml_to_matches



def _extract_id_number(native_id: str) -> Optional[int]:
    """Extract the numeric value from a native ID string."""
    m = re.search(r'(?:scan|index|spectrum)=(\d+)', native_id)
    if m:
        return int(m.group(1))
    m = re.search(r'(\d+)', native_id)
    if m:
        return int(m.group(1))
    return None


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """
    Open mzML with targeted (not full) index building.

    Optimizations over open_indexed_mzml():
    - For index=N format: O(1) lookup, no scanning at all
    - For scan=N with known needed_ids: only look up those spectra O(k)
    - Falls back to full scan for unknown formats

    Args:
        mzml_path: Path to mzML file
        needed_ids: Optional list of spectrum IDs we need to find.
                    If provided, enables targeted lookup instead of full scan.

    Returns:
        Tuple of (OnDiscMSExperiment, id_to_index_dict)
    """
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))

    n_spectra = indexed_exp.getNrSpectra()
    id_to_index = {}

    # Detect format and index offset from first spectrum
    id_format = "unknown"
    index_offset = 0
    if n_spectra > 0:
        first_id = indexed_exp.getSpectrum(0).getNativeID()
        if "index=" in first_id:
            id_format = "index"
            # Detect 0-based vs 1-based indexing
            first_num = _extract_id_number(first_id)
            if first_num is not None:
                index_offset = first_num  # 0 for 0-based, 1 for 1-based
        elif "scan=" in first_id:
            id_format = "scan"
        elif "spectrum=" in first_id:
            id_format = "spectrum"

    id_to_index["__format__"] = id_format
    id_to_index["__n_spectra__"] = n_spectra
    id_to_index["__index_offset__"] = index_offset

    # FAST PATH: For "index=N" format, apply offset to convert native ID number
    # to array position. E.g., 1-based "index=1" → array position 0.
    if id_format == "index":
        logger.debug(f"Fast path: index format detected (offset={index_offset}), "
                     f"skipping full scan ({n_spectra} spectra)")
        return indexed_exp, id_to_index

    # FAST PATH: If we know needed IDs and format is numeric, do targeted lookup
    if needed_ids and id_format in ("scan", "spectrum"):
        found = 0
        for sid in needed_ids:
            clean_id = clean_spectrum_id(sid)
            num = _extract_id_number(clean_id)
            if num is None:
                continue

            # Try direct index first (common case: scan numbers are near indices)
            candidates = []
            if 0 <= num < n_spectra:
                candidates.append(num)
            # Also try num-1 (1-based scan to 0-based index)
            if num > 0 and (num - 1) < n_spectra:
                candidates.append(num - 1)

            matched = False
            for candidate_idx in candidates:
                try:
                    check_spec = indexed_exp.getSpectrum(candidate_idx)
                    check_id = check_spec.getNativeID()
                    if clean_id == clean_spectrum_id(check_id):
                        id_to_index[clean_id] = candidate_idx
                        check_num = _extract_id_number(check_id)
                        if check_num is not None:
                            id_to_index[f"_num_{check_num}"] = candidate_idx
                        found += 1
                        matched = True
                        break
                except Exception:
                    continue

            if not matched:
                # Binary search: scan numbers are typically monotonically increasing
                lo, hi = 0, n_spectra - 1
                while lo <= hi:
                    mid = (lo + hi) // 2
                    try:
                        mid_spec = indexed_exp.getSpectrum(mid)
                        mid_id = mid_spec.getNativeID()
                        mid_num = _extract_id_number(mid_id)
                        if mid_num is None:
                            break
                        if mid_num == num:
                            id_to_index[clean_id] = mid
                            id_to_index[f"_num_{mid_num}"] = mid
                            found += 1
                            break
                        elif mid_num < num:
                            lo = mid + 1
                        else:
                            hi = mid - 1
                    except Exception:
                        break

        logger.debug(f"Targeted lookup: found {found}/{len(needed_ids)} spectra in {Path(mzml_path).name}")

        # If we found all needed IDs, skip full scan
        if found >= len(needed_ids):
            return indexed_exp, id_to_index

        # If we found most, return what we have (find_spectrum_index handles misses)
        if found > 0:
            return indexed_exp, id_to_index

    # SLOW FALLBACK: Full scan (unknown format or no needed_ids for scan/spectrum)
    logger.debug(f"Full scan fallback for {Path(mzml_path).name} ({n_spectra} spectra)")
    for i in range(n_spectra):
        try:
            spec = indexed_exp.getSpectrum(i)
            native_id = spec.getNativeID()
            id_to_index[native_id] = i
            num_match = re.search(r'(\d+)', native_id)
            if num_match:
                id_to_index[f"_num_{num_match.group(1)}"] = i
        except Exception:
            continue

    return indexed_exp, id_to_index


def clean_spectrum_id(spectrum_id: str) -> str:
    """Clean a spectrum ID for lookup."""
    if not spectrum_id:
        return ""

    if isinstance(spectrum_id, bytes):
        spectrum_id = spectrum_id.decode()

    clean_id = spectrum_id.strip()

    # If it contains "controllerType" (Thermo format), keep it
    if "controllerType" in clean_id:
        return clean_id

    # For simple scan=N format, ensure clean extraction
    for tag in ["scan=", "index=", "spectrum="]:
        if tag in clean_id:
            idx = clean_id.find(tag)
            rest = clean_id[idx:]
            parts = rest.split()
            return parts[0]

    return clean_id


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int
) -> Optional[int]:
    """
    Find spectrum index using the pre-built lookup.

    Args:
        spectrum_id: Native ID to find
        id_to_index: Dictionary mapping native IDs to indices
        n_spectra: Total number of spectra in file

    Returns:
        Spectrum index or None if not found
    """
    clean_id = clean_spectrum_id(spectrum_id)
    if not clean_id:
        return None

    id_format = id_to_index.get("__format__", "unknown")
    index_offset = id_to_index.get("__index_offset__", 0)

    # Strategy 1: Direct lookup
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    # Strategy 2: Extract number and try direct index access
    scan_match = re.search(r'(?:scan|index|spectrum)=(\d+)', clean_id)
    if scan_match:
        num = int(scan_match.group(1))

        # If format is "index", apply offset to convert native ID number
        # to array position. E.g., 1-based "index=1" → position 0.
        if id_format == "index":
            arr_pos = num - index_offset
            if 0 <= arr_pos < n_spectra:
                return arr_pos

        # Try numeric lookup
        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

    # Strategy 3: Try just the number
    num_match = re.search(r'(\d+)', clean_id)
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


def _load_spectra_from_mzml(
    mzml_path_str: str,
    display_name: str,
    matches_for_file: List[Tuple["PeptideMatch", str]],
    is_reference: bool,
) -> Tuple[str, bool, List[Tuple[str, "LoadedSpectrum"]], int, int, List[str]]:
    """
    Worker: open one mzML file and extract all matched spectra.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Args:
        mzml_path_str: Full path to mzML file
        display_name: Display name for logging
        matches_for_file: List of (PeptideMatch, variant_key) tuples
        is_reference: Whether this is a reference file

    Returns:
        Tuple of (display_name, is_reference, loaded_spectra, loaded_count,
                  not_found_count, errors)
        where loaded_spectra is a list of (variant_key, LoadedSpectrum) tuples
    """
    loaded: List[Tuple[str, LoadedSpectrum]] = []
    errors: List[str] = []
    loaded_count = 0
    not_found_count = 0

    try:
        needed_ids = [match.spectrum_id for match, _ in matches_for_file]
        indexed_exp, id_to_index = open_indexed_mzml_fast(mzml_path_str, needed_ids=needed_ids)
        n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

        for match, variant_key in matches_for_file:
            try:
                idx = find_spectrum_index(match.spectrum_id, id_to_index, n_spectra)

                if idx is None:
                    not_found_count += 1
                    errors.append(f"Spectrum not found: {match.spectrum_id} in {display_name}")
                    continue

                spectrum = indexed_exp.getSpectrum(idx)

                mz_array, intensity_array = spectrum.get_peaks()

                precursor_mz = 0.0
                precursor_charge = match.charge
                if spectrum.getPrecursors():
                    prec = spectrum.getPrecursors()[0]
                    precursor_mz = prec.getMZ()
                    if prec.getCharge() > 0:
                        precursor_charge = prec.getCharge()

                loaded_spec = LoadedSpectrum(
                    spectrum_id=match.spectrum_id,
                    mz=mz_array,
                    intensity=intensity_array,
                    precursor_mz=precursor_mz,
                    precursor_charge=precursor_charge,
                    peptidoform=match.peptidoform,
                    base_sequence=match.base_sequence,
                    source_file=display_name,
                    is_reference=is_reference,
                    rt=spectrum.getRT(),
                    match=match
                )

                loaded.append((variant_key, loaded_spec))
                loaded_count += 1

            except Exception as e:
                not_found_count += 1
                errors.append(f"Error loading {match.spectrum_id}: {e}")

    except Exception as e:
        errors.append(f"Error opening {display_name}: {e}")

    return display_name, is_reference, loaded, loaded_count, not_found_count, errors



@dataclass
class StreamingResult:
    """Result from streaming spectrum loader - one query spectrum with its matching references."""
    query_spectrum: LoadedSpectrum
    reference_spectra: List[LoadedSpectrum]
    variant_key: str
    base_sequence: str
    modification_key: str
    query_index: int  # Index of this query in its variant group
    total_queries: int  # Total queries expected for this variant


def load_spectra_streaming(
    variant_groups: Dict[str, VariantGroup],
    mzml_paths: List[Path],
    reference_substring: str = "",
    on_reference_loaded: Optional[Callable[[str, int], None]] = None,
    on_query_ready: Optional[Callable[[StreamingResult], None]] = None,
    on_file_start: Optional[Callable[[str, bool], None]] = None,
    on_file_complete: Optional[Callable[[str, int, int], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None,
    workers: int = 1,
    skip_ids: Optional[set] = None,
) -> Dict[str, List[LoadedSpectrum]]:
    """
    Load spectra with streaming: references first, then yield each query spectrum.

    This approach allows immediate processing of each query spectrum as it's loaded,
    enabling real-time visualization and export.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup with matches
        mzml_paths: List of available mzML paths
        reference_substring: Substring to identify reference files (case-insensitive)
        on_reference_loaded: Callback(variant_key, count) when all references for a variant are loaded
        on_query_ready: Callback(StreamingResult) for each query spectrum ready for processing
        on_file_start: Callback(filename, is_reference) when starting to process a file
        on_file_complete: Callback(filename, loaded_count, not_found_count) when file is done
        error_callback: Callback(error_message) for errors
        workers: Number of parallel workers for mzML loading. 1 = sequential (default).
        skip_ids: Optional set of (variant_key, spectrum_id) tuples to skip.
            Used for resuming from a partial cache — spectra already loaded are skipped.

    Returns:
        Dict of variant_key -> list of reference spectra (for summary purposes)
    """
    # Build available mzML lookup keyed by resolved lowercase path (collision-free)
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}

    # Compute disambiguated display names
    mzml_disambig = disambiguate_paths([str(p) for p in mzml_paths])

    # Separate reference and query files
    reference_files = {}
    query_files = {}

    for path_key, path in mzml_path_map.items():
        is_ref = reference_substring and reference_substring.lower() in str(path).lower()
        if is_ref:
            reference_files[path_key] = path
        else:
            query_files[path_key] = path

    # Build mapping using resolved paths
    mzml_to_matches = build_mzml_to_matches_map(variant_groups, mzml_path_map)

    # Filter out already-loaded spectra when resuming from partial cache
    skipped_count = 0
    if skip_ids:
        filtered_map = {}
        for path_key, matches_for_file in mzml_to_matches.items():
            remaining = [
                (match, variant_key) for match, variant_key in matches_for_file
                if (variant_key, match.spectrum_id) not in skip_ids
            ]
            skipped_count += len(matches_for_file) - len(remaining)
            if remaining:
                filtered_map[path_key] = remaining
        mzml_to_matches = filtered_map
        if skipped_count > 0:
            logger.info(f"Resume: skipping {skipped_count} already-cached spectra, "
                        f"{sum(len(v) for v in mzml_to_matches.values())} remaining")

    # Storage for reference spectra by variant
    reference_spectra: Dict[str, List[LoadedSpectrum]] = {k: [] for k in variant_groups.keys()}

    # Track query counts per variant for progress
    query_counts: Dict[str, int] = {k: 0 for k in variant_groups.keys()}
    for path_key, matches in mzml_to_matches.items():
        if path_key not in reference_files:
            for match, variant_key in matches:
                query_counts[variant_key] += 1

    # Determine effective worker count for file-level parallelism
    total_mzml_files = sum(1 for pk in list(reference_files) + list(query_files) if pk in mzml_to_matches)
    effective_workers = max(1, min(workers, total_mzml_files))

    # =========================================================================
    # Phase 1: Load all reference spectra
    # =========================================================================
    logger.info(f"Loading reference spectra from {len(reference_files)} files...")

    # Build list of reference file tasks
    ref_file_tasks = []
    for path_key in reference_files:
        if path_key not in mzml_to_matches:
            continue
        mzml_path = mzml_path_map[path_key]
        display_name = mzml_disambig.get(str(mzml_path), mzml_path.name)
        matches_for_file = mzml_to_matches[path_key]
        ref_file_tasks.append((str(mzml_path), display_name, matches_for_file))

    if effective_workers <= 1 or len(ref_file_tasks) <= 1:
        # Sequential reference loading
        for mzml_path_str, display_name, matches_for_file in ref_file_tasks:
            if on_file_start:
                on_file_start(display_name, True)

            _, _, loaded_list, loaded_count, not_found_count, errs = _load_spectra_from_mzml(
                mzml_path_str, display_name, matches_for_file, True
            )
            for variant_key, spec in loaded_list:
                reference_spectra[variant_key].append(spec)
            for err in errs:
                if error_callback:
                    error_callback(err)

            if on_file_complete:
                on_file_complete(display_name, loaded_count, not_found_count)
    else:
        # Parallel reference loading
        logger.info(f"Parallel reference loading: {len(ref_file_tasks)} files, {effective_workers} workers")

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            future_map = {
                pool.submit(_load_spectra_from_mzml, mzml_path_str, display_name, matches_for_file, True): display_name
                for mzml_path_str, display_name, matches_for_file in ref_file_tasks
            }

            for future in as_completed(future_map):
                try:
                    display_name, _, loaded_list, loaded_count, not_found_count, errs = future.result()
                    for variant_key, spec in loaded_list:
                        reference_spectra[variant_key].append(spec)
                    for err in errs:
                        if error_callback:
                            error_callback(err)
                    if on_file_complete:
                        on_file_complete(display_name, loaded_count, not_found_count)
                except Exception as e:
                    display_name = future_map[future]
                    if error_callback:
                        error_callback(f"Worker crashed loading reference {display_name}: {e}")

    # Notify that references are loaded for each variant
    if on_reference_loaded:
        for variant_key, refs in reference_spectra.items():
            if refs:
                on_reference_loaded(variant_key, len(refs))

    logger.info(f"Reference loading complete: {sum(len(v) for v in reference_spectra.values())} spectra")

    # =========================================================================
    # Phase 2: Stream query spectra one at a time
    # =========================================================================
    logger.info(f"Streaming query spectra from {len(query_files)} files...")

    # Track query index per variant
    query_indices: Dict[str, int] = {k: 0 for k in variant_groups.keys()}

    # Build list of query file tasks
    query_file_tasks = []
    for path_key in query_files:
        if path_key not in mzml_to_matches:
            continue
        mzml_path = mzml_path_map[path_key]
        display_name = mzml_disambig.get(str(mzml_path), mzml_path.name)
        matches_for_file = mzml_to_matches[path_key]
        query_file_tasks.append((str(mzml_path), display_name, matches_for_file))

    def _dispatch_query_results(loaded_list, file_display_name):
        """Dispatch loaded query spectra through on_query_ready callback (main process)."""
        for variant_key, loaded in loaded_list:
            variant_group = variant_groups[variant_key]
            current_idx = query_indices[variant_key]
            query_indices[variant_key] += 1

            if on_query_ready:
                result = StreamingResult(
                    query_spectrum=loaded,
                    reference_spectra=reference_spectra[variant_key],
                    variant_key=variant_key,
                    base_sequence=variant_group.base_sequence,
                    modification_key=variant_group.modification_key,
                    query_index=current_idx,
                    total_queries=query_counts[variant_key]
                )
                on_query_ready(result)

    if effective_workers <= 1 or len(query_file_tasks) <= 1:
        # Sequential query loading
        for mzml_path_str, display_name, matches_for_file in query_file_tasks:
            if on_file_start:
                on_file_start(display_name, False)

            _, _, loaded_list, loaded_count, not_found_count, errs = _load_spectra_from_mzml(
                mzml_path_str, display_name, matches_for_file, False
            )
            for err in errs:
                if error_callback:
                    error_callback(err)
            _dispatch_query_results(loaded_list, display_name)

            if on_file_complete:
                on_file_complete(display_name, loaded_count, not_found_count)
    else:
        # Parallel query loading — load files in parallel, dispatch callbacks sequentially
        logger.info(f"Parallel query loading: {len(query_file_tasks)} files, {effective_workers} workers")

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            future_map = {
                pool.submit(_load_spectra_from_mzml, mzml_path_str, display_name, matches_for_file, False): display_name
                for mzml_path_str, display_name, matches_for_file in query_file_tasks
            }

            for future in as_completed(future_map):
                try:
                    display_name, _, loaded_list, loaded_count, not_found_count, errs = future.result()
                    for err in errs:
                        if error_callback:
                            error_callback(err)
                    _dispatch_query_results(loaded_list, display_name)
                    if on_file_complete:
                        on_file_complete(display_name, loaded_count, not_found_count)
                except Exception as e:
                    display_name = future_map[future]
                    if error_callback:
                        error_callback(f"Worker crashed loading query {display_name}: {e}")

    logger.info("Query streaming complete")

    return reference_spectra


def spectrum_lookup_key(spec: "LoadedSpectrum") -> str:
    """
    Return a stable key for a LoadedSpectrum used to match it to a PSM table row.

    The key is (spectrum_id, peptidoform) which is unique per PSM: the same
    physical spectrum may appear with different peptidoforms in different PSMs,
    and spectrum IDs alone can collide across mzML files.

    Use this key together with build_psm_table()'s psm_lookup dict to resolve
    psm_id when building the peaks and sequence DataFrames in parquet_builder.
    """
    return f"{spec.spectrum_id}|{spec.peptidoform}"

