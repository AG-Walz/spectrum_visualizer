"""
Spectrum loading from mzML files and Parquet caches.

Provides memory-efficient indexed access to spectra using PyOpenMS.
Uses pre-resolved mzML mappings from the search phase.
Supports fast Parquet-based loading for pre-converted datasets.
"""
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable

import numpy as np
import pyopenms as oms

from src.peptide_search import PeptideMatch, VariantGroup

logger = logging.getLogger(__name__)


@dataclass
class LoadingProgress:
    """Progress tracking object for spectrum loading."""
    files_done: int = 0
    total_files: int = 0
    spectra_loaded: int = 0
    spectra_expected: int = 0
    spectra_not_found: int = 0
    current_file: str = ""
    current_file_spectra: int = 0
    current_file_total: int = 0
    errors: List[str] = field(default_factory=list)

    def add_error(self, error: str):
        """Add an error message."""
        self.errors.append(error)

    def file_complete(self):
        """Mark current file as complete."""
        self.files_done += 1

    def spectrum_loaded(self):
        """Mark a spectrum as successfully loaded."""
        self.spectra_loaded += 1
        self.current_file_spectra += 1

    def spectrum_not_found(self):
        """Mark a spectrum as not found."""
        self.spectra_not_found += 1

    def start_file(self, filename: str, expected_spectra: int):
        """Start processing a new file."""
        self.current_file = filename
        self.current_file_spectra = 0
        self.current_file_total = expected_spectra

    @property
    def total_processed(self) -> int:
        """Total spectra processed (found + not found)."""
        return self.spectra_loaded + self.spectra_not_found

    @property
    def progress_fraction(self) -> float:
        """Overall progress as a fraction 0-1."""
        if self.spectra_expected == 0:
            return 0.0
        return self.total_processed / self.spectra_expected


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

    for variant_key, vg in variant_groups.items():
        for match in vg.matches:
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
    from src.path_utils import disambiguate_paths

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
                (match, vk) for match, vk in matches_for_file
                if (vk, match.spectrum_id) not in skip_ids
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

            _, _, loaded_list, lc, nf, errs = _load_spectra_from_mzml(
                mzml_path_str, display_name, matches_for_file, True
            )
            for variant_key, spec in loaded_list:
                reference_spectra[variant_key].append(spec)
            for err in errs:
                if error_callback:
                    error_callback(err)

            if on_file_complete:
                on_file_complete(display_name, lc, nf)
    else:
        # Parallel reference loading
        from concurrent.futures import ProcessPoolExecutor, as_completed

        logger.info(f"Parallel reference loading: {len(ref_file_tasks)} files, {effective_workers} workers")

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            future_map = {
                pool.submit(_load_spectra_from_mzml, mzml_path_str, dn, mff, True): dn
                for mzml_path_str, dn, mff in ref_file_tasks
            }

            for future in as_completed(future_map):
                try:
                    dn, _, loaded_list, lc, nf, errs = future.result()
                    for variant_key, spec in loaded_list:
                        reference_spectra[variant_key].append(spec)
                    for err in errs:
                        if error_callback:
                            error_callback(err)
                    if on_file_complete:
                        on_file_complete(dn, lc, nf)
                except Exception as e:
                    dn = future_map[future]
                    if error_callback:
                        error_callback(f"Worker crashed loading reference {dn}: {e}")

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

    def _dispatch_query_results(loaded_list, dn):
        """Dispatch loaded query spectra through on_query_ready callback (main process)."""
        for variant_key, loaded in loaded_list:
            vg = variant_groups[variant_key]
            current_idx = query_indices[variant_key]
            query_indices[variant_key] += 1

            if on_query_ready:
                result = StreamingResult(
                    query_spectrum=loaded,
                    reference_spectra=reference_spectra[variant_key],
                    variant_key=variant_key,
                    base_sequence=vg.base_sequence,
                    modification_key=vg.modification_key,
                    query_index=current_idx,
                    total_queries=query_counts[variant_key]
                )
                on_query_ready(result)

    if effective_workers <= 1 or len(query_file_tasks) <= 1:
        # Sequential query loading
        for mzml_path_str, display_name, matches_for_file in query_file_tasks:
            if on_file_start:
                on_file_start(display_name, False)

            _, _, loaded_list, lc, nf, errs = _load_spectra_from_mzml(
                mzml_path_str, display_name, matches_for_file, False
            )
            for err in errs:
                if error_callback:
                    error_callback(err)
            _dispatch_query_results(loaded_list, display_name)

            if on_file_complete:
                on_file_complete(display_name, lc, nf)
    else:
        # Parallel query loading — load files in parallel, dispatch callbacks sequentially
        from concurrent.futures import ProcessPoolExecutor, as_completed

        logger.info(f"Parallel query loading: {len(query_file_tasks)} files, {effective_workers} workers")

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            future_map = {
                pool.submit(_load_spectra_from_mzml, mzml_path_str, dn, mff, False): dn
                for mzml_path_str, dn, mff in query_file_tasks
            }

            for future in as_completed(future_map):
                try:
                    dn, _, loaded_list, lc, nf, errs = future.result()
                    for err in errs:
                        if error_callback:
                            error_callback(err)
                    _dispatch_query_results(loaded_list, dn)
                    if on_file_complete:
                        on_file_complete(dn, lc, nf)
                except Exception as e:
                    dn = future_map[future]
                    if error_callback:
                        error_callback(f"Worker crashed loading query {dn}: {e}")

    logger.info("Query streaming complete")

    return reference_spectra


# =============================================================================
# Parquet-based spectrum loading
# =============================================================================

# Parquet schema field names (matches tools/convert_mzml_to_parquet.py)
PARQUET_COLUMNS = [
    "source_mzml", "spectrum_id", "peptidoform", "base_sequence",
    "variant_key", "charge", "precursor_mz", "rt", "is_reference",
    "mz", "intensity", "score", "idxml_file",
]


def load_spectra_from_parquet(
    parquet_path: str,
    reference_substring: str = "",
    on_reference_loaded: Optional[Callable[[str, int], None]] = None,
    on_query_ready: Optional[Callable[["StreamingResult"], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None,
) -> Dict[str, List[LoadedSpectrum]]:
    """
    Load pre-converted spectra from a Parquet file.

    Replaces the entire mzML loading pipeline when Parquet is available.
    ~100x faster than mzML loading (no XML parsing, no index building).

    The Parquet file must have been created by tools/convert_mzml_to_parquet.py
    with the expected schema (PARQUET_COLUMNS).

    Args:
        parquet_path: Path to .parquet file
        reference_substring: Substring to identify reference files (used as
            fallback if is_reference column is missing)
        on_reference_loaded: Callback(variant_key, count) when references are ready
        on_query_ready: Callback(StreamingResult) for each query spectrum
        error_callback: Callback(error_message) for errors

    Returns:
        Dict of variant_key -> list of reference LoadedSpectrum objects
    """
    import pyarrow.parquet as pq

    logger.info(f"Loading spectra from Parquet: {parquet_path}")

    try:
        table = pq.read_table(str(parquet_path))
    except Exception as e:
        msg = f"Failed to read Parquet file: {e}"
        logger.error(msg)
        if error_callback:
            error_callback(msg)
        return {}

    df = table.to_pandas()
    logger.info(f"Parquet loaded: {len(df)} rows, columns: {list(df.columns)}")

    # Validate required columns
    required = {"source_mzml", "spectrum_id", "peptidoform", "base_sequence",
                "variant_key", "charge", "precursor_mz", "mz", "intensity"}
    missing = required - set(df.columns)
    if missing:
        msg = f"Parquet file missing required columns: {missing}"
        logger.error(msg)
        if error_callback:
            error_callback(msg)
        return {}

    # Determine reference/query flag
    if "is_reference" in df.columns:
        pass  # Use stored flag
    elif reference_substring:
        df["is_reference"] = df["source_mzml"].str.contains(
            reference_substring, case=False, na=False
        )
    else:
        df["is_reference"] = False

    # Build VariantGroup-like structures from the Parquet data.
    # Extract columns upfront to avoid per-row boxing overhead of iterrows().
    reference_spectra: Dict[str, List[LoadedSpectrum]] = {}
    query_spectra_by_variant: Dict[str, List[LoadedSpectrum]] = {}
    variant_info: Dict[str, Tuple[str, str]] = {}  # variant_key -> (base_sequence, mod_key)

    col_variant_key = df["variant_key"].to_numpy()
    col_base_seq    = df["base_sequence"].to_numpy()
    col_spec_id     = df["spectrum_id"].to_numpy()
    col_peptidoform = df["peptidoform"].to_numpy()
    col_source_mzml = df["source_mzml"].to_numpy()
    col_precursor_mz = df["precursor_mz"].to_numpy(dtype=float)
    col_charge      = df["charge"].to_numpy(dtype=int)
    col_is_ref      = df["is_reference"].to_numpy(dtype=bool)
    col_rt          = df["rt"].to_numpy(dtype=float) if "rt" in df.columns else np.zeros(len(df))
    col_mz          = df["mz"].to_numpy()
    col_int         = df["intensity"].to_numpy()

    for i in range(len(df)):
        variant_key = col_variant_key[i]
        base_seq    = col_base_seq[i]

        # Derive modification_key from variant_key
        if "|" in variant_key:
            mod_key = variant_key.split("|", 1)[1]
        else:
            mod_key = ""

        if variant_key not in variant_info:
            variant_info[variant_key] = (base_seq, mod_key)
            reference_spectra[variant_key] = []
            query_spectra_by_variant[variant_key] = []

        mz_data = col_mz[i]
        int_data = col_int[i]
        mz_arr  = np.array(mz_data,  dtype=np.float64) if mz_data  is not None else np.array([])
        int_arr = np.array(int_data, dtype=np.float64) if int_data is not None else np.array([])

        loaded = LoadedSpectrum(
            spectrum_id=col_spec_id[i],
            mz=mz_arr,
            intensity=int_arr,
            precursor_mz=col_precursor_mz[i],
            precursor_charge=int(col_charge[i]),
            peptidoform=col_peptidoform[i],
            base_sequence=base_seq,
            source_file=col_source_mzml[i],
            is_reference=bool(col_is_ref[i]),
            rt=float(col_rt[i]),
        )

        if loaded.is_reference:
            reference_spectra[variant_key].append(loaded)
        else:
            query_spectra_by_variant[variant_key].append(loaded)

    # Phase 1: Notify references
    if on_reference_loaded:
        for vk, refs in reference_spectra.items():
            if refs:
                on_reference_loaded(vk, len(refs))

    total_refs = sum(len(v) for v in reference_spectra.values())
    total_queries = sum(len(v) for v in query_spectra_by_variant.values())
    logger.info(f"Parquet: {total_refs} reference, {total_queries} query spectra across {len(variant_info)} variants")

    # Phase 2: Stream query spectra
    if on_query_ready:
        for vk, queries in query_spectra_by_variant.items():
            base_seq, mod_key = variant_info[vk]
            refs = reference_spectra.get(vk, [])
            for i, query in enumerate(queries):
                result = StreamingResult(
                    query_spectrum=query,
                    reference_spectra=refs,
                    variant_key=vk,
                    base_sequence=base_seq,
                    modification_key=mod_key,
                    query_index=i,
                    total_queries=len(queries),
                )
                on_query_ready(result)

    return reference_spectra


# =============================================================================
# Polars conversion helpers
# =============================================================================

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


def spectra_metadata_to_polars(
    spectra_with_ids: List[Tuple[int, "LoadedSpectrum"]],
) -> "pl.DataFrame":
    """
    Convert a list of (psm_id, LoadedSpectrum) tuples to a Polars metadata DataFrame.

    Contains spectrum-level information (no peak arrays).  Useful for debugging
    or for building the sequence DataFrame independently of parquet_builder.

    Args:
        spectra_with_ids: List of (psm_id, LoadedSpectrum) tuples.

    Returns:
        Polars DataFrame with columns:
        psm_id, spectrum_id, peptidoform, base_sequence,
        source_file, precursor_mz, precursor_charge, rt, is_reference
    """
    import polars as pl

    if not spectra_with_ids:
        return pl.DataFrame(schema={
            "psm_id": pl.Int32, "spectrum_id": pl.String, "peptidoform": pl.String,
            "base_sequence": pl.String, "source_file": pl.String,
            "precursor_mz": pl.Float32, "precursor_charge": pl.Int32,
            "rt": pl.Float32, "is_reference": pl.Boolean,
        })

    rows = [
        {
            "psm_id": psm_id,
            "spectrum_id": spec.spectrum_id,
            "peptidoform": spec.peptidoform,
            "base_sequence": spec.base_sequence,
            "source_file": spec.source_file,
            "precursor_mz": float(spec.precursor_mz),
            "precursor_charge": int(spec.precursor_charge),
            "rt": float(spec.rt),
            "is_reference": bool(spec.is_reference),
        }
        for psm_id, spec in spectra_with_ids
    ]

    return pl.DataFrame(rows).with_columns([
        pl.col("psm_id").cast(pl.Int32),
        pl.col("precursor_mz").cast(pl.Float32),
        pl.col("precursor_charge").cast(pl.Int32),
        pl.col("rt").cast(pl.Float32),
    ])
