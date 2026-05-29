"""Streaming mzML loader: per-spectrum lookup + per-file worker + orchestration.

Combines three layers in one file:

- mzML spectrum-ID resolution (`open_indexed_mzml_fast`, `find_spectrum_index`):
  builds an ID -> array-position index for a single mzML so callers can fetch a
  spectrum by its native ID (scan=N / index=N / spectrum=N).
- Per-file worker (`_load_spectra_from_mzml`): the unit of work farmed out by
  the streaming stage. Opens one mzML and extracts the spectra named by a batch
  of PSMs, returning a `LoadResult` that bundles the wrapped `LoadedSpectrum`
  objects, per-spectrum counts and any errors. Must stay picklable and free of
  process-global state — runs in a ProcessPoolExecutor worker.
- Orchestration (`load_spectra_streaming`): drives spectrum extraction in
  reference-first order so each query can be compared against its references
  the moment it arrives. References are loaded and bucketed by peptide variant
  first; only then are queries streamed, each one handed to `on_query_ready`
  packaged as a `StreamingResult` carrying its matching reference spectra.

Called by core.pipeline.streaming_orchestration, which wires the
metric-computing callback into `on_query_ready`; also exercised directly by
the integration tests.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import pyopenms as oms
from psm_utils import PSM, PSMList

from core.helpers import run_tasks
from core.spectra.loaded_spectrum import (
    LoadResult,
    LoadedSpectrum,
    StreamingResult,
)
from core.spectra.mzml_native_ids import (
    INDEX_TAG,
    SCAN_TAG,
    SPECTRUM_TAG,
    clean_spectrum_id,
    extract_id_number,
)

logger = logging.getLogger(__name__)


# ============================================================================
# mzML spectrum-ID resolution: native ID -> array position
# ============================================================================

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
    needed_ids: Optional[List[str]] = None,
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """Open one mzML and build its native-ID -> array-position map. Returns (experiment, id_to_index).
    Three strategies, cheapest first: (1) for the index format the map needs no per-spectrum entries
    since find_spectrum_index derives positions from the offset directly;
    (2) when needed_ids is supplied for a scan / spectrum file, only those IDs are resolved;
    (3) otherwise every spectrum is walked.
    """
    # Open the file on disc and probe its native-ID format and index offset.
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(mzml_path)
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


# ============================================================================
# Per-file mzML worker (picklable; runs in ProcessPoolExecutor)
# ============================================================================

def _build_loaded_spectrum(
    spectrum, psm: PSM, display_name: str, is_reference: bool,
) -> LoadedSpectrum:
    """Wrap a pyopenms MSSpectrum and its matched PSM into a single LoadedSpectrum.
    Pulls the peak arrays, precursor m/z and charge, retention time, and
    peptidoform off the raw objects. When the spectrum carries no precursor
    block — or its charge is unset — the precursor charge falls back to the
    PSM's charge, and precursor_mz defaults to 0.0. Returns the fully
    populated LoadedSpectrum; runs inside a worker process, so the result
    must stay picklable.
    """
    mz_array, intensity_array = spectrum.get_peaks()
    # PSM charge is the fallback precursor charge when the spectrum lacks one
    psm_charge = psm.get_precursor_charge() or 0
    # Precursor: pull m/z and charge from the spectrum, falling back to psm_charge
    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        precursor_mz = prec.getMZ()
        precursor_charge = prec.getCharge() if prec.getCharge() > 0 else psm_charge
    else:
        precursor_mz = 0.0
        precursor_charge = psm_charge
    return LoadedSpectrum(
        spectrum_id=psm.spectrum_id,
        mz=mz_array,
        intensity=intensity_array,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        peptidoform=psm.peptidoform.modified_sequence,
        base_sequence=psm.peptidoform.sequence.upper(),
        source_file=display_name,
        is_reference=is_reference,
        rt=spectrum.getRT(),
        match=psm,
    )


def _load_spectra_from_mzml(
    mzml_path_str: str,
    display_name: str,
    matches_for_file: List[Tuple[PSM, str]],
    is_reference: bool,
) -> LoadResult:
    """Worker entry point — open one mzML and load every spectrum named in matches_for_file.
    matches_for_file is a list of (PSM, variant_key) pairs; each PSM's
    spectrum_id is resolved to an array position via a targeted index, the
    spectrum is fetched, and it is wrapped into a LoadedSpectrum tagged with
    its variant key. Errors are collected per spectrum (and once for the whole
    file if it fails to open) rather than raised — a failure on one spectrum or
    file never aborts the batch.
    Returns a LoadResult carrying loaded as a list of
    (variant_key, LoadedSpectrum) pairs, plus the loaded / not-found counts
    and the error list. Runs inside a ProcessPoolExecutor worker, so the
    return value must remain picklable.
    """
    # Per-file accumulators — successfully loaded spectra and error strings
    loaded: List[Tuple[str, LoadedSpectrum]] = []
    errors: List[str] = []

    try:
        # Open the file once and build an ID -> array-position map limited to
        # the spectra this batch actually needs — avoids a full-file scan.
        needed_ids = [psm.spectrum_id for psm, _ in matches_for_file]
        indexed_exp, id_to_index = open_indexed_mzml_fast(mzml_path_str, needed_ids=needed_ids)
        n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

        # Resolve each PSM's spectrum ID to a position, fetch it, and wrap it
        for psm, variant_key in matches_for_file:
            try:
                spectrum_index = find_spectrum_index(psm.spectrum_id, id_to_index, n_spectra)
                # Missing ID: record and skip — one absent spectrum is not fatal
                if spectrum_index is None:
                    errors.append(f"Spectrum not found: {psm.spectrum_id} in {display_name}")
                    continue
                spectrum = indexed_exp.getSpectrum(spectrum_index)
                loaded.append((variant_key, _build_loaded_spectrum(
                    spectrum, psm, display_name, is_reference,
                )))
            # Per-spectrum failure: keep going so the rest of the batch loads
            except Exception as error:
                errors.append(f"Error loading {psm.spectrum_id}: {error}")

    # Whole-file failure (e.g. unreadable mzML): record and return an empty result
    except Exception as error:
        errors.append(f"Error opening {display_name}: {error}")

    # Bundle counts and errors into a LoadResult for the streaming orchestrator
    return LoadResult(
        display_name=display_name,
        is_reference=is_reference,
        loaded=loaded,
        loaded_count=len(loaded),
        not_found_count=len(matches_for_file) - len(loaded),
        errors=errors,
    )


# ============================================================================
# Streaming orchestration: reference-first dispatch, query callback streaming
# ============================================================================

def build_mzml_to_matches_map(
    variant_groups: Dict[str, PSMList],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Group PSMs by the mzML file they live in; return {path_key: [(psm, variant_key), ...]}.
    path_key is the resolved, lower-cased absolute path of an mzML file. Each
    PSM is filed under exactly one file: the resolver prefers the pre-computed
    psm.metadata["mzml_resolved_path"] and falls back to matching psm.run
    against the mzML filename stem. PSMs that resolve to no available file (or
    carry an empty psm.run) are silently dropped.
    """
    # Build stem index: cheap reverse-lookup keyed on the mzML filename stem
    stem_to_key = {path.stem.lower(): key for key, path in available_mzml_paths.items()}

    # Resolver: prefer the PSM's pre-resolved path, fall back to stem match
    def resolve(psm: PSM) -> Optional[str]:
        resolved = psm.metadata.get("mzml_resolved_path", "")
        if resolved:
            key = str(Path(resolved).resolve()).lower()
            if key in available_mzml_paths:
                return key
        return stem_to_key.get(psm.run.lower())

    # Bucket: file each PSM under its resolved mzML path key
    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]] = {}
    dropped_examples: list[str] = []
    dropped_count = 0
    for variant_key, psm_list in variant_groups.items():
        for psm in psm_list:
            if psm.run and (path_key := resolve(psm)):
                mzml_to_matches.setdefault(path_key, []).append((psm, variant_key))
            else:
                dropped_count += 1
                if len(dropped_examples) < 5:
                    resolved = psm.metadata.get("mzml_resolved_path", "")
                    key = str(Path(resolved).resolve()).lower() if resolved else ""
                    dropped_examples.append(
                        f"psm.run={psm.run!r} resolved_path={resolved!r} "
                        f"key={key!r} in_map={key in available_mzml_paths} "
                        f"stem_hit={stem_to_key.get((psm.run or '').lower())!r}"
                    )
    if dropped_count:
        logger.warning(
            "build_mzml_to_matches_map dropped %d/%d PSMs. "
            "Available stems sample: %s. First dropped: %s",
            dropped_count,
            sum(len(g) for g in variant_groups.values()),
            list(stem_to_key.keys())[:3],
            dropped_examples,
        )
    return mzml_to_matches


def _build_file_tasks(file_keys, mzml_to_matches, mzml_path_map):
    """Pack (path_str, display_name, matches) task tuples for every file in
    file_keys that has at least one PSM match.
    file_keys is the references or queries subset; files with no matches are
    skipped so the worker is never dispatched on a file it has nothing to do in.
    The resulting tuples are the per-file work units handed to the mzML worker.
    """
    return [
        (str(mzml_path_map[key]), mzml_path_map[key].name, mzml_to_matches[key])
        for key in file_keys if key in mzml_to_matches
    ]


def _execute_load_tasks(file_tasks, is_reference, workers,
                        on_file_start, on_file_complete, error_callback) -> List[LoadResult]:
    """Run the mzML worker over file_tasks (sequential or parallel) and return
    a list of per-file LoadResult objects, firing the progress callbacks along the way.
    workers selects sequential (1) versus ProcessPoolExecutor dispatch
    via run_tasks. A crashed worker is reported through error_callback and
    its file is skipped rather than aborting the whole run; per-spectrum errors
    inside a successful result are forwarded the same way. is_reference is
    stamped onto every task so the worker knows which kind of file it is loading.
    """
    # Setup: short-circuit on empty, wire optional callbacks
    if not file_tasks:
        return []
    kind = "reference" if is_reference else "query"
    on_start = on_file_start or (lambda *_: None)
    on_complete = on_file_complete or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Build & announce: stamp is_reference onto each task; fire start callbacks
    tasks = [(p, d, m, is_reference) for p, d, m in file_tasks]
    for _path, display_name, _matches, _is_ref in tasks:
        on_start(display_name, is_reference)

    # Dispatch & collect: surface worker crashes, accumulate per-file results
    all_results: List[LoadResult] = []
    for task, ok, result in run_tasks(
        tasks=tasks, worker=_load_spectra_from_mzml, workers=workers,
    ):
        if not ok:
            on_error(f"Worker crashed loading {kind} {task[1]}: {result}")
            continue
        for error in result.errors:
            on_error(error)
        all_results.append(result)
        on_complete(result.display_name, result.loaded_count, result.not_found_count)
    return all_results


def _split_reference_and_query_keys(
    mzml_path_map: Dict[str, Path], reference_substring: str,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """Split the mzML path map into (references, queries) by case-insensitive
    substring match against the full path string.
    A file lands in references when reference_substring appears anywhere
    in its path; everything else is a query. An empty or missing
    reference_substring short-circuits the test so every file is a query —
    the no-experimental-reference case, where the pipeline later falls back to
    Koina-predicted references.
    """
    needle = reference_substring.lower() if reference_substring else None
    references: Dict[str, Path] = {}
    queries: Dict[str, Path] = {}
    for path_key, path in mzml_path_map.items():
        bucket = references if needle and needle in str(path).lower() else queries
        bucket[path_key] = path
    return references, queries


def _filter_already_loaded(
    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]],
    skip_ids: set,
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Drop matches whose (variant_key, spectrum_id) pair is in skip_ids,
    returning the trimmed mzml_to_matches map.
    Powers resuming a partial run: spectra already persisted in an earlier pass
    are skipped so only the remainder is reloaded. Files left with no surviving
    matches are dropped entirely (the walrus guard keeps only non-empty lists).
    Logs how many spectra were skipped versus how many remain.
    """
    before = sum(len(matches) for matches in mzml_to_matches.values())
    filtered = {
        path_key: kept for path_key, matches in mzml_to_matches.items()
        if (kept := [(psm, vk) for psm, vk in matches if (vk, psm.spectrum_id) not in skip_ids])
    }
    after = sum(len(matches) for matches in filtered.values())
    if before > after:
        logger.info(f"Resume: skipping {before - after} cached spectra, {after} remaining")
    return filtered


def _collect_reference_spectra(
    ref_results: List[LoadResult], variant_groups: Dict[str, PSMList],
    on_reference_loaded: Optional[Callable[[str, int], None]],
) -> Dict[str, List[LoadedSpectrum]]:
    """Bucket every loaded reference spectrum by variant_key and fire the
    per-variant callback once references are in; return {variant_key: [LoadedSpectrum]}.
    The returned map has one slot for every known variant — including variants
    that got zero reference spectra — so a later query lookup never raises a
    KeyError. on_reference_loaded is invoked only for variants that
    actually received references.
    """
    # Bucket: one slot per known variant, fan out the loaded spectra
    reference_spectra: Dict[str, List[LoadedSpectrum]] = {variant_key: [] for variant_key in variant_groups}
    for result in ref_results:
        for variant_key, spectrum in result.loaded:
            reference_spectra[variant_key].append(spectrum)

    # Notify: per-variant callback only fires for variants that actually got refs
    notify = on_reference_loaded or (lambda *_: None)
    for variant_key, references in reference_spectra.items():
        if references:
            notify(variant_key, len(references))

    logger.info(f"References loaded: {sum(len(r) for r in reference_spectra.values())} spectra")
    return reference_spectra


def _stream_query_results(
    query_results: List[LoadResult],
    reference_spectra: Dict[str, List[LoadedSpectrum]],
    on_query_ready: Optional[Callable[[StreamingResult], None]],
) -> None:
    """Fire on_query_ready once per loaded query spectrum, each call carrying a
    StreamingResult that pairs the query with its variant's reference spectra.
    Returns nothing — the loaded spectra reach the caller purely through the
    callback. When no callback is supplied the query results are simply
    discarded, since this function exists only to replay them downstream.
    """
    # Nothing to replay to: bail before touching any results
    if not on_query_ready:
        return
    for result in query_results:
        for variant_key, loaded in result.loaded:
            on_query_ready(StreamingResult(
                query_spectrum=loaded,
                reference_spectra=reference_spectra[variant_key],
                variant_key=variant_key,
            ))


def load_spectra_streaming(
    variant_groups: Dict[str, PSMList],
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
    """Load references first, then stream each query through on_query_ready;
    return {variant_key: [reference_spectra]}.
    The public entry point of the streaming loader. reference_substring
    (case-insensitive) partitions mzml_paths into reference and query files;
    references are loaded and bucketed by peptide variant up front so that each
    query, on arrival, can already be compared against its references. Spectra
    are surfaced through the optional callbacks rather than collected in bulk:
    on_query_ready receives a StreamingResult per query, the on_file_*
    and error_callback hooks drive progress reporting. workers controls
    sequential versus parallel per-file dispatch; skip_ids resumes a partial
    run by dropping already-loaded spectra. The returned reference map is the
    only bulk result — every query spectrum has already flowed through the
    callback by the time this returns.
    """
    # Partition inputs: resolve paths to a canonical key map, split into
    # references/queries, then group every PSM under the mzML it belongs to
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}
    reference_keys, query_keys = _split_reference_and_query_keys(mzml_path_map, reference_substring)

    mzml_to_matches = build_mzml_to_matches_map(variant_groups, mzml_path_map)
    # Resume support: trim matches already persisted by an earlier pass
    if skip_ids:
        mzml_to_matches = _filter_already_loaded(mzml_to_matches, skip_ids)

    # References first: load + bucket so queries can be compared on arrival
    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map)
    ref_results = _execute_load_tasks(
        ref_tasks, True, workers, on_file_start, on_file_complete, error_callback,
    )
    reference_spectra = _collect_reference_spectra(ref_results, variant_groups, on_reference_loaded)

    # Queries next: stream each loaded query through on_query_ready
    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map)
    query_results = _execute_load_tasks(
        query_tasks, False, workers, on_file_start, on_file_complete, error_callback,
    )
    _stream_query_results(query_results, reference_spectra, on_query_ready)

    logger.info("Query streaming complete")
    return reference_spectra
