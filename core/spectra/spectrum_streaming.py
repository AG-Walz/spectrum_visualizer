"""Streaming mzML loader: per-file worker + reference-first orchestration.

Combines two layers in one file:

- Per-file worker (`_load_spectra_from_mzml`): the unit of work farmed out by
  the streaming stage. Opens one mzML and extracts the spectra named by a batch
  of PSMs — resolving each PSM's native ID with pyopenms'
  OnDiscMSExperiment.getSpectrumByNativeId (an indexed exact-string lookup) and
  falling back to SpectrumLookup.findByScanNumber for references that carry only
  the scan number. Returns a `LoadResult` that bundles the prepared PSM carrier
  objects, per-spectrum counts and any errors. Must stay picklable and free of
  process-global state — runs in a ProcessPoolExecutor worker.
- Orchestration (`load_spectra_streaming`): drives spectrum extraction in
  reference-first order so each query can be compared against its references
  the moment it arrives. References are loaded and bucketed by peptidoform
  first; only then are queries streamed, each one handed to `on_query_ready`
  packaged as a `StreamingResult` carrying its matching reference spectra.

Called by core.pipeline.load_annotate.run_streaming_annotation, which wires the
metric-computing callback into `on_query_ready`; also exercised directly by
the integration tests.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import pyopenms as oms
from psm_utils import PSM, PSMList

from core.helpers import run_tasks
from core.spectra.streaming_types import (
    LoadResult,
    StreamingResult,
)
from core.spectra.mzml_native_ids import extract_id_number

logger = logging.getLogger(__name__)


# ============================================================================
# Per-file mzML worker (picklable; runs in ProcessPoolExecutor)
# ============================================================================

# SpectrumLookup.readSpectra needs a "?<SCAN>" named group; mirror the keyed
# pattern extract_id_number uses so the fallback resolves on the same number.
SCAN_NUMBER_REGEX = r"(?:scan|index|spectrum)=(?<SCAN>\d+)"


def _prepare_loaded_psm(spectrum, psm: PSM) -> PSM:
    """Stamp the loaded spectrum's peaks and acquisition metadata onto the matched PSM
    and return it — the PSM is the carrier the pipeline streams. Attaches the peak arrays
    to psm.spectrum and pulls precursor m/z, retention time, and ion mobility (each None
    when absent). The precursor charge stays the PSM's identification charge (peptidoform).
    Runs inside a worker process, so the result must stay picklable.
    """
    mz_array, intensity_array = spectrum.get_peaks()
    precursors = spectrum.getPrecursors()
    precursor_mz = precursors[0].getMZ() if precursors else 0.0
    drift_time = spectrum.getDriftTime()

    psm.spectrum = {"mz": mz_array, "intensity": intensity_array}
    # precursor m/z is None when the spectrum carries no precursor block
    psm.precursor_mz = precursor_mz if precursor_mz > 0 else None
    psm.retention_time = spectrum.getRT()
    # ion mobility is None when unset — getDriftTime() returns -1.0 for non-IM data
    psm.ion_mobility = drift_time if drift_time and drift_time > 0 else None
    return psm


def _load_spectra_from_mzml(
    mzml_path_str: str,
    display_name: str,
    matches_for_file: List[Tuple[PSM, str]],
    is_reference: bool,
) -> LoadResult:
    """Worker entry point — open one mzML and load every spectrum named in matches_for_file.
    Each (PSM, peptidoform) pair has its native ID resolved to a spectrum, fetched, and prepared as a
    PSM carrier tagged with its peptidoform. Errors are collected per spectrum (and once per file if it
    fails to open) rather than raised, so one failure never aborts the batch. Returns a picklable
    LoadResult with the (peptidoform, PSM) pairs, counts, and errors.
    """
    # Per-file accumulators — successfully loaded PSM carriers and error strings
    peptidoform_psm_pairs: List[Tuple[str, PSM]] = []
    errors: List[str] = []

    try:
        # Open the file once on disc; spectra are fetched lazily by native ID below.
        indexed_exp = oms.OnDiscMSExperiment()
        indexed_exp.openFile(mzml_path_str)
        # Number-keyed fallback index, built lazily only when an exact match misses.
        scan_lookup: Optional[oms.SpectrumLookup] = None

        # Resolve each PSM's spectrum, fetch it, and wrap it as a carrier
        for psm, peptidoform in matches_for_file:
            try:
                # Exact path: pyopenms resolves the full native-ID string (scan= / index= /
                # spectrum= / raw) straight off the indexed mzML, raising RuntimeError on a miss.
                try:
                    spectrum = indexed_exp.getSpectrumByNativeId(psm.spectrum_id)
                # Fallback: the reference carries only the scan number, not the exact native-ID
                # string. Build a SpectrumLookup once (from metadata — no peaks) and match by number.
                except RuntimeError:
                    scan_number = extract_id_number(psm.spectrum_id)
                    if scan_number is None:
                        errors.append(f"Spectrum not found: {psm.spectrum_id} in {display_name}")
                        continue
                    if scan_lookup is None:
                        scan_lookup = oms.SpectrumLookup()
                        scan_lookup.readSpectra(indexed_exp.getMetaData(), SCAN_NUMBER_REGEX)
                    try:
                        spectrum = indexed_exp.getSpectrum(scan_lookup.findByScanNumber(scan_number))
                    except RuntimeError:
                        errors.append(f"Spectrum not found: {psm.spectrum_id} in {display_name}")
                        continue
                peptidoform_psm_pairs.append((peptidoform, _prepare_loaded_psm(spectrum, psm)))
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
        peptidoform_psm_pairs=peptidoform_psm_pairs,
        loaded_count=len(peptidoform_psm_pairs),
        not_found_count=len(matches_for_file) - len(peptidoform_psm_pairs),
        errors=errors,
    )


# ============================================================================
# Streaming orchestration: reference-first dispatch, query callback streaming
# ============================================================================

def build_mzml_to_matches_map(
    peptidoform_psms: Dict[str, PSMList],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Group PSMs by the mzML file they live in; return {path_key: [(psm, peptidoform), ...]}.
    path_key is an mzML's resolved, lower-cased absolute path. PSMs are bucketed by run, and each
    run is resolved to one available path — preferring psm.metadata["mzml_resolved_path"], else
    matching psm.run against the filename stem. PSMs whose run resolves to no available file (or
    carry an empty run) are dropped, with a warning logged for the dropped count.
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
        return stem_to_key.get((psm.run or "").lower())

    # Bucket PSMs by run (all of a run's PSMs share .run, so resolve once per run).
    by_run: Dict[str, List[PSM]] = {}
    total_psms = 0
    for psm_list in peptidoform_psms.values():
        for psm in psm_list:
            total_psms += 1
            by_run.setdefault(psm.run or "", []).append(psm)

    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]] = {}
    dropped_count = 0
    for run, run_psms in by_run.items():
        path_key = resolve(run_psms[0]) if run else None
        if path_key:
            # Pair each PSM with its modified-sequence peptidoform — the same key its
            # peptidoform_psms entry is keyed on.
            mzml_to_matches.setdefault(path_key, []).extend(
                (psm, psm.peptidoform.modified_sequence) for psm in run_psms
            )
        else:
            dropped_count += len(run_psms)
    if dropped_count:
        logger.warning(
            "build_mzml_to_matches_map dropped %d/%d PSMs that resolved to no available mzML file.",
            dropped_count,
            total_psms,
        )
    return mzml_to_matches


def _build_file_tasks(file_keys, mzml_to_matches, mzml_path_map):
    """Pack (path_str, display_name, matches) task tuples for every file in file_keys
    that has at least one PSM match. file_keys is the references or queries subset;
    files with no matches are skipped so the worker is never dispatched on an empty
    file. The tuples are the per-file work units handed to the mzML worker.
    """
    return [
        (str(mzml_path_map[key]), mzml_path_map[key].name, mzml_to_matches[key])
        for key in file_keys if key in mzml_to_matches
    ]


def _execute_load_tasks(file_tasks, is_reference, workers,
                        on_file_start, on_file_complete, error_callback) -> List[LoadResult]:
    """Run the mzML worker over file_tasks (sequential or parallel) and return the per-file
    LoadResult list, firing the progress callbacks along the way. workers selects sequential (1)
    versus ProcessPoolExecutor dispatch. A crashed worker is reported through error_callback and
    its file skipped rather than aborting the whole run; per-spectrum errors are forwarded the same way.
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
        # Forward this file's per-spectrum errors through the same callback
        for error in result.errors:
            on_error(error)
        all_results.append(result)
        on_complete(result.display_name, result.loaded_count, result.not_found_count)
    return all_results


def _split_reference_and_query_keys(
    mzml_path_map: Dict[str, Path], reference_substring: str,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """Split the mzML path map into (references, queries) by a case-insensitive substring
    match against the full path. A file lands in references when reference_substring appears
    anywhere in its path; everything else is a query. An empty reference_substring makes every
    file a query — the no-experimental-reference case, where the pipeline later falls back to Koina references.
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
    """Drop matches whose (peptidoform, spectrum_id) pair is in skip_ids, returning the trimmed
    mzml_to_matches map. Powers resuming a partial run: spectra already persisted in an earlier
    pass are skipped so only the remainder is reloaded. Files left with no surviving matches are
    dropped entirely, and the skipped-versus-remaining counts are logged.
    """
    before = sum(len(matches) for matches in mzml_to_matches.values())
    # The walrus guard keeps only files whose `kept` list is non-empty (still has matches).
    filtered = {
        path_key: kept for path_key, matches in mzml_to_matches.items()
        if (kept := [(psm, peptidoform) for psm, peptidoform in matches if (peptidoform, psm.spectrum_id) not in skip_ids])
    }
    after = sum(len(matches) for matches in filtered.values())
    if before > after:
        logger.info(f"Resume: skipping {before - after} cached spectra, {after} remaining")
    return filtered


def _collect_reference_spectra(
    ref_results: List[LoadResult], peptidoform_psms: Dict[str, PSMList],
    on_reference_loaded: Optional[Callable[[str, int], None]],
) -> Dict[str, List[PSM]]:
    """Bucket every loaded reference spectrum by peptidoform and fire the per-peptidoform callback;
    return {peptidoform: [PSM]}. The map has one slot for every known peptidoform — including those
    that got zero reference spectra — so a later query lookup never raises a KeyError.
    on_reference_loaded is invoked only for peptidoforms that actually received references.
    """
    # Bucket: one slot per known peptidoform, fan out the loaded spectra
    reference_spectra: Dict[str, List[PSM]] = {peptidoform: [] for peptidoform in peptidoform_psms}
    for load_result in ref_results:
        for peptidoform, reference_psm in load_result.peptidoform_psm_pairs:
            reference_spectra[peptidoform].append(reference_psm)

    # Notify: per-peptidoform callback only fires for peptidoforms that actually got refs
    notify = on_reference_loaded or (lambda *_: None)
    for peptidoform, references in reference_spectra.items():
        if references:
            notify(peptidoform, len(references))

    logger.info(f"References loaded: {sum(len(r) for r in reference_spectra.values())} spectra")
    return reference_spectra


def _stream_query_results(
    query_results: List[LoadResult],
    reference_spectra: Dict[str, List[PSM]],
    on_query_ready: Optional[Callable[[StreamingResult], None]],
) -> None:
    """Fire on_query_ready once per loaded query spectrum, each call carrying a
    StreamingResult that pairs the query with its peptidoform's reference spectra.
    The loaded spectra reach the caller purely through the callback; with no callback
    the results are discarded.
    """
    # Nothing to replay to: bail before touching any results
    if not on_query_ready:
        return
    for load_result in query_results:
        for peptidoform, query_psm in load_result.peptidoform_psm_pairs:
            on_query_ready(StreamingResult(
                query_psm=query_psm,
                reference_psms=reference_spectra[peptidoform],
                peptidoform=peptidoform,
            ))


def load_spectra_streaming(
    peptidoform_psms: Dict[str, PSMList],
    mzml_paths: List[Path],
    reference_substring: str = "",
    on_reference_loaded: Optional[Callable[[str, int], None]] = None,
    on_query_ready: Optional[Callable[[StreamingResult], None]] = None,
    on_file_start: Optional[Callable[[str, bool], None]] = None,
    on_file_complete: Optional[Callable[[str, int, int], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None,
    workers: int = 1,
    skip_ids: Optional[set] = None,
) -> Dict[str, List[PSM]]:
    """Load references first, then stream each query through on_query_ready; return
    {peptidoform: [reference_spectra]}. reference_substring (case-insensitive) partitions
    mzml_paths into reference and query files; references are bucketed by peptidoform up front so
    each query is compared against its references on arrival. Query spectra reach the caller only
    through the callbacks — the returned reference map is the one bulk result; workers sets dispatch, skip_ids resumes.
    """
    # Partition inputs: resolve paths to a canonical key map, split into
    # references/queries, then group every PSM under the mzML it belongs to
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}
    reference_keys, query_keys = _split_reference_and_query_keys(mzml_path_map, reference_substring)

    mzml_to_matches = build_mzml_to_matches_map(peptidoform_psms, mzml_path_map)
    # Resume support: trim matches already persisted by an earlier pass
    if skip_ids:
        mzml_to_matches = _filter_already_loaded(mzml_to_matches, skip_ids)

    # References first: load + bucket so queries can be compared on arrival
    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map)
    ref_results = _execute_load_tasks(
        ref_tasks, True, workers, on_file_start, on_file_complete, error_callback,
    )
    reference_spectra = _collect_reference_spectra(ref_results, peptidoform_psms, on_reference_loaded)

    # Queries next: stream each loaded query through on_query_ready
    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map)
    query_results = _execute_load_tasks(
        query_tasks, False, workers, on_file_start, on_file_complete, error_callback,
    )
    _stream_query_results(query_results, reference_spectra, on_query_ready)

    logger.info("Query streaming complete")
    return reference_spectra
