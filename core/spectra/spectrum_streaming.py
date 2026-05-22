"""Streaming mzML loader — orchestration layer that drives spectrum extraction in
reference-first order so each query can be compared against its references the
moment it arrives.
The single public entry point is load_spectra_streaming: it partitions the
mzML files into references and queries, dispatches the picklable per-file worker
(spectra.mzml_workers._load_spectra_from_mzml) through helpers.run_tasks,
and replays the loaded spectra to caller-supplied callbacks. References are loaded
and bucketed by peptide variant first; only then are queries streamed, each one
handed to on_query_ready packaged as a StreamingResult carrying its
matching reference spectra.
Called by pipeline.streaming_orchestration, which wires the metric-computing
callback into on_query_ready; also exercised directly by the integration
tests. Consumes the LoadResult / LoadedSpectrum / StreamingResult
dataclasses from spectra.loaded_spectrum.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from psm_utils import PSM, PSMList

from core.helpers import run_tasks
from core.spectra.loaded_spectrum import LoadResult, LoadedSpectrum, StreamingResult
from core.spectra.mzml_workers import _load_spectra_from_mzml

logger = logging.getLogger(__name__)


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
    for variant_key, psm_list in variant_groups.items():
        for psm in psm_list:
            if psm.run and (path_key := resolve(psm)):
                mzml_to_matches.setdefault(path_key, []).append((psm, variant_key))
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
