"""Top-level streaming-annotation orchestrator for the load stage of the pipeline.
Drives core.spectra.load_spectra_streaming to stream every mzML, hands each
arriving query spectrum to core.pipeline.process_query_result for annotation
and metric scoring, and buffers the outcomes into three accumulators returned to
the caller: annotated_queries (per-PSM annotated spectra), reference_store
(experimental reference cache, populated as a side effect), and comparison_rows
(one metrics dict per query/reference pair).
(called by cli.main and streamlit_app.views.load_annotate; the returned
buffers are consumed downstream by core.dataframes to build the openms_insight LazyFrames.)
"""
from pathlib import Path
from typing import Callable, Optional

from core.pipeline.query_processing import process_query_result
from core.spectra import load_spectra_streaming


def make_query_ready_callback(
    *,
    psm_lookup: dict,
    annotation_kwargs: dict,
    reference_store: dict,
    koina_kwargs: dict,
    annotated_queries: list,
    comparison_rows: list,
    on_progress: Callable[[int], None],
):
    """Build the per-PSM callback fed to load_spectra_streaming as on_query_ready.
    The returned closure annotates each arriving query, appends results to the
    caller-owned annotated_queries and comparison_rows buffers, and ticks
    on_progress with the running count of successfully loaded PSMs. Returns the
    nested on_query_ready callable.
    """
    # Mutable-cell closure idiom: a one-element list whose single slot the nested
    # closure can mutate (a plain int would be read-only from inside on_query_ready).
    loaded_count = [0]

    def on_query_ready(streaming_result):
        """Annotate one streamed query spectrum and record it into the shared buffers.
        Invoked once per query spectrum by load_spectra_streaming. Delegates to
        process_query_result; a None result means the spectrum has no matching
        PSM in the lookup results, so it is silently skipped. On success, the per-PSM
        annotated spectrum is appended to annotated_queries, its metric rows are
        extended into comparison_rows, and the progress counter advances. Returns None.
        """
        # Annotate the query, resolve its references, and score every pair.
        processed = process_query_result(
            streaming_result, psm_lookup, annotation_kwargs, reference_store, **koina_kwargs,
        )
        # Skip spectra with no matching PSM in the lookup results.
        if processed is None:
            return
        # Buffer the annotated query and its per-pair metric rows for the caller.
        psm_id, query_psm, query_frame, rows = processed
        annotated_queries.append((psm_id, query_psm, query_frame))
        comparison_rows.extend(rows)
        # Advance the running count and report progress to the frontend.
        loaded_count[0] += 1
        on_progress(loaded_count[0])

    return on_query_ready


def run_streaming_annotation(
    *,
    variant_groups: dict,
    mzml_paths,
    reference_substring: str,
    psm_lookup: dict,
    annotation_kwargs: dict,
    koina_kwargs: dict,
    workers: int,
    on_progress: Optional[Callable[[int], None]] = None,
    on_file_start: Optional[Callable] = None,
    on_file_complete: Optional[Callable] = None,
    on_error: Optional[Callable[[str], None]] = None,
) -> tuple[list, dict, list]:
    """Stream-load and annotate every PSM, returning the three load-stage buffers.
    Wires a per-PSM callback into load_spectra_streaming and runs the full
    streaming load. The optional on_* callbacks let the frontend track progress,
    file lifecycle, and errors; any left None fall back to a no-op. Returns
    (annotated_queries, reference_store, comparison_rows) — respectively the
    (psm_id, query_psm, query_frame) tuples, the experimental reference
    cache keyed by variant, and the flat list of comparison-metric dicts.
    """
    # Init buffers: shared accumulators the on_query_ready callback writes into;
    # passing them by reference lets the callback mutate them as spectra arrive.
    annotated_queries: list = []
    reference_store: dict = {}
    comparison_rows: list = []
    # No-op stand-in for any optional callback the caller did not supply.
    noop = lambda *_: None

    # Build callback: per-PSM closure that annotates and appends to buffers
    on_query_ready = make_query_ready_callback(
        psm_lookup=psm_lookup,
        annotation_kwargs=annotation_kwargs,
        reference_store=reference_store,
        koina_kwargs=koina_kwargs,
        annotated_queries=annotated_queries,
        comparison_rows=comparison_rows,
        on_progress=on_progress or noop,
    )

    # Dispatch: stream every mzML, fanning lifecycle callbacks to the frontend
    load_spectra_streaming(
        variant_groups=variant_groups,
        mzml_paths=[Path(path) for path in mzml_paths],
        reference_substring=reference_substring,
        on_reference_loaded=noop,
        on_query_ready=on_query_ready,
        on_file_start=on_file_start or noop,
        on_file_complete=on_file_complete or noop,
        error_callback=on_error or noop,
        workers=workers,
    )

    return annotated_queries, reference_store, comparison_rows
