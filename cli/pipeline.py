"""
CLI pipeline phases — the five-step analysis run.

Split out of ``cli.main`` so ``main.py`` can stay focused on argparse wiring.
Each ``_run_*`` helper is a phase of :func:`run_pipeline`; export-related
phases live in ``cli.exports``. They all share ``CLIConfig`` and the
progress/printing helpers from ``cli.output``.
"""
from __future__ import annotations

import copy
import logging
import sys
import time
from pathlib import Path

from cli.config import CLIConfig
from cli.exports import export_dataframes, export_plots, prepare_export_dir
from cli.output import TerminalProgress, print_completion_summary, print_config_summary
from mspecviz.dataframes import build_psm_table
from mspecviz.pipeline import process_query_result
from mspecviz.search import search_peptides_in_idxmls
from mspecviz.spectra import load_spectra_streaming

logger = logging.getLogger("mspecviz.cli")


def _validate_inputs(config: CLIConfig) -> int:
    """Return exit code 0 if inputs are valid, 2 otherwise."""
    if not config.idxml_files:
        print("Error: No idXML files found. Check your config or --idxml argument.", file=sys.stderr)
        return 2
    if not config.mzml_files:
        print("Error: No mzML files found. Check your config or --mzml argument.", file=sys.stderr)
        return 2

    for label, files in [("idXML", config.idxml_files), ("mzML", config.mzml_files)]:
        for f in files:
            if not Path(f).exists():
                print(f"Error: {label} file not found: {f}", file=sys.stderr)
                return 2
    return 0


def _run_search(config: CLIConfig, quiet: bool):
    """Phase 1: peptide search. Returns (variant_groups, total_psms)."""
    if not quiet:
        print("Step 1: Searching peptides in idXML files...")

    progress = TerminalProgress("Searching", len(config.idxml_files), quiet)
    errors = []

    variant_groups = search_peptides_in_idxmls(
        idxml_paths=config.idxml_files,
        target_peptides=config.peptides,
        mzml_paths=config.mzml_files,
        progress_callback=lambda current, total, message: progress.update(current, message),
        error_callback=lambda message: (errors.append(message), logger.warning(message)),
        workers=min(config.workers, len(config.idxml_files)),
    )
    progress.finish()

    if errors and not quiet:
        for error in errors:
            print(f"  Warning: {error}")

    total_psms = sum(len(g.matches) for g in variant_groups.values())
    if not quiet:
        print(f"  Found {total_psms} PSMs across {len(variant_groups)} variant(s)\n")

    return variant_groups, total_psms


def _run_load(config: CLIConfig, variant_groups_copy, psm_lookup, total_psms: int, quiet: bool):
    """Phase 2: stream spectra + annotate. Returns (annotated_queries, reference_store, comparison_rows, normalize)."""
    if not quiet:
        print("Step 2: Loading and annotating spectra...")

    normalize = config.intensity_mode == "Relative"
    annotation_kwargs = dict(
        ion_types=config.ion_types,
        mz_tolerance=config.mz_tolerance,
        normalize=normalize,
        add_neutral_losses=config.add_neutral_losses,
        add_first_prefix_ion=config.add_first_prefix_ion,
        remove_precursor_peaks=config.remove_precursor_peaks,
    )

    annotated_queries = []   # (psm_id, LoadedSpectrum, annotated_df)
    reference_store = {}     # variant_key -> [(LoadedSpectrum, annotated_df)]
    comparison_rows = []

    progress = TerminalProgress("Loading", total_psms, quiet)
    loaded_count = [0]

    def on_query_ready(result):
        processed = process_query_result(
            result, psm_lookup, annotation_kwargs, reference_store,
            koina_enabled=config.koina_enabled,
            koina_model=config.koina_model,
            koina_collision_energy=config.koina_collision_energy,
            koina_instrument_type=config.koina_instrument_type,
        )
        if processed is None:
            return
        psm_id, query_spec, query_df, comparison_rows_for_psm = processed
        annotated_queries.append((psm_id, query_spec, query_df))
        comparison_rows.extend(comparison_rows_for_psm)
        loaded_count[0] += 1
        progress.update(loaded_count[0])

    load_spectra_streaming(
        variant_groups=variant_groups_copy,
        mzml_paths=[Path(m) for m in config.mzml_files],
        reference_substring=config.reference_substring,
        on_reference_loaded=lambda _v, _c: None,
        on_query_ready=on_query_ready,
        on_file_start=lambda p, is_ref: logger.info(
            f"Processing {'reference' if is_ref else 'query'}: {Path(p).name}"),
        on_file_complete=lambda p, loaded, total: logger.info(
            f"Completed {Path(p).name}: {loaded}/{total} spectra"),
        error_callback=lambda e: logger.error(f"Streaming error: {e}"),
        workers=config.workers,
    )
    progress.finish()

    return annotated_queries, reference_store, comparison_rows, normalize


def run_pipeline(config: CLIConfig, quiet: bool = False) -> int:
    """Execute the full MSpecViz analysis pipeline. Returns exit code."""
    pipeline_start = time.time()

    return_code = _validate_inputs(config)
    if return_code != 0:
        return return_code

    print_config_summary(config, quiet)

    variant_groups, total_psms = _run_search(config, quiet)
    if total_psms == 0:
        print("Error: No PSMs found. Check your peptide list and idXML files.", file=sys.stderr)
        return 1

    variant_groups_copy = copy.deepcopy(variant_groups)
    psm_frame, psm_lookup = build_psm_table(variant_groups_copy)

    annotated_queries, reference_store, comparison_rows, normalize = _run_load(
        config, variant_groups_copy, psm_lookup, total_psms, quiet,
    )
    if not annotated_queries:
        print("Error: No spectra could be loaded or annotated.", file=sys.stderr)
        return 1

    if not quiet:
        print(f"  Annotated {len(annotated_queries)} spectra, "
              f"{len(comparison_rows)} comparison(s)\n")
        print("Step 3: Building DataFrames...")

    root_dir, png_dir, svg_dir, pdf_dir = prepare_export_dir(config)

    if not quiet:
        print("Step 4: Exporting results...")

    _, _, _, comparison_df = export_dataframes(
        annotated_queries, comparison_rows, psm_frame, root_dir, quiet,
    )

    if not config.no_plots:
        export_plots(config, annotated_queries, reference_store, variant_groups_copy,
                     comparison_df, root_dir, png_dir, svg_dir, pdf_dir,
                     normalize, quiet)

    print_completion_summary(
        variant_count=len(variant_groups),
        total_psms=total_psms,
        annotated_count=len(annotated_queries),
        comparison_count=len(comparison_rows),
        root_dir=root_dir,
        elapsed=time.time() - pipeline_start,
        quiet=quiet,
    )
    return 0
