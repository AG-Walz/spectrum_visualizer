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


def _validate_inputs(cfg: CLIConfig) -> int:
    """Return exit code 0 if inputs are valid, 2 otherwise."""
    if not cfg.idxml_files:
        print("Error: No idXML files found. Check your config or --idxml argument.", file=sys.stderr)
        return 2
    if not cfg.mzml_files:
        print("Error: No mzML files found. Check your config or --mzml argument.", file=sys.stderr)
        return 2

    for label, files in [("idXML", cfg.idxml_files), ("mzML", cfg.mzml_files)]:
        for f in files:
            if not Path(f).exists():
                print(f"Error: {label} file not found: {f}", file=sys.stderr)
                return 2
    return 0


def _run_search(cfg: CLIConfig, quiet: bool):
    """Phase 1: peptide search. Returns (variant_groups, total_psms)."""
    if not quiet:
        print("Step 1: Searching peptides in idXML files...")

    progress = TerminalProgress("Searching", len(cfg.idxml_files), quiet)
    errors = []

    variant_groups = search_peptides_in_idxmls(
        idxml_paths=cfg.idxml_files,
        target_peptides=cfg.peptides,
        mzml_paths=cfg.mzml_files,
        progress_callback=lambda c, t, msg: progress.update(c, msg),
        error_callback=lambda msg: (errors.append(msg), logger.warning(msg)),
        workers=min(cfg.workers, len(cfg.idxml_files)),
    )
    progress.finish()

    if errors and not quiet:
        for err in errors:
            print(f"  Warning: {err}")

    total_psms = sum(len(g.matches) for g in variant_groups.values())
    if not quiet:
        print(f"  Found {total_psms} PSMs across {len(variant_groups)} variant(s)\n")

    return variant_groups, total_psms


def _run_load(cfg: CLIConfig, vg_copy, psm_lookup, total_psms: int, quiet: bool):
    """Phase 2: stream spectra + annotate. Returns (annotated_queries, reference_store, comparison_rows, normalize)."""
    if not quiet:
        print("Step 2: Loading and annotating spectra...")

    normalize = cfg.intensity_mode == "Relative"
    annotation_kwargs = dict(
        ion_types=cfg.ion_types,
        mz_tolerance=cfg.mz_tolerance,
        normalize=normalize,
        add_neutral_losses=cfg.add_neutral_losses,
        add_first_prefix_ion=cfg.add_first_prefix_ion,
        remove_precursor_peaks=cfg.remove_precursor_peaks,
    )

    annotated_queries = []   # (psm_id, LoadedSpectrum, annotated_df)
    reference_store = {}     # variant_key -> [(LoadedSpectrum, annotated_df)]
    comparison_rows = []

    progress = TerminalProgress("Loading", total_psms, quiet)
    loaded_count = [0]

    def on_query_ready(result):
        processed = process_query_result(
            result, psm_lookup, annotation_kwargs, reference_store,
            koina_enabled=cfg.koina_enabled,
            koina_model=cfg.koina_model,
            koina_collision_energy=cfg.koina_collision_energy,
            koina_instrument_type=cfg.koina_instrument_type,
        )
        if processed is None:
            return
        psm_id, query_spec, query_df, comp_rows = processed
        annotated_queries.append((psm_id, query_spec, query_df))
        comparison_rows.extend(comp_rows)
        loaded_count[0] += 1
        progress.update(loaded_count[0])

    load_spectra_streaming(
        variant_groups=vg_copy,
        mzml_paths=[Path(m) for m in cfg.mzml_files],
        reference_substring=cfg.reference_substring,
        on_reference_loaded=lambda _v, _c: None,
        on_query_ready=on_query_ready,
        on_file_start=lambda p, is_ref: logger.info(
            f"Processing {'reference' if is_ref else 'query'}: {Path(p).name}"),
        on_file_complete=lambda p, loaded, total: logger.info(
            f"Completed {Path(p).name}: {loaded}/{total} spectra"),
        error_callback=lambda e: logger.error(f"Streaming error: {e}"),
        workers=cfg.workers,
    )
    progress.finish()

    return annotated_queries, reference_store, comparison_rows, normalize


def run_pipeline(cfg: CLIConfig, quiet: bool = False) -> int:
    """Execute the full MSpecViz analysis pipeline. Returns exit code."""
    pipeline_start = time.time()

    rc = _validate_inputs(cfg)
    if rc != 0:
        return rc

    print_config_summary(cfg, quiet)

    variant_groups, total_psms = _run_search(cfg, quiet)
    if total_psms == 0:
        print("Error: No PSMs found. Check your peptide list and idXML files.", file=sys.stderr)
        return 1

    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)

    annotated_queries, reference_store, comparison_rows, normalize = _run_load(
        cfg, vg_copy, psm_lookup, total_psms, quiet,
    )
    if not annotated_queries:
        print("Error: No spectra could be loaded or annotated.", file=sys.stderr)
        return 1

    if not quiet:
        print(f"  Annotated {len(annotated_queries)} spectra, "
              f"{len(comparison_rows)} comparison(s)\n")
        print("Step 3: Building DataFrames...")

    root_dir, png_dir, svg_dir, pdf_dir = prepare_export_dir(cfg)

    if not quiet:
        print("Step 4: Exporting results...")

    _, _, _, comparison_df = export_dataframes(
        annotated_queries, comparison_rows, psm_lf, root_dir, quiet,
    )

    if not cfg.no_plots:
        export_plots(cfg, annotated_queries, reference_store, vg_copy,
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
