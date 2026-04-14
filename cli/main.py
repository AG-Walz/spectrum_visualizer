#!/usr/bin/env python3
"""
MSpecViz CLI — Terminal-based analysis pipeline.

Runs the full MSpecViz analysis from the terminal without opening a browser.
All inputs are provided via a YAML config file (with optional CLI argument
overrides). Reuses the exact same backend functions as the Streamlit app;
config/parsing, progress/printing, and plot dispatch live in sibling modules
(``cli.config``, ``cli.output``, ``cli.plots``).

Usage:
    python -m cli --config config.yaml
    python -m cli --config config.yaml --peptides peptides.txt --verbose
    python -m cli --config config.yaml --mz-tolerance 0.05 --no-plots --quiet

Exit codes:
    0 = success, 1 = pipeline failure, 2 = configuration / input error
"""
import copy
import logging
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cli.config import CLIConfig, build_config, build_parser, load_peptides, load_yaml_config
from cli.output import TerminalProgress, print_completion_summary, print_config_summary
from cli.plots import export_mirror_plots
from mspecviz.export_dirs import commit_export_directory, find_or_create_export_directory
from mspecviz.insight_dataframes import (
    build_heatmap_dataframe, build_peaks_dataframe,
    build_psm_table, build_sequence_dataframe,
)
from mspecviz.peptide_search import search_peptides_in_idxmls
from mspecviz.query_processing import process_query_result
from mspecviz.spectrum_loader import load_spectra_streaming
from mspecviz.summary_plots import generate_all_summary_plots

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
    """Phase 2: stream spectra + annotate. Returns (annotated_queries, reference_store, comparison_rows)."""
    if not quiet:
        print("Step 2: Loading and annotating spectra...")

    normalize = cfg.intensity_mode == "Relative"
    annotation_kwargs = dict(
        ion_types=cfg.ion_types,
        mz_tolerance=cfg.mz_tolerance,
        normalize=normalize,
        add_losses=cfg.add_losses,
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
            koina_ce=cfg.koina_ce,
            koina_instrument=cfg.koina_instrument,
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


def _export_dataframes(annotated_queries, comparison_rows, psm_lf, root_dir: Path, quiet: bool):
    """Build peaks/seq/heatmap LazyFrames, write all four Parquets + comparison CSV/XLSX."""
    peaks_lf = build_peaks_dataframe([(pid, df) for pid, _, df in annotated_queries])
    seq_lf = build_sequence_dataframe([(pid, s) for pid, s, _ in annotated_queries])
    heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated_queries])
    comparison_df = pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()

    if not quiet:
        peaks_count = peaks_lf.select("x").collect().height
        print(f"  PSMs: {psm_lf.select('id').collect().height}, "
              f"Peaks: {peaks_count}, Comparisons: {len(comparison_df)}\n")

    for name, lf in [("psm", psm_lf), ("peaks", peaks_lf),
                     ("seq", seq_lf), ("heatmap", heatmap_lf)]:
        lf.collect().write_parquet(root_dir / f"{name}.parquet")

    if not quiet:
        print("  Saved 4 Parquet files")

    if not comparison_df.empty:
        comparison_df.to_csv(root_dir / "comparison_summary.csv", index=False)
        xlsx_saved = False
        try:
            comparison_df.to_excel(root_dir / "comparison_summary.xlsx",
                                   index=False, engine="openpyxl")
            xlsx_saved = True
        except ImportError:
            logger.warning("openpyxl not installed -- skipping XLSX export")
        except Exception as e:
            logger.warning(f"XLSX export failed: {e}")
        if not quiet:
            print(f"  Saved comparison summary (CSV{'/XLSX' if xlsx_saved else ''})")

    return peaks_lf, seq_lf, heatmap_lf, comparison_df


def run_pipeline(cfg: CLIConfig, quiet: bool = False) -> int:
    """Execute the full MSpecViz analysis pipeline. Returns exit code."""
    pipeline_start = time.time()

    rc = _validate_inputs(cfg)
    if rc != 0:
        return rc

    print_config_summary(cfg, quiet)

    # Phase 1: search peptides
    variant_groups, total_psms = _run_search(cfg, quiet)
    if total_psms == 0:
        print("Error: No PSMs found. Check your peptide list and idXML files.", file=sys.stderr)
        return 1

    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)

    # Phase 2: load + annotate spectra
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

    # Phase 3: prepare export directory
    run_settings = {
        "idxml_files": cfg.idxml_files, "mzml_files": cfg.mzml_files,
        "peptides": cfg.peptides, "reference_substring": cfg.reference_substring,
        "ion_types": cfg.ion_types, "mz_tolerance": cfg.mz_tolerance,
        "add_losses": cfg.add_losses, "add_first_prefix_ion": cfg.add_first_prefix_ion,
        "remove_precursor_peaks": cfg.remove_precursor_peaks,
        "intensity_mode": cfg.intensity_mode,
        "koina_enabled": cfg.koina_enabled, "koina_model": cfg.koina_model,
        "koina_ce": cfg.koina_ce, "koina_instrument": cfg.koina_instrument,
        "export_formats": cfg.export_formats, "dpi": cfg.dpi,
        "figsize": list(cfg.figsize),
    }
    root_dir, png_dir, svg_dir, pdf_dir = find_or_create_export_directory(
        cfg.export_dir, cfg.project_name, run_settings,
    )
    commit_export_directory(root_dir, run_settings)
    for fmt in cfg.export_formats:
        fmt_dir = {"png": png_dir, "svg": svg_dir, "pdf": pdf_dir}.get(fmt)
        if fmt_dir:
            fmt_dir.mkdir(parents=True, exist_ok=True)

    if not quiet:
        print("Step 4: Exporting results...")

    # Phase 4: write DataFrames + comparison summary
    _, _, _, comparison_df = _export_dataframes(
        annotated_queries, comparison_rows, psm_lf, root_dir, quiet,
    )

    # Phase 5: plots (optional)
    if not cfg.no_plots:
        n_plots = export_mirror_plots(
            annotated_queries, reference_store, vg_copy,
            cfg, root_dir, png_dir, svg_dir, pdf_dir,
            normalize, quiet,
        )
        if not quiet:
            print(f"  Generated {n_plots} plot file(s)")

        if not comparison_df.empty:
            summary_formats = [f for f in cfg.export_formats if f in ("png", "svg", "pdf")]
            if not summary_formats:
                summary_formats = ["png"]
            results = generate_all_summary_plots(
                source=comparison_df, output_dir=str(root_dir),
                formats=summary_formats,
            )
            n_summary = sum(len(v) for v in results.values())
            if not quiet:
                print(f"  Generated {n_summary} summary plot(s)")

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


def main():
    args = build_parser().parse_args()

    log_level = logging.ERROR if args.quiet else (
        logging.DEBUG if args.verbose else logging.WARNING)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    yaml_cfg = load_yaml_config(args.config)
    cfg = build_config(yaml_cfg, args)

    peptide_spec = args.peptides if args.peptides else yaml_cfg.get("peptides", [])
    if not peptide_spec:
        print("Error: No peptides specified. Use 'peptides' in config or --peptides.", file=sys.stderr)
        sys.exit(2)
    cfg.peptides, cfg.peptides_mode = load_peptides(
        peptide_spec, cfg.idxml_files, cfg.workers, args.quiet,
    )
    if not cfg.peptides:
        print("Error: Peptide list is empty.", file=sys.stderr)
        sys.exit(2)

    sys.exit(run_pipeline(cfg, quiet=args.quiet))


if __name__ == "__main__":
    main()
