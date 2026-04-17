"""CLI export helpers — dataframe writes, export-dir resolution, plot dispatch."""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from cli.config import CLIConfig
from cli.plots import export_mirror_plots
from mspecviz.dataframes import (
    build_heatmap_dataframe,
    build_peaks_dataframe,
    build_sequence_dataframe,
)
from mspecviz.export_dirs import commit_export_directory, find_or_create_export_directory
from mspecviz.plotting import generate_all_summary_plots

logger = logging.getLogger("mspecviz.cli")


def export_dataframes(annotated_queries, comparison_rows, psm_frame, root_dir: Path, quiet: bool):
    """Build peaks/sequence/heatmap LazyFrames, write all four Parquets + comparison CSV/XLSX."""
    peaks_frame = build_peaks_dataframe([(psm_id, df) for psm_id, _, df in annotated_queries])
    sequence_frame = build_sequence_dataframe([(psm_id, s) for psm_id, s, _ in annotated_queries])
    heatmap_frame = build_heatmap_dataframe([s for _, s, _ in annotated_queries])
    comparison_frame = pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()

    if not quiet:
        peaks_count = peaks_frame.select("x").collect().height
        print(f"  PSMs: {psm_frame.select('id').collect().height}, "
              f"Peaks: {peaks_count}, Comparisons: {len(comparison_frame)}\n")

    for name, frame in [("psm", psm_frame), ("peaks", peaks_frame),
                        ("sequence", sequence_frame), ("heatmap", heatmap_frame)]:
        frame.collect().write_parquet(root_dir / f"{name}.parquet")

    if not quiet:
        print("  Saved 4 Parquet files")

    if not comparison_frame.empty:
        comparison_frame.to_csv(root_dir / "comparison_summary.csv", index=False)
        xlsx_saved = False
        try:
            comparison_frame.to_excel(root_dir / "comparison_summary.xlsx",
                                      index=False, engine="openpyxl")
            xlsx_saved = True
        except ImportError:
            logger.warning("openpyxl not installed -- skipping XLSX export")
        except Exception as e:
            logger.warning(f"XLSX export failed: {e}")
        if not quiet:
            print(f"  Saved comparison summary (CSV{'/XLSX' if xlsx_saved else ''})")

    return peaks_frame, sequence_frame, heatmap_frame, comparison_frame


def _build_run_settings(config: CLIConfig) -> dict:
    """Serialise ``config`` into the run_settings dict stored alongside exports."""
    return {
        "idxml_files": config.idxml_files, "mzml_files": config.mzml_files,
        "peptides": config.peptides, "reference_substring": config.reference_substring,
        "ion_types": config.ion_types, "mz_tolerance": config.mz_tolerance,
        "add_neutral_losses": config.add_neutral_losses,
        "add_first_prefix_ion": config.add_first_prefix_ion,
        "remove_precursor_peaks": config.remove_precursor_peaks,
        "intensity_mode": config.intensity_mode,
        "koina_enabled": config.koina_enabled, "koina_model": config.koina_model,
        "koina_collision_energy": config.koina_collision_energy,
        "koina_instrument_type": config.koina_instrument_type,
        "export_formats": config.export_formats, "export_resolution": config.export_resolution,
        "figure_size": list(config.figure_size),
        "summary_metric": config.summary_metric,
        "summary_plot_types": list(config.summary_plot_types),
    }


def prepare_export_dir(config: CLIConfig):
    """Resolve/create the versioned export directory and its format subfolders."""
    run_settings = _build_run_settings(config)
    root_dir, png_dir, svg_dir, pdf_dir = find_or_create_export_directory(
        config.export_directory, config.project_name, run_settings,
    )
    commit_export_directory(root_dir, run_settings)
    for fmt in config.export_formats:
        fmt_dir = {"png": png_dir, "svg": svg_dir, "pdf": pdf_dir}.get(fmt)
        if fmt_dir:
            fmt_dir.mkdir(parents=True, exist_ok=True)
    return root_dir, png_dir, svg_dir, pdf_dir


def export_plots(config: CLIConfig, annotated_queries, reference_store, variant_groups_copy,
                 comparison_frame, root_dir, png_dir, svg_dir, pdf_dir,
                 normalize: bool, quiet: bool) -> None:
    """Generate mirror + summary plots (skipped when ``--no-plots`` is set)."""
    plot_file_count = export_mirror_plots(
        annotated_queries, reference_store, variant_groups_copy,
        config, root_dir, png_dir, svg_dir, pdf_dir,
        normalize, quiet,
    )
    if not quiet:
        print(f"  Generated {plot_file_count} plot file(s)")

    if not comparison_frame.empty:
        summary_formats = [f for f in config.export_formats if f in ("png", "svg", "pdf")]
        if not summary_formats:
            summary_formats = ["png"]
        results = generate_all_summary_plots(
            source=comparison_frame, output_dir=str(root_dir),
            formats=summary_formats,
        )
        summary_plot_count = sum(len(v) for v in results.values())
        if not quiet:
            print(f"  Generated {summary_plot_count} summary plot(s)")
