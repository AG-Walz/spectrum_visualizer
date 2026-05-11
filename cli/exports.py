"""CLI export wrapper — builds frames + dispatches :func:`mspecviz.export.run_export`."""
from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from cli.config import CLIConfig
from cli.output import TerminalProgress
from mspecviz.dataframes import build_visualization_frames
from mspecviz.export import (
    ExportCallbacks,
    ExportResult,
    ExportSelection,
    PlotSettings,
    build_run_settings,
    run_export,
)

logger = logging.getLogger("mspecviz.cli")


def _make_callbacks(quiet: bool) -> ExportCallbacks:
    """Wire CLI TerminalProgress + stdout banners into the shared callback shape."""
    # Plot progress needs a TerminalProgress bar lazily — we only know ``total``
    # once ``run_tasks`` emits its first progress call.
    progress_holder: dict = {"bar": None}

    def _on_step(message: str) -> None:
        if not quiet:
            print(f"  {message}")

    def _on_plot_progress(completed: int, total: int) -> None:
        bar: Optional[TerminalProgress] = progress_holder["bar"]
        if bar is None:
            bar = TerminalProgress("Generating plots", total, quiet)
            progress_holder["bar"] = bar
        bar.update(completed)
        if completed >= total:
            bar.finish()
            progress_holder["bar"] = None

    def _on_warning(message: str) -> None:
        logger.warning(message)

    return ExportCallbacks(
        on_step=_on_step,
        on_plot_progress=_on_plot_progress,
        on_warning=_on_warning,
    )


def _print_export_preamble(frames: dict, comparison_frame) -> None:
    """Echo PSM/Peaks/Comparison row counts before kicking off the export."""
    peaks_height = frames["peaks"].select("x").collect().height
    psm_height = frames["psm"].select("id").collect().height
    print(f"  PSMs: {psm_height}, Peaks: {peaks_height}, "
          f"Comparisons: {len(comparison_frame)}\n")


def _build_export_selection(config: CLIConfig) -> ExportSelection:
    """Decide which artifact groups the CLI run will write."""
    return ExportSelection(
        parquet=True,
        csv=True,
        mirror_plots=not config.no_plots,
        summary_plots=not config.no_plots,
        html="html" in config.export_formats,
    )


def _build_plot_settings(config: CLIConfig) -> PlotSettings:
    """Translate the CLI config into a :class:`PlotSettings` for the runner."""
    return PlotSettings(
        export_formats=list(config.export_formats),
        figure_size=config.figure_size,
        export_resolution=config.export_resolution,
        plot_unpaired_queries=config.plot_unpaired_queries,
        workers=config.workers,
    )


def _print_export_summary(result: ExportResult) -> None:
    """Echo plot counts after the export completes."""
    print(f"  Generated {result.mirror_plot_count} plot file(s)")
    if result.summary_plot_count:
        print(f"  Generated {result.summary_plot_count} summary plot(s)")


def run_cli_export(
    config: CLIConfig,
    annotated_queries: list,
    reference_store: dict,
    variant_groups: dict,
    comparison_rows: list,
    psm_frame,
    quiet: bool,
) -> ExportResult:
    """Full CLI export — delegates to :func:`mspecviz.export.run_export`."""
    frames = build_visualization_frames(annotated_queries, psm_frame)
    comparison_frame = pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()

    if not quiet:
        _print_export_preamble(frames, comparison_frame)

    result = run_export(
        export_directory=config.export_directory,
        project_name=config.project_name,
        run_settings=build_run_settings(
            lambda key, default: getattr(config, key, default),
            peptides=config.peptides,
        ),
        frames=frames,
        annotated_queries=annotated_queries,
        reference_store=reference_store,
        variant_groups=variant_groups,
        comparison_frame=comparison_frame,
        selection=_build_export_selection(config),
        plot_settings=_build_plot_settings(config),
        callbacks=_make_callbacks(quiet),
    )

    if not quiet:
        _print_export_summary(result)
    return result
