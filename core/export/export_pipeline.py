"""Top-level export orchestrator that turns analysis results into an on-disk run folder.
run_export is the single entry point: it plans (and commits) a versioned export
folder, writes the four parquet LazyFrames, the comparison table as CSV/XLSX, the
per-PSM mirror plots (rendered across a worker pool), and the dual-axis / ridgeline /
heatmap summary plots — returning an ExportResult that lists every saved file.
The ExportSelection / PlotSettings / ExportCallbacks dataclasses let each
frontend describe one export run without passing a long argument list.
(Called by cli.main._run_export and streamlit_app.views.export.export_actions;
re-exported by core.export. Dispatches to core.pipeline.build_plot_tasks,
core.plotting.process_single_task / generate_all_summary_plots, and the
folder helpers in core.export.export_directories.)
"""
from __future__ import annotations
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable
import pandas as pd

from core.helpers import run_tasks
from core.export.export_directories import (
    commit_export_directory,
    find_or_create_export_directory,
    sanitize_filename,
)
from core.pipeline import build_plot_tasks
from core.plotting import generate_all_summary_plots, process_single_task

logger = logging.getLogger(__name__)


# Two shapes a frontend may hand in: a list [(psm_id, spec, df), ...] or a
# dict {psm_id: (spec, df)} — run_export normalises both into a list.
AnnotatedQueries = list[tuple] | dict[str, tuple]
# The four openms_insight LazyFrames written one parquet file each.
_PARQUET_FRAMES = ("psm", "peaks", "sequence", "heatmap")
# Raster/vector image formats accepted for plots — distinct from html,
# which is handled on its own track because it is interactive, not a picture.
_PICTURE_FORMATS = ("png", "svg", "pdf")


@dataclass
class ExportSelection:
    """Which artefact groups to write (parquet / CSV / mirror plots / summaries / HTML).
    """
    parquet: bool = True
    csv: bool = True
    mirror_plots: bool = True
    summary_plots: bool = True
    html: bool = False


@dataclass
class PlotSettings:
    """Figure size, image formats, DPI, and plot-worker count.
    """
    export_formats: list[str] = field(default_factory=lambda: ["png", "svg"])
    figure_size: tuple[int, int] = (12, 8)
    export_resolution: int = 300
    plot_unpaired_queries: bool = False
    workers: int = 1


@dataclass
class ExportCallbacks:
    """Optional UI hooks for step messages, plot progress, and warnings (no-op when None).
    """
    on_step: Callable[[str], None] | None = None
    on_plot_progress: Callable[[int, int], None] | None = None
    on_warning: Callable[[str], None] | None = None


@dataclass
class ExportResult:
    """Final outcome: the chosen folder, every saved file, and per-group plot counts.
    """
    export_root: Path
    saved_files: list[str]
    mirror_plot_count: int = 0
    summary_plot_count: int = 0


def _resolve_format_dirs(
    root_dir: Path, export_formats: list[str],
) -> tuple[dict[str, str], str]:
    """Map each requested format to its subfolder path under root_dir.
    Returns (format_dirs, html_dir) where format_dirs is {fmt: path} for
    every picture format present in export_formats and html_dir is the html
    subfolder path (empty string when html is not requested). Paths are strings only —
    the directories themselves are created later by the plot workers.
    """
    format_dirs = {fmt: str(root_dir / fmt) for fmt in _PICTURE_FORMATS if fmt in export_formats}
    html_dir = str(root_dir / "html") if "html" in export_formats else ""
    return format_dirs, html_dir


def _run_plot_tasks(
    annotated_queries: list[tuple],
    reference_store: dict,
    variant_groups: dict,
    root_dir: Path,
    plot_settings: PlotSettings,
    do_html: bool,
    saved_files: list[str],
    on_progress: Callable[[int, int], None],
    on_warning: Callable[[str], None],
) -> int:
    """Build one mirror plot per (query, reference) pair and render via a worker pool.
    Each successfully rendered task contributes its non-empty per-format output paths
    to saved_files (mutated in place). Returns the count of files written; 0 when
    no tasks were produced. Worker failures are reported through on_warning rather
    than raised, so a single bad PSM does not abort the whole export.
    """
    # Resolve format set: append html when requested, map each to a subdir
    base_formats = list(plot_settings.export_formats)
    if do_html and "html" not in base_formats:
        base_formats.append("html")
    format_dirs, html_dir = _resolve_format_dirs(root_dir, base_formats)

    # Build tasks: one per (query, reference) pair with all per-format output paths
    tasks = build_plot_tasks(
        annotated_queries, reference_store, variant_groups,
        sanitize_fn=sanitize_filename,
        format_dirs=format_dirs,
        figure_size=plot_settings.figure_size,
        export_resolution=plot_settings.export_resolution,
        html_dir=html_dir,
        plot_unpaired_queries=plot_settings.plot_unpaired_queries,
    )
    if not tasks:
        return 0

    # Dispatch: tally saved paths, surface worker errors, tick progress
    total_files = 0
    total = len(tasks)
    for completed, (_task, ok, result) in enumerate(run_tasks(
        tasks=tasks, worker=process_single_task, workers=plot_settings.workers,
    ), start=1):
        if not ok:
            on_warning(f"Plot generation failed: {result}")
        elif not result.success:
            on_warning(f"Plot generation failed: {result.error}")
        else:
            for path in (result.png_path, result.svg_path, result.pdf_path, result.html_path):
                if path:
                    saved_files.append(path)
                    total_files += 1
        on_progress(completed, total)
    return total_files


def _run_summary_plots(
    comparison_frame: pd.DataFrame | None,
    root_dir: Path,
    export_formats: list[str],
    saved_files: list[str],
) -> int:
    """Render dual-axis / ridgeline / heatmap summaries from the comparison table.
    Summaries are picture-only — html is filtered out, and an empty format set falls
    back to ["png"]. Appends every written path to saved_files and returns the
    count; returns 0 when there is no comparison frame or it is empty.
    """
    if comparison_frame is None or comparison_frame.empty:
        return 0
    summary_formats = [f for f in export_formats if f in _PICTURE_FORMATS] or ["png"]
    results = generate_all_summary_plots(
        source=comparison_frame,
        output_dir=str(root_dir),
        formats=summary_formats,
    )
    paths = [p for fmt_paths in results.values() for p in fmt_paths.values()]
    saved_files.extend(paths)
    return len(paths)


def run_export(
    *,
    export_directory: str,
    project_name: str,
    run_settings: dict,
    frames: dict[str, Any],
    annotated_queries: AnnotatedQueries,
    reference_store: dict,
    variant_groups: dict,
    comparison_frame: pd.DataFrame | None = None,
    selection: ExportSelection | None = None,
    plot_settings: PlotSettings | None = None,
    callbacks: ExportCallbacks | None = None,
) -> ExportResult:
    """Run a full export end-to-end: versioned dir, parquets, CSV/XLSX, plots, summaries.
    Keyword-only entry point shared by both frontends. Allocates (or reuses) a versioned
    folder for project_name under export_directory, writes run_settings.json,
    then runs each artefact stage gated by selection. run_settings.json is always
    counted as saved even when every stage is disabled.
    Returns an ExportResult whose saved_files lists every written path and whose
    mirror_plot_count / summary_plot_count report per-group totals. Stage failures
    that are recoverable (a bad plot worker, a missing XLSX engine) are surfaced through
    callbacks.on_warning rather than raised.
    """
    # Normalise inputs: defaults for missing dataclasses, no-op for missing callbacks
    selection = selection or ExportSelection()
    plot_settings = plot_settings or PlotSettings()
    callbacks = callbacks or ExportCallbacks()
    on_step = callbacks.on_step or (lambda *_: None)
    on_progress = callbacks.on_plot_progress or (lambda *_: None)
    on_warning = callbacks.on_warning or (lambda *_: None)

    # Stage 1 — Plan & commit folder: reuse on fingerprint match, else allocate next NNN
    root_dir = find_or_create_export_directory(export_directory, project_name, run_settings)
    commit_export_directory(root_dir, run_settings)

    saved_files: list[str] = [str(root_dir / "run_settings.json")]

    # Stage 2 — Parquet frames: write each of the four LazyFrames as {name}.parquet
    if selection.parquet:
        on_step("Saving Parquet DataFrames...")
        for name in _PARQUET_FRAMES:
            frame = frames.get(name)
            if frame is None:
                continue
            out = root_dir / f"{name}.parquet"
            frame.collect().write_parquet(out)
            saved_files.append(str(out))

    # Stage 3 — Comparison table: CSV unconditional; XLSX best-effort (needs openpyxl)
    if selection.csv and comparison_frame is not None and not comparison_frame.empty:
        on_step("Saving comparison summary...")
        csv_path = root_dir / "comparison_summary.csv"
        comparison_frame.to_csv(csv_path, index=False)
        saved_files.append(str(csv_path))
        xlsx_path = root_dir / "comparison_summary.xlsx"
        try:
            comparison_frame.to_excel(xlsx_path, index=False, engine="openpyxl")
            saved_files.append(str(xlsx_path))
        except ImportError:
            on_warning("openpyxl not installed -- skipping XLSX export")
        except Exception as error:
            on_warning(f"XLSX export failed: {error}")

    # Stage 4 — Mirror plots: one per PSM, rendered across the worker pool
    mirror_plot_count = 0
    if selection.mirror_plots or selection.html:
        on_step("Generating mirror plots...")
        # Normalise annotated_queries: dict {psm_id: (spec, df)} or list [(psm_id, spec, df)]
        if isinstance(annotated_queries, dict):
            queries_list = [(psm_id, spec, df) for psm_id, (spec, df) in annotated_queries.items()]
        else:
            queries_list = list(annotated_queries)
        mirror_plot_count = _run_plot_tasks(
            queries_list, reference_store, variant_groups, root_dir,
            plot_settings, selection.html, saved_files, on_progress, on_warning,
        )

    # Stage 5 — Summary plots: dual-axis / ridgeline / heatmap from the comparison frame
    summary_plot_count = 0
    if selection.summary_plots:
        on_step("Generating summary plots...")
        summary_plot_count = _run_summary_plots(
            comparison_frame, root_dir, plot_settings.export_formats, saved_files,
        )

    return ExportResult(
        export_root=root_dir,
        saved_files=saved_files,
        mirror_plot_count=mirror_plot_count,
        summary_plot_count=summary_plot_count,
    )
