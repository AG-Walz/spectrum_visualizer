"""Shared export orchestrator: directory versioning, parquets, plots, summaries.

``run_export`` is the single entry point. Both ``cli.exports`` and
``streamlit_app.views.export.pipeline`` call it with their own progress
callbacks so the on-disk layout stays identical across frontends.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import pandas as pd

from mspecviz.common.parallel import run_tasks
from mspecviz.export.directories import sanitize_filename
from mspecviz.export.directories import (
    commit_export_directory,
    find_or_create_export_directory,
)
from mspecviz.pipeline import build_plot_tasks
from mspecviz.plotting import generate_all_summary_plots, process_single_task
from mspecviz.plotting.worker_tasks import ComparisonResult

logger = logging.getLogger(__name__)


AnnotatedQueries = Union[List[tuple], Dict[str, tuple]]


@dataclass
class ExportSelection:
    """Which artifact groups the caller wants written."""
    parquet: bool = True
    csv: bool = True
    mirror_plots: bool = True
    summary_plots: bool = True
    html: bool = False


@dataclass
class PlotSettings:
    """Plot task generation parameters."""
    export_formats: List[str] = field(default_factory=lambda: ["png", "svg"])
    figure_size: Tuple[int, int] = (12, 8)
    export_resolution: int = 300
    plot_unpaired_queries: bool = False
    normalize: bool = False
    workers: int = 1


@dataclass
class ExportCallbacks:
    """UI hooks. All optional; no-op when ``None``."""
    on_step: Optional[Callable[[str], None]] = None
    on_plot_progress: Optional[Callable[[int, int], None]] = None
    on_warning: Optional[Callable[[str], None]] = None


@dataclass
class ExportResult:
    export_root: Path
    saved_files: List[str]
    mirror_plot_count: int = 0
    summary_plot_count: int = 0


def _normalize_annotated_queries(annotated_queries: AnnotatedQueries) -> List[tuple]:
    """Accept dict ``{psm_id: (spec, df)}`` or list ``[(psm_id, spec, df)]``; return list."""
    if isinstance(annotated_queries, dict):
        return [(psm_id, spec, df) for psm_id, (spec, df) in annotated_queries.items()]
    return list(annotated_queries)


def _write_parquets(
    frames: Dict[str, Any], root_dir: Path, saved_files: List[str],
) -> None:
    """Collect the four LazyFrames and write ``{name}.parquet``."""
    for name in ("psm", "peaks", "sequence", "heatmap"):
        frame = frames.get(name)
        if frame is None:
            continue
        out = root_dir / f"{name}.parquet"
        frame.collect().write_parquet(out)
        saved_files.append(str(out))


def _write_comparison_summary(
    comparison_frame: Optional[pd.DataFrame],
    root_dir: Path,
    saved_files: List[str],
    on_warning: Optional[Callable[[str], None]],
) -> None:
    """Write ``comparison_summary.csv`` and (if openpyxl is available) ``.xlsx``."""
    if comparison_frame is None or comparison_frame.empty:
        return

    csv_path = root_dir / "comparison_summary.csv"
    comparison_frame.to_csv(csv_path, index=False)
    saved_files.append(str(csv_path))

    xlsx_path = root_dir / "comparison_summary.xlsx"
    try:
        comparison_frame.to_excel(xlsx_path, index=False, engine="openpyxl")
        saved_files.append(str(xlsx_path))
    except ImportError:
        if on_warning:
            on_warning("openpyxl not installed -- skipping XLSX export")
    except Exception as e:
        if on_warning:
            on_warning(f"XLSX export failed: {e}")


def _count_written(result: ComparisonResult) -> int:
    """Tally files a single plot worker actually wrote."""
    return sum(
        1 for p in (result.png_path, result.svg_path, result.pdf_path, result.html_path)
        if p
    )


def _resolve_format_dirs(
    root_dir: Path, export_formats: List[str],
) -> Tuple[Dict[str, str], str]:
    """Return ``(fmt_dirs, html_dir)`` without creating any folders.

    ``build_plot_tasks`` uses the paths verbatim; ``process_single_task``
    creates subfolders lazily when it writes a file, so we don't leave
    empty format folders behind for plots that weren't generated.
    """
    fmt_dirs: Dict[str, str] = {}
    for fmt in ("png", "svg", "pdf"):
        if fmt in export_formats:
            fmt_dirs[fmt] = str(root_dir / fmt)
    html_dir = str(root_dir / "html") if "html" in export_formats else ""
    return fmt_dirs, html_dir


def _plot_worker(task):
    """ProcessPool-friendly worker: run ``process_single_task`` + surface errors."""
    result = process_single_task(task)
    if not result.success:
        raise RuntimeError(result.error or "plot generation failed")
    return result


def _effective_export_formats(picture_formats: List[str], do_html: bool) -> List[str]:
    """Add ``"html"`` to *picture_formats* when HTML output was requested."""
    if do_html and "html" not in picture_formats:
        return list(picture_formats) + ["html"]
    return list(picture_formats)


def _record_plot_result(
    ok: bool, result, saved_files: List[str], callbacks: ExportCallbacks,
) -> int:
    """Process one ``run_tasks`` outcome; returns count of files written by it."""
    if not ok:
        msg = f"Plot generation failed: {result}"
        logger.error(msg)
        if callbacks.on_warning:
            callbacks.on_warning(msg)
        return 0
    written = _count_written(result)
    for path in (result.png_path, result.svg_path, result.pdf_path, result.html_path):
        if path:
            saved_files.append(path)
    return written


def _run_plot_tasks(
    annotated_queries: List[tuple],
    reference_store: dict,
    variant_groups: dict,
    root_dir: Path,
    plot_settings: PlotSettings,
    do_html: bool,
    saved_files: List[str],
    callbacks: ExportCallbacks,
) -> int:
    """Build tasks via :func:`build_plot_tasks` and dispatch through :func:`run_tasks`."""
    fmt_dirs, html_dir = _resolve_format_dirs(
        root_dir, _effective_export_formats(plot_settings.export_formats, do_html),
    )
    tasks = build_plot_tasks(
        annotated_queries, reference_store, variant_groups,
        normalize=plot_settings.normalize,
        sanitize_fn=sanitize_filename,
        fmt_dirs=fmt_dirs,
        figure_size=plot_settings.figure_size,
        export_resolution=plot_settings.export_resolution,
        html_dir=html_dir,
        plot_unpaired_queries=plot_settings.plot_unpaired_queries,
    )
    if not tasks:
        return 0

    total_files = 0
    total = len(tasks)
    completed = 0
    for _task, ok, result in run_tasks(
        tasks=tasks, worker=_plot_worker, workers=plot_settings.workers,
    ):
        completed += 1
        total_files += _record_plot_result(ok, result, saved_files, callbacks)
        if callbacks.on_plot_progress:
            callbacks.on_plot_progress(completed, total)

    return total_files


def _run_summary_plots(
    comparison_frame: Optional[pd.DataFrame],
    root_dir: Path,
    export_formats: List[str],
    saved_files: List[str],
) -> int:
    """Delegate to :func:`generate_all_summary_plots` with a picture-format filter."""
    if comparison_frame is None or comparison_frame.empty:
        return 0

    summary_formats = [f for f in export_formats if f in ("png", "svg", "pdf")]
    if not summary_formats:
        summary_formats = ["png"]

    results = generate_all_summary_plots(
        source=comparison_frame,
        output_dir=str(root_dir),
        formats=summary_formats,
    )

    count = 0
    for paths_by_fmt in results.values():
        for path in paths_by_fmt.values():
            saved_files.append(path)
            count += 1
    return count


def _announce(callbacks: ExportCallbacks, message: str) -> None:
    """Push a step message through ``on_step`` if the caller registered one."""
    if callbacks.on_step:
        callbacks.on_step(message)


def _maybe_write_parquets(selection: ExportSelection, frames, root_dir: Path,
                          saved_files: List[str], callbacks: ExportCallbacks) -> None:
    if not selection.parquet:
        return
    _announce(callbacks, "Saving Parquet DataFrames...")
    _write_parquets(frames, root_dir, saved_files)


def _maybe_write_summary_csv(selection: ExportSelection, comparison_frame,
                             root_dir: Path, saved_files: List[str],
                             callbacks: ExportCallbacks) -> None:
    if not selection.csv:
        return
    _announce(callbacks, "Saving comparison summary...")
    _write_comparison_summary(comparison_frame, root_dir, saved_files, callbacks.on_warning)


def _maybe_run_mirror_plots(
    selection: ExportSelection, annotated_queries, reference_store, variant_groups,
    root_dir: Path, plot_settings: PlotSettings, saved_files: List[str],
    callbacks: ExportCallbacks,
) -> int:
    if not (selection.mirror_plots or selection.html):
        return 0
    _announce(callbacks, "Generating mirror plots...")
    return _run_plot_tasks(
        _normalize_annotated_queries(annotated_queries),
        reference_store, variant_groups, root_dir,
        plot_settings, selection.html, saved_files, callbacks,
    )


def _maybe_run_summary_plots(
    selection: ExportSelection, comparison_frame, root_dir: Path,
    plot_settings: PlotSettings, saved_files: List[str], callbacks: ExportCallbacks,
) -> int:
    if not selection.summary_plots:
        return 0
    _announce(callbacks, "Generating summary plots...")
    return _run_summary_plots(
        comparison_frame, root_dir, plot_settings.export_formats, saved_files,
    )


def run_export(
    *,
    export_directory: str,
    project_name: str,
    run_settings: dict,
    frames: Dict[str, Any],
    annotated_queries: AnnotatedQueries,
    reference_store: dict,
    variant_groups: dict,
    comparison_frame: Optional[pd.DataFrame] = None,
    selection: Optional[ExportSelection] = None,
    plot_settings: Optional[PlotSettings] = None,
    callbacks: Optional[ExportCallbacks] = None,
) -> ExportResult:
    """Orchestrate a full export: versioned dir, parquets, CSV/XLSX, plots, summaries.

    Plans (or reuses) ``{project}_{NNN}_{YYYYMMDD}/`` via
    :func:`find_or_create_export_directory`, writes ``run_settings.json``, then
    dispatches each optional group in ``selection``. Progress and warnings are
    surfaced through ``callbacks`` so the CLI can render a terminal bar and
    Streamlit can render ``st.caption`` / ``st.warning``.
    """
    selection = selection or ExportSelection()
    plot_settings = plot_settings or PlotSettings()
    callbacks = callbacks or ExportCallbacks()

    root_dir, _png, _svg, _pdf = find_or_create_export_directory(
        export_directory, project_name, run_settings,
    )
    commit_export_directory(root_dir, run_settings)

    saved_files: List[str] = [str(root_dir / "run_settings.json")]
    _maybe_write_parquets(selection, frames, root_dir, saved_files, callbacks)
    _maybe_write_summary_csv(
        selection, comparison_frame, root_dir, saved_files, callbacks,
    )
    mirror_plot_count = _maybe_run_mirror_plots(
        selection, annotated_queries, reference_store, variant_groups,
        root_dir, plot_settings, saved_files, callbacks,
    )
    summary_plot_count = _maybe_run_summary_plots(
        selection, comparison_frame, root_dir, plot_settings, saved_files, callbacks,
    )

    return ExportResult(
        export_root=root_dir,
        saved_files=saved_files,
        mirror_plot_count=mirror_plot_count,
        summary_plot_count=summary_plot_count,
    )
