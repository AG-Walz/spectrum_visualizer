"""
CLI plot generation — thin orchestration over :mod:`mspecviz.plotting`.

``build_plot_tasks`` emits :class:`ComparisonTask` instances; this module
just selects the fmt dirs, dispatches via :func:`run_tasks`, and tallies
how many files hit disk.
"""
import logging
from pathlib import Path
from typing import Dict

import matplotlib
matplotlib.use("Agg")

from mspecviz.common.parallel import run_tasks
from mspecviz.export_dirs import sanitize_filename
from mspecviz.plotting import ComparisonResult, ComparisonTask, process_single_task
from mspecviz.pipeline import build_plot_tasks

logger = logging.getLogger("mspecviz.cli")


def _count_written(result: ComparisonResult) -> int:
    """Count how many output files a worker actually wrote."""
    return sum(
        1 for p in (result.png_path, result.svg_path, result.pdf_path, result.html_path)
        if p
    )


def _generate_single_plot(task: ComparisonTask) -> int:
    """Worker: run the shared ``process_single_task`` and tally written files."""
    result = process_single_task(task)
    if not result.success:
        # run_tasks' on_error handler expects an exception to surface via the
        # worker — raise so it can be logged there.
        raise RuntimeError(result.error or "plot generation failed")
    return _count_written(result)


def export_mirror_plots(
    annotated_queries, reference_store, variant_groups,
    cfg, root_dir: Path, png_dir: Path, svg_dir: Path, pdf_dir: Path,
    normalize: bool, quiet: bool,
) -> int:
    """Build plot tasks then dispatch via ``run_tasks`` (sequential or parallel)."""
    from cli.output import TerminalProgress
    fmt_dirs: Dict[str, str] = {
        fmt: str(d)
        for fmt, d in [("png", png_dir), ("svg", svg_dir), ("pdf", pdf_dir)]
        if fmt in cfg.export_formats
    }
    html_dir = str(root_dir / "html") if "html" in cfg.export_formats else ""

    tasks = build_plot_tasks(
        annotated_queries, reference_store, variant_groups,
        normalize=normalize,
        sanitize_fn=sanitize_filename,
        fmt_dirs=fmt_dirs,
        figsize=cfg.figsize,
        dpi=cfg.dpi,
        html_dir=html_dir,
        plot_unpaired=cfg.plot_unpaired,
    )
    if not tasks:
        return 0

    progress = TerminalProgress("Generating plots", len(tasks), quiet)
    total_files = 0

    def _on_result(_task, count):
        nonlocal total_files
        total_files += count

    def _on_error(_task, exc):
        logger.error(f"Plot generation failed: {exc}")

    def _on_progress(completed, _total):
        progress.update(completed)

    run_tasks(
        tasks=tasks,
        worker=_generate_single_plot,
        workers=cfg.workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
    )

    progress.finish()
    return total_files
