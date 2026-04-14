"""
CLI plot generation — mirror / single-spectrum worker plus parallel orchestration.

``_generate_single_plot`` is module-level so it can be submitted to
``ProcessPoolExecutor`` through ``run_tasks``. ``export_mirror_plots`` handles
fmt-dir selection + HTML routing and returns the total file count.
"""
import logging
from pathlib import Path
from typing import Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mspecviz.common.parallel import run_tasks
from mspecviz.export_dirs import sanitize_filename
from mspecviz.mirror_plot import (
    create_mirror_plot, create_plotly_mirror_plot,
    create_single_spectrum_plot, figure_to_bytes,
)
from mspecviz.query_processing import build_plot_tasks

logger = logging.getLogger("mspecviz.cli")


def _generate_single_plot(task: dict) -> int:
    """Worker: generate one mirror/single plot, write all requested formats."""
    count = 0

    if task["type"] == "mirror":
        fig = create_mirror_plot(
            upper_spectrum=task["query_df"],
            lower_spectrum=task["ref_df"],
            title=task["title"],
            upper_label=task["query_label"],
            lower_label=task["ref_label"],
            metrics_text=task["metrics_text"],
            normalize=task["normalize"],
            figsize=task.get("figsize", (12, 8)),
            dpi=task.get("dpi", 300),
        )
    else:
        fig = create_single_spectrum_plot(
            spectrum=task["query_df"],
            title=task["title"],
            label=task.get("query_label", ""),
            normalize=task["normalize"],
            figsize=task.get("figsize", (12, 8)),
            dpi=task.get("dpi", 300),
        )

    for fmt_dir_key, fmt in [("png_dir", "png"), ("svg_dir", "svg"), ("pdf_dir", "pdf")]:
        fmt_dir = task.get(fmt_dir_key, "")
        if not fmt_dir:
            continue
        out_dir = Path(fmt_dir) / task["folder_name"]
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"{task['safe_name']}.{fmt}"
        out_file.write_bytes(figure_to_bytes(fig, format=fmt))
        count += 1

    plt.close(fig)

    if task.get("html_dir") and task["type"] == "mirror":
        plotly_fig = create_plotly_mirror_plot(
            upper_spectrum=task["query_df"],
            lower_spectrum=task["ref_df"],
            title=task["title"],
            upper_label=task["query_label"],
            lower_label=task["ref_label"],
            metrics_text=task["metrics_text"],
            normalize=task["normalize"],
        )
        html_out_dir = Path(task["html_dir"]) / task["folder_name"]
        html_out_dir.mkdir(parents=True, exist_ok=True)
        plotly_fig.write_html(
            str(html_out_dir / f"{task['safe_name']}.html"),
            include_plotlyjs="cdn",
        )
        count += 1

    return count


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
        tasks=tasks,  # each task is a dict → run_tasks auto-wraps to (dict,)
        worker=_generate_single_plot,
        workers=cfg.workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
    )

    progress.finish()
    return total_files
