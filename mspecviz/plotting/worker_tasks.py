"""Task / result dataclasses for parallel plot generation.

Kept separate from ``worker_execution`` so they remain lightweight to import
(no matplotlib / plotly dependency) for callers that only need to construct
tasks — e.g. the CLI builder in ``mspecviz.pipeline``.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class ComparisonTask:
    """All data needed to produce one plot in a worker process."""

    task_id: int
    task_type: str  # "mirror" or "single"

    # Annotated spectrum DataFrames (picklable)
    query_df: pd.DataFrame
    ref_df: Optional[pd.DataFrame] = None

    # Labels for the plot
    query_label: str = ""
    ref_label: str = ""
    plot_title: str = ""

    # Export paths (strings -- workers create subdirs as needed)
    folder_name: str = ""
    safe_filename: str = ""
    png_dir: str = ""
    svg_dir: str = ""
    pdf_dir: str = ""

    # Plot/export parameters
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    figure_size: Tuple[int, int] = (12, 8)
    export_dpi: int = 300
    normalize: bool = True

    # Output path for interactive Plotly HTML (empty string = skip HTML output)
    html_plot_path: str = ""

    # Skip plot generation (metrics only -- plot already exists on disk)
    skip_plot: bool = False

    # Pre-formatted metrics box text. When non-empty, the worker uses this
    # directly instead of recomputing via ``format_metrics_text`` — lets the
    # caller (e.g. the CLI) compute metrics once in the main process and
    # avoid re-doing the work per worker.
    metrics_text: str = ""

    # Summary metadata (for CSV row)
    query_file: str = ""
    ref_file: str = ""
    peptide: str = ""
    modification: str = ""
    charge: int = 0

    def outputs_exist(self) -> bool:
        """Return True if every requested output format already exists on disk."""
        fmt_dirs = {"png": self.png_dir, "svg": self.svg_dir, "pdf": self.pdf_dir}
        for fmt, fmt_dir in fmt_dirs.items():
            if fmt in self.export_formats and fmt_dir:
                out = Path(fmt_dir) / self.folder_name / f"{self.safe_filename}.{fmt}"
                if not out.exists():
                    return False
        if "html" in self.export_formats and self.html_plot_path:
            if not Path(self.html_plot_path).exists():
                return False
        return True


@dataclass
class ComparisonResult:
    """Output from one worker invocation."""

    task_id: int
    task_type: str
    success: bool
    error: str = ""

    png_path: str = ""
    svg_path: str = ""
    pdf_path: str = ""
    html_path: str = ""

    # Dict with metrics for summary CSV (only for mirror tasks)
    summary_row: Optional[Dict] = None
