"""Shared export orchestration used by both ``cli`` and ``streamlit_app``.

Both frontends call :func:`run_export` with their own progress callbacks;
every file layout decision (directory versioning, format subdirs, parquet
names, CSV/XLSX naming, plot task dispatch, summary plots) lives here so
the two frontends produce byte-identical outputs given the same inputs.
"""
from mspecviz.export.export_pipeline import (
    ExportCallbacks,
    ExportResult,
    ExportSelection,
    PlotSettings,
    run_export,
)

__all__ = [
    "ExportCallbacks",
    "ExportResult",
    "ExportSelection",
    "PlotSettings",
    "run_export",
]