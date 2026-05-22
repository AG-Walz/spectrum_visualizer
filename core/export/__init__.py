"""Shared export orchestration (called by cli.main and streamlit_app.views.export).
On-disk layout decisions live here so both frontends produce byte-identical outputs — see run_export.
"""
# The export pipeline itself — run_export plus the request/result/callback
# dataclasses that frontends populate to drive a single export run.
from core.export.export_pipeline import (
    ExportCallbacks,
    ExportResult,
    ExportSelection,
    PlotSettings,
    run_export,
)
# Reproducibility metadata — the run_settings.json schema and its builder,
# written alongside every export so a past run can be inspected or reloaded.
from core.export.run_settings import RUN_SETTINGS_FIELDS, build_run_settings

__all__ = [
    "ExportCallbacks",
    "ExportResult",
    "ExportSelection",
    "PlotSettings",
    "RUN_SETTINGS_FIELDS",
    "build_run_settings",
    "run_export",
]