"""Import smoke tests — verify all critical modules load without error.

Each test is independent so a broken import in one module doesn't mask others.
"""
import importlib
import sys
from pathlib import Path

import pytest

# Ensure project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


@pytest.mark.parametrize("module_name", [
    "src.spectrum_metrics",
    "src.spectrum_loader",
    "src.parquet_cache",
    "src.peak_annotation",
    "src.mirror_plot",
    "src.export",
    "src.peptide_search",
    "src.path_utils",
    "src.parallel_processor",
    "src.summary_plots",
    "src.common.common",
    "cli",
])
def test_import(module_name):
    """Module should import without raising."""
    importlib.import_module(module_name)
