"""Import smoke tests -- verify all critical modules load without error.

Each test is independent so a broken import in one module doesn't mask others.
"""
import importlib
import sys
from pathlib import Path

import pytest

# Ensure project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


@pytest.mark.parametrize("module_name", [
    "src.peak_annotation",
    "src.peptide_search",
    "src.spectrum_loader",
    "src.spectrum_metrics",
    "src.mirror_plot",
    "src.koina_prediction",
    "src.parallel_processor",
    "src.export",
    "src.summary_plots",
    "src.path_utils",
    "src.parquet_builder",
    "src.common.common",
    "utils.helpers",
])
def test_import(module_name):
    """Module should import without raising."""
    importlib.import_module(module_name)
