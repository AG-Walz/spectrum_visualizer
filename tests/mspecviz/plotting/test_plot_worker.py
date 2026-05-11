"""
Unit tests for mspecviz/plot_worker.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mspecviz.plotting import ComparisonTask, process_single_task
from mspecviz.annotation import prepare_spectrum_for_plot

# Pinned singly-charged b2 / y2 m/z values for PEPTIDEK (computed via
# spectrum_utils.MsmsSpectrum.annotate_proforma) — keeps the test from
# coupling to whichever theoretical-generation library we use today.
PEPTIDEK_B2_PLUS = 227.1026
PEPTIDEK_Y2_PLUS = 276.1554


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_annotated_frame(peptidoform="PEPTIDEK", charge=2):
    """Build a real annotated DataFrame using prepare_spectrum_for_plot."""
    mz = np.array([PEPTIDEK_B2_PLUS, PEPTIDEK_Y2_PLUS, 500.0, 600.0])
    intensity = np.array([1000.0, 800.0, 200.0, 150.0])
    return prepare_spectrum_for_plot(mz, intensity, peptidoform, charge=charge)


def _make_single_task(tmp_path, skip_plot=False):
    """Create a single-spectrum ComparisonTask."""
    png_dir = str(tmp_path / "png")
    (tmp_path / "png").mkdir(parents=True, exist_ok=True)
    return ComparisonTask(
        task_id=1,
        task_type="single",
        query_frame=_make_annotated_frame(),
        query_label="Test Query",
        plot_title="PEPTIDEK",
        folder_name="PEPTIDEK",
        safe_filename="test_single",
        png_dir=png_dir,
        svg_dir="",
        pdf_dir="",
        export_formats=["png"],
        figure_size=(8, 5),
        export_resolution=72,
        skip_plot=skip_plot,
    )


def _make_mirror_task(tmp_path, skip_plot=False):
    """Create a mirror ComparisonTask with both query and reference frames."""
    png_dir = str(tmp_path / "png")
    (tmp_path / "png").mkdir(parents=True, exist_ok=True)
    return ComparisonTask(
        task_id=2,
        task_type="mirror",
        query_frame=_make_annotated_frame(),
        reference_frame=_make_annotated_frame(),
        query_label="Query",
        reference_label="Reference",
        plot_title="PEPTIDEK",
        folder_name="PEPTIDEK",
        safe_filename="test_mirror",
        png_dir=png_dir,
        svg_dir="",
        pdf_dir="",
        export_formats=["png"],
        figure_size=(8, 5),
        export_resolution=72,
        skip_plot=skip_plot,
        query_file="query.mzML",
        reference_file="ref.mzML",
        peptide="PEPTIDEK",
        modification="",
        charge=2,
    )


# ---------------------------------------------------------------------------
# TestProcessSingleTask
# ---------------------------------------------------------------------------

class TestProcessSingleTask:
    def test_mirror_task_produces_png(self, tmp_path):
        task = _make_mirror_task(tmp_path)
        result = process_single_task(task)
        assert result.success is True
        assert result.png_path != ""
        assert Path(result.png_path).exists()

    def test_single_spectrum_task_produces_png(self, tmp_path):
        task = _make_single_task(tmp_path)
        result = process_single_task(task)
        assert result.success is True
        assert result.png_path != ""
        assert Path(result.png_path).exists()

    def test_skip_plot_does_not_write_file(self, tmp_path):
        task = _make_mirror_task(tmp_path, skip_plot=True)
        result = process_single_task(task)
        assert result.success is True
        # skip_plot=True means no file written
        assert result.png_path == ""
        # But metrics should still be computed for mirror tasks
        assert result.summary_row is not None

    def test_result_has_expected_fields(self, tmp_path):
        task = _make_mirror_task(tmp_path)
        result = process_single_task(task)
        assert hasattr(result, "task_id")
        assert hasattr(result, "task_type")
        assert hasattr(result, "success")
        assert hasattr(result, "png_path")
        assert hasattr(result, "summary_row")
