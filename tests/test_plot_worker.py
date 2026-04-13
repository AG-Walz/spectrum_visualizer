"""
Unit tests for src/plot_worker.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.plot_worker import ComparisonTask, process_single_task
from src.peak_annotation import prepare_spectrum_for_plot, generate_theoretical_fragments


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_annotated_df(peptidoform="PEPTIDEK", charge=2):
    """Build a real annotated DataFrame using prepare_spectrum_for_plot."""
    frags = generate_theoretical_fragments(peptidoform, charge=charge, ion_types=["b", "y"])
    b2_mz = float(frags[frags["annotation"] == "b2+"]["mz"].values[0])
    y2_mz = float(frags[frags["annotation"] == "y2+"]["mz"].values[0])
    mz = np.array([b2_mz, y2_mz, 500.0, 600.0])
    intensity = np.array([1000.0, 800.0, 200.0, 150.0])
    return prepare_spectrum_for_plot(mz, intensity, peptidoform, charge=charge)


def _make_single_task(tmp_path, skip_plot=False):
    """Create a single-spectrum ComparisonTask."""
    png_dir = str(tmp_path / "png")
    (tmp_path / "png").mkdir(parents=True, exist_ok=True)
    return ComparisonTask(
        task_id=1,
        task_type="single",
        query_df=_make_annotated_df(),
        query_label="Test Query",
        plot_title="PEPTIDEK",
        folder_name="PEPTIDEK",
        safe_filename="test_single",
        png_dir=png_dir,
        svg_dir="",
        pdf_dir="",
        export_formats=["png"],
        figsize=(8, 5),
        dpi=72,
        normalize=True,
        skip_plot=skip_plot,
    )


def _make_mirror_task(tmp_path, skip_plot=False):
    """Create a mirror ComparisonTask with both query and reference DFs."""
    png_dir = str(tmp_path / "png")
    (tmp_path / "png").mkdir(parents=True, exist_ok=True)
    return ComparisonTask(
        task_id=2,
        task_type="mirror",
        query_df=_make_annotated_df(),
        ref_df=_make_annotated_df(),
        query_label="Query",
        ref_label="Reference",
        plot_title="PEPTIDEK",
        folder_name="PEPTIDEK",
        safe_filename="test_mirror",
        png_dir=png_dir,
        svg_dir="",
        pdf_dir="",
        export_formats=["png"],
        figsize=(8, 5),
        dpi=72,
        normalize=True,
        skip_plot=skip_plot,
        query_file="query.mzML",
        ref_file="ref.mzML",
        peptide="PEPTIDEK",
        modification="",
        charge=2,
    )


# ---------------------------------------------------------------------------
# TestOutputsExist
# ---------------------------------------------------------------------------

class TestOutputsExist:
    def test_all_formats_present_returns_true(self, tmp_path):
        task = ComparisonTask(
            task_id=0,
            task_type="single",
            query_df=pd.DataFrame(),
            folder_name="PEPTIDEK",
            safe_filename="test_file",
            png_dir=str(tmp_path / "png"),
            svg_dir=str(tmp_path / "svg"),
            pdf_dir="",
            export_formats=["png", "svg"],
        )
        # Create both expected output files
        (tmp_path / "png" / "PEPTIDEK").mkdir(parents=True)
        (tmp_path / "svg" / "PEPTIDEK").mkdir(parents=True)
        (tmp_path / "png" / "PEPTIDEK" / "test_file.png").write_text("fake")
        (tmp_path / "svg" / "PEPTIDEK" / "test_file.svg").write_text("fake")
        assert task.outputs_exist() is True

    def test_missing_format_returns_false(self, tmp_path):
        task = ComparisonTask(
            task_id=0,
            task_type="single",
            query_df=pd.DataFrame(),
            folder_name="PEPTIDEK",
            safe_filename="test_file",
            png_dir=str(tmp_path / "png"),
            svg_dir=str(tmp_path / "svg"),
            pdf_dir="",
            export_formats=["png", "svg"],
        )
        # Only create the PNG, not the SVG
        (tmp_path / "png" / "PEPTIDEK").mkdir(parents=True)
        (tmp_path / "png" / "PEPTIDEK" / "test_file.png").write_text("fake")
        assert task.outputs_exist() is False

    def test_empty_formats_returns_true(self, tmp_path):
        task = ComparisonTask(
            task_id=0,
            task_type="single",
            query_df=pd.DataFrame(),
            folder_name="PEPTIDEK",
            safe_filename="test_file",
            png_dir="",
            svg_dir="",
            pdf_dir="",
            export_formats=[],
        )
        assert task.outputs_exist() is True


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
