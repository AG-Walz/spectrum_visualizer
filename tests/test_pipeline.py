"""
Unit tests for src/pipeline.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.spectrum_loader import LoadedSpectrum
from src.peak_annotation import generate_theoretical_fragments
from src.pipeline import (
    accumulate_plot_result,
    annotate_spectrum,
    build_comparison_tasks,
    mark_skipped_tasks,
)
from src.parallel_processor import ComparisonResult, ComparisonTask


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_loaded_spectrum(
    peptidoform="PEPTIDEK",
    charge=2,
    mz=None,
    intensity=None,
    source_file="/data/query.mzML",
    is_reference=False,
    spectrum_id="scan=1",
    rt=10.0,
):
    """Build a real LoadedSpectrum with PEPTIDEK b/y m/z values."""
    if mz is None:
        frags = generate_theoretical_fragments(peptidoform, charge=charge, ion_types=["b", "y"])
        b_mzs = frags[frags["annotation"].str.match(r"^b\d\+$")]["mz"].values[:3]
        y_mzs = frags[frags["annotation"].str.match(r"^y\d\+$")]["mz"].values[:3]
        mz = np.concatenate([b_mzs, y_mzs, [500.0, 600.0]])
        intensity = np.array([1000.0, 900.0, 800.0, 750.0, 700.0, 650.0, 200.0, 150.0])

    return LoadedSpectrum(
        spectrum_id=spectrum_id,
        mz=mz,
        intensity=intensity,
        precursor_mz=500.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence="PEPTIDEK",
        source_file=source_file,
        is_reference=is_reference,
        rt=rt,
    )


# Default annotation kwargs used in annotate_spectrum calls
_ANN_KWARGS = dict(
    ion_types=["b", "y"],
    mz_tolerance=0.02,
    add_losses=False,
    add_first_prefix_ion=True,
    remove_precursor_peaks=False,
    normalize=True,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_query():
    """LoadedSpectrum for PEPTIDEK, charge=2."""
    return _make_loaded_spectrum()


@pytest.fixture
def simple_ref():
    """LoadedSpectrum for PEPTIDEK, charge=2, is_reference=True."""
    return _make_loaded_spectrum(
        source_file="/data/reference.mzML",
        is_reference=True,
        spectrum_id="scan=2",
    )


# ---------------------------------------------------------------------------
# TestAnnotateSpectrum
# ---------------------------------------------------------------------------

class TestAnnotateSpectrum:
    def test_valid_spectrum_returns_dataframe(self, simple_query):
        errors = []
        result = annotate_spectrum(simple_query, errors=errors, **_ANN_KWARGS)
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert len(errors) == 0

    def test_empty_mz_returns_none(self):
        spec = _make_loaded_spectrum(mz=np.array([]), intensity=np.array([]))
        errors = []
        result = annotate_spectrum(spec, errors=errors, **_ANN_KWARGS)
        assert result is None
        assert len(errors) == 1

    def test_zero_charge_returns_none(self):
        spec = _make_loaded_spectrum(charge=0)
        errors = []
        result = annotate_spectrum(spec, errors=errors, **_ANN_KWARGS)
        assert result is None
        assert len(errors) == 1

    def test_missing_peptidoform_returns_none(self):
        spec = _make_loaded_spectrum(peptidoform="")
        errors = []
        result = annotate_spectrum(spec, errors=errors, **_ANN_KWARGS)
        assert result is None
        assert len(errors) == 1


# ---------------------------------------------------------------------------
# TestBuildComparisonTasks
# ---------------------------------------------------------------------------

class TestBuildComparisonTasks:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.tmp_path = tmp_path
        self.png_dir = str(tmp_path / "png")
        self.svg_dir = str(tmp_path / "svg")
        self.pdf_dir = str(tmp_path / "pdf")
        self.task_counter = [0]

    def _annotate(self, spec):
        errors = []
        return annotate_spectrum(spec, errors=errors, **_ANN_KWARGS)

    def test_no_refs_unpaired_false_returns_empty(self, simple_query):
        query_df = self._annotate(simple_query)
        tasks = build_comparison_tasks(
            query=simple_query,
            ref_data=[],
            query_df=query_df,
            base_sequence="PEPTIDEK",
            modification_key="",
            png_dir=self.png_dir,
            svg_dir=self.svg_dir,
            pdf_dir=self.pdf_dir,
            export_formats=["png"],
            figsize=(10, 6),
            dpi=72,
            normalize=True,
            plot_unpaired=False,
            task_counter=self.task_counter,
        )
        assert tasks == []

    def test_no_refs_unpaired_true_returns_task(self, simple_query):
        query_df = self._annotate(simple_query)
        tasks = build_comparison_tasks(
            query=simple_query,
            ref_data=[],
            query_df=query_df,
            base_sequence="PEPTIDEK",
            modification_key="",
            png_dir=self.png_dir,
            svg_dir=self.svg_dir,
            pdf_dir=self.pdf_dir,
            export_formats=["png"],
            figsize=(10, 6),
            dpi=72,
            normalize=True,
            plot_unpaired=True,
            task_counter=self.task_counter,
        )
        assert len(tasks) == 1
        assert tasks[0].task_type == "single"

    def test_one_ref_returns_mirror_task(self, simple_query, simple_ref):
        query_df = self._annotate(simple_query)
        ref_df = self._annotate(simple_ref)
        tasks = build_comparison_tasks(
            query=simple_query,
            ref_data=[(simple_ref, ref_df)],
            query_df=query_df,
            base_sequence="PEPTIDEK",
            modification_key="",
            png_dir=self.png_dir,
            svg_dir=self.svg_dir,
            pdf_dir=self.pdf_dir,
            export_formats=["png"],
            figsize=(10, 6),
            dpi=72,
            normalize=True,
            plot_unpaired=False,
            task_counter=self.task_counter,
        )
        assert len(tasks) == 1
        assert tasks[0].task_type == "mirror"

    def test_two_refs_returns_two_tasks(self, simple_query):
        query_df = self._annotate(simple_query)
        ref1 = _make_loaded_spectrum(source_file="/data/ref1.mzML", is_reference=True)
        ref2 = _make_loaded_spectrum(source_file="/data/ref2.mzML", is_reference=True)
        ref1_df = self._annotate(ref1)
        ref2_df = self._annotate(ref2)
        tasks = build_comparison_tasks(
            query=simple_query,
            ref_data=[(ref1, ref1_df), (ref2, ref2_df)],
            query_df=query_df,
            base_sequence="PEPTIDEK",
            modification_key="",
            png_dir=self.png_dir,
            svg_dir=self.svg_dir,
            pdf_dir=self.pdf_dir,
            export_formats=["png"],
            figsize=(10, 6),
            dpi=72,
            normalize=True,
            plot_unpaired=False,
            task_counter=self.task_counter,
        )
        assert len(tasks) == 2


# ---------------------------------------------------------------------------
# TestMarkSkippedTasks
# ---------------------------------------------------------------------------

class TestMarkSkippedTasks:
    def _make_task(self, tmp_path, create_output=False):
        png_dir = str(tmp_path / "png")
        (tmp_path / "png").mkdir(parents=True, exist_ok=True)
        task = ComparisonTask(
            task_id=0,
            task_type="single",
            query_df=pd.DataFrame({"mz": [100.0], "intensity": [1.0], "annotation": ["b1"]}),
            folder_name="PEPTIDEK",
            safe_filename="test_file",
            png_dir=png_dir,
            svg_dir="",
            pdf_dir="",
            export_formats=["png"],
        )
        if create_output:
            out_dir = tmp_path / "png" / "PEPTIDEK"
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "test_file.png").write_text("fake png")
        return task

    def test_existing_file_marked_skip(self, tmp_path):
        task = self._make_task(tmp_path, create_output=True)
        assert task.skip_plot is False
        mark_skipped_tasks([task])
        assert task.skip_plot is True

    def test_missing_file_not_skipped(self, tmp_path):
        task = self._make_task(tmp_path, create_output=False)
        mark_skipped_tasks([task])
        assert task.skip_plot is False


# ---------------------------------------------------------------------------
# TestAccumulatePlotResult
# ---------------------------------------------------------------------------

class TestAccumulatePlotResult:
    def _init_summary(self):
        return {
            "plots_generated": 0,
            "files_exported": {"png": [], "svg": []},
            "errors": [],
        }

    def test_success_increments_plots_generated(self):
        results_summary = self._init_summary()
        summary_rows = []
        result = ComparisonResult(
            task_id=0,
            task_type="mirror",
            success=True,
            png_path="/some/file.png",
            summary_row={"peptide": "PEPTIDEK"},
        )
        accumulate_plot_result(result, results_summary, summary_rows)
        assert results_summary["plots_generated"] == 1

    def test_error_appended_to_errors(self):
        results_summary = self._init_summary()
        summary_rows = []
        result = ComparisonResult(
            task_id=1,
            task_type="mirror",
            success=False,
            error="Something went wrong",
        )
        accumulate_plot_result(result, results_summary, summary_rows)
        assert len(results_summary["errors"]) == 1
        assert "Something went wrong" in results_summary["errors"][0]

    def test_summary_row_added(self):
        results_summary = self._init_summary()
        summary_rows = []
        row = {"peptide": "PEPTIDEK", "charge": 2}
        result = ComparisonResult(
            task_id=0,
            task_type="mirror",
            success=True,
            summary_row=row,
        )
        accumulate_plot_result(result, results_summary, summary_rows)
        assert len(summary_rows) == 1
        assert summary_rows[0]["peptide"] == "PEPTIDEK"
