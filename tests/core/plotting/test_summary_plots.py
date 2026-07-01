"""Unit tests for core.plotting.summary (ridgeline + multi-boxplot renderers)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from core.plotting import (
    generate_all_summary_plots,
    prepare_summary_data,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_summary_frame(n_peptides=3, n_reps=4):
    """Build a synthetic comparison_summary.csv-shaped DataFrame (3 peptides × 4 replicates)."""
    rng = np.random.RandomState(42)
    rows = []
    peptides = ["PEPTIDEK", "AAAVAAAK", "ACDEFGHIK"][:n_peptides]
    for pep in peptides:
        for i in range(n_reps):
            rows.append({
                "Peptide": pep,
                "Modification": "",
                "Charge": 2,
                "Query": f"query_{i}.mzML",
                "Reference": "ref.mzML",
                "R\u00b2 binned": float(rng.uniform(0.5, 1.0)),
                "R\u00b2 match": float(rng.uniform(0.5, 1.0)),
                "Cosine": float(rng.uniform(0.6, 1.0)),
                "Entropy": float(rng.uniform(0.6, 1.0)),
                "Dice": float(rng.uniform(0.4, 1.0)),
                "Jaccard": float(rng.uniform(0.3, 1.0)),
                "SharedRatio": float(rng.uniform(0.4, 1.0)),
            })
    return pd.DataFrame(rows)


@pytest.fixture
def summary_csv(tmp_path):
    """Write a minimal comparison_summary.csv to tmp_path and return its path."""
    frame = _make_summary_frame()
    csv_path = tmp_path / "comparison_summary.csv"
    frame.to_csv(csv_path, index=False)
    return csv_path


# ---------------------------------------------------------------------------
# TestPrepareSummaryData
# ---------------------------------------------------------------------------

class TestPrepareSummaryData:
    def test_aggregated_frame_has_expected_columns(self):
        frame = _make_summary_frame()
        _, aggregated_frame = prepare_summary_data(frame, metric="R\u00b2 binned")
        for col in ["peptidoform", "count", "median", "mean", "std"]:
            assert col in aggregated_frame.columns

    def test_peptidoform_column_built(self):
        frame = _make_summary_frame(n_peptides=2, n_reps=2)
        raw_frame, _ = prepare_summary_data(frame, metric="R\u00b2 binned")
        assert "peptidoform" in raw_frame.columns
        assert all("/" in peptidoform for peptidoform in raw_frame["peptidoform"])

    def test_invalid_metric_raises_valueerror(self):
        frame = _make_summary_frame()
        with pytest.raises(ValueError, match="not found"):
            prepare_summary_data(frame, metric="NonExistentMetric")

    def test_default_sort_is_count_ascending(self):
        """Omitting sort_by keeps the historical ascending-count ordering."""
        frame = _make_summary_frame(n_peptides=3, n_reps=4)
        _, default = prepare_summary_data(frame, metric="R² binned")
        _, count_asc = prepare_summary_data(frame, metric="R² binned", sort_by="count_asc")
        assert list(default["peptidoform"]) == list(count_asc["peptidoform"])
        assert list(count_asc["count"]) == sorted(count_asc["count"])

    def test_sort_alphabetical(self):
        frame = _make_summary_frame(n_peptides=3, n_reps=2)
        _, asc = prepare_summary_data(frame, metric="R² binned", sort_by="alpha_asc")
        _, desc = prepare_summary_data(frame, metric="R² binned", sort_by="alpha_desc")
        assert list(asc["peptidoform"]) == sorted(asc["peptidoform"])
        assert list(desc["peptidoform"]) == sorted(desc["peptidoform"], reverse=True)

    def test_sort_by_count_both_directions(self):
        """Peptidoforms with distinct replicate counts order by count each way."""
        rows = []
        for peptide, replicates in [("AAAAK", 1), ("BBBBK", 2), ("CCCCK", 3)]:
            for i in range(replicates):
                rows.append({
                    "Peptide": peptide, "Modification": "", "Charge": 2,
                    "Query": f"q{i}.mzML", "Reference": "r.mzML",
                    "R² binned": 0.5 + 0.1 * i,
                })
        frame = pd.DataFrame(rows)
        _, asc = prepare_summary_data(frame, metric="R² binned", sort_by="count_asc")
        _, desc = prepare_summary_data(frame, metric="R² binned", sort_by="count_desc")
        assert list(asc["count"]) == [1, 2, 3]
        assert list(desc["count"]) == [3, 2, 1]

    def test_sort_by_median_both_directions(self):
        frame = _make_summary_frame(n_peptides=3, n_reps=4)
        _, asc = prepare_summary_data(frame, metric="R² binned", sort_by="median_asc")
        _, desc = prepare_summary_data(frame, metric="R² binned", sort_by="median_desc")
        assert list(asc["median"]) == sorted(asc["median"])
        assert list(desc["median"]) == sorted(desc["median"], reverse=True)


# ---------------------------------------------------------------------------
# TestSummaryPlotsFromCsv
# ---------------------------------------------------------------------------

class TestSummaryPlotsFromCsv:
    def test_generates_files(self, summary_csv, tmp_path):
        results = generate_all_summary_plots(
            source=str(summary_csv),
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["ridgeline"],
        )
        assert len(results) > 0
        for plot_type, paths in results.items():
            for fmt, path in paths.items():
                assert Path(path).exists()

    def test_all_plot_types_from_frame(self, tmp_path):
        frame = _make_summary_frame(n_peptides=3, n_reps=5)
        results = generate_all_summary_plots(
            source=frame,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["ridgeline", "multi_boxplot"],
        )
        assert "ridgeline" in results
        assert "multi_boxplot" in results

    def test_multiple_metrics_write_one_file_each(self, tmp_path):
        # One boxplot per measure (issue #82): two metrics -> two distinct files.
        frame = _make_summary_frame(n_peptides=3, n_reps=5)
        results = generate_all_summary_plots(
            source=frame,
            output_dir=str(tmp_path),
            formats=["png"],
            metrics=["R² binned", "Cosine"],
            plot_types=["multi_boxplot"],
        )
        boxplots = sorted((tmp_path / "summary").glob("multi_boxplot_*.png"))
        assert {p.name for p in boxplots} == {"multi_boxplot_r2_binned.png", "multi_boxplot_cosine.png"}
        assert len(results["multi_boxplot"]) == 2

    def test_unknown_metric_skipped(self, tmp_path):
        # A metric absent from the table is logged and skipped, not raised.
        frame = _make_summary_frame()
        results = generate_all_summary_plots(
            source=frame,
            output_dir=str(tmp_path),
            formats=["png"],
            metrics=["Cosine", "NopeMetric"],
            plot_types=["ridgeline"],
        )
        assert len(results["ridgeline"]) == 1

    def test_summary_subdir_created(self, tmp_path):
        frame = _make_summary_frame()
        generate_all_summary_plots(
            source=frame,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["ridgeline"],
        )
        assert (tmp_path / "summary").exists()
