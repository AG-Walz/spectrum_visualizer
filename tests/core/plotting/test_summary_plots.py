"""Unit tests for core.plotting.summary (dual-axis + ridgeline summary renderers)."""
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


# ---------------------------------------------------------------------------
# TestSummaryPlotsFromCsv
# ---------------------------------------------------------------------------

class TestSummaryPlotsFromCsv:
    def test_generates_files(self, summary_csv, tmp_path):
        results = generate_all_summary_plots(
            source=str(summary_csv),
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["dual_axis"],
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
            plot_types=["dual_axis", "ridgeline"],
        )
        assert "dual_axis" in results
        assert "ridgeline" in results

    def test_summary_subdir_created(self, tmp_path):
        frame = _make_summary_frame()
        generate_all_summary_plots(
            source=frame,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["dual_axis"],
        )
        assert (tmp_path / "summary").exists()
