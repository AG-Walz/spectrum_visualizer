"""
Unit tests for src/summary_plots.py

MSInsight changes from MSpecViz:
- Removed: plot_dumbbell, generate_peptide_detail_plots, _shorten_labels,
           prepare_peptide_detail_data
- Added: sanitize_filename, shorten_stem (tested in test_export.py)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.summary_plots import (
    generate_all_summary_plots,
    prepare_summary_data,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_summary_df(n_peptides=3, n_reps=4):
    """Build a realistic synthetic comparison_summary DataFrame."""
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
    """Write a minimal comparison_summary.csv and return its path."""
    df = _make_summary_df()
    csv_path = tmp_path / "comparison_summary.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


# NOTE: TestSanitizeFilename and TestShortenStem are in test_export.py
# (sanitize_filename and shorten_stem are shared utilities tested there)

# ---------------------------------------------------------------------------
# TestPrepareSummaryData
# ---------------------------------------------------------------------------

class TestPrepareSummaryData:
    def test_returns_two_dataframes(self):
        df = _make_summary_df()
        raw_df, agg_df = prepare_summary_data(df, metric="R\u00b2 binned")
        assert isinstance(raw_df, pd.DataFrame)
        assert isinstance(agg_df, pd.DataFrame)

    def test_agg_df_has_expected_columns(self):
        df = _make_summary_df()
        _, agg_df = prepare_summary_data(df, metric="R\u00b2 binned")
        for col in ["peptidoform", "count", "median", "mean", "std"]:
            assert col in agg_df.columns

    def test_peptidoform_column_built(self):
        df = _make_summary_df(n_peptides=2, n_reps=2)
        raw_df, _ = prepare_summary_data(df, metric="R\u00b2 binned")
        assert "peptidoform" in raw_df.columns
        assert all("/" in pf for pf in raw_df["peptidoform"])

    def test_invalid_metric_raises_valueerror(self):
        df = _make_summary_df()
        with pytest.raises(ValueError, match="not found"):
            prepare_summary_data(df, metric="NonExistentMetric")


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

    def test_single_peptide_does_not_crash(self, tmp_path):
        # Only one peptide/peptidoform -- some plots may not have enough data
        df = _make_summary_df(n_peptides=1, n_reps=3)
        results = generate_all_summary_plots(
            source=df,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["dual_axis"],
        )
        # Should not raise -- results may be empty or have one entry
        assert isinstance(results, dict)

    def test_all_plot_types_from_df(self, tmp_path):
        df = _make_summary_df(n_peptides=3, n_reps=5)
        results = generate_all_summary_plots(
            source=df,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["dual_axis", "ridgeline"],
        )
        assert "dual_axis" in results
        assert "ridgeline" in results

    def test_empty_csv_raises_or_returns_empty(self, tmp_path):
        # CSV with required column headers but no data rows
        empty_df = pd.DataFrame(columns=["Peptide", "Modification", "Charge",
                                          "R\u00b2 binned", "R\u00b2 match", "Cosine",
                                          "Entropy", "Dice", "Jaccard", "SharedRatio"])
        csv_path = tmp_path / "empty.csv"
        empty_df.to_csv(csv_path, index=False)
        try:
            results = generate_all_summary_plots(
                source=str(csv_path),
                output_dir=str(tmp_path),
                formats=["png"],
                plot_types=["dual_axis"],
            )
            # If it doesn't raise, result must be dict (possibly empty)
            assert isinstance(results, dict)
        except Exception:
            # Acceptable to raise on empty data
            pass

    def test_summary_subdir_created(self, tmp_path):
        df = _make_summary_df()
        generate_all_summary_plots(
            source=df,
            output_dir=str(tmp_path),
            formats=["png"],
            plot_types=["dual_axis"],
        )
        assert (tmp_path / "summary").exists()
