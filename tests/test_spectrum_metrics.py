"""
Unit tests for mspecviz/spectrum_metrics.py
"""
import math

import numpy as np
import pandas as pd
import pytest

from mspecviz.spectrum_metrics import (
    compute_all_metrics,
    compute_cosine_similarity,
    compute_r2_binned,
    compute_r2_matching,
    compute_shared_peak_dice,
    compute_shared_peak_jaccard,
    compute_shared_peak_ratio,
    compute_spectral_entropy_similarity,
    format_metrics_text,
    shorten_filenames,
)


# ---------------------------------------------------------------------------
# Helpers to build annotated DataFrames
# ---------------------------------------------------------------------------

def _make_df(mz, intensity, annotations, is_annotated=None):
    """Build a DataFrame mimicking prepare_spectrum_for_plot output."""
    n = len(mz)
    if is_annotated is None:
        is_annotated = [a != "unannotated" for a in annotations]
    return pd.DataFrame({
        "mz": mz,
        "intensity": intensity,
        "annotation": annotations,
        "is_annotated": is_annotated,
    })


def _identical_pair():
    """Two identical annotated spectra."""
    mz = [100.0, 200.0, 300.0, 400.0, 500.0]
    intensity = [0.5, 1.0, 0.3, 0.8, 0.6]
    annotations = ["b1", "b2", "y1", "y2", "y3"]
    df = _make_df(mz, intensity, annotations)
    return df, df.copy()


def _disjoint_pair():
    """Two spectra with no overlapping annotations."""
    q = _make_df(
        [100.0, 200.0],
        [1.0, 0.5],
        ["b1", "b2"],
    )
    r = _make_df(
        [300.0, 400.0],
        [0.8, 0.4],
        ["y1", "y2"],
    )
    return q, r


def _partial_overlap_pair():
    """Two spectra sharing some annotations."""
    q = _make_df(
        [100.0, 200.0, 300.0],
        [1.0, 0.5, 0.3],
        ["b1", "b2", "unannotated"],
    )
    r = _make_df(
        [100.0, 200.0, 400.0],
        [0.9, 0.6, 0.4],
        ["b1", "y1", "unannotated"],
    )
    return q, r


# ---------------------------------------------------------------------------
# shorten_filenames
# ---------------------------------------------------------------------------

class TestShortenFilenames:
    @pytest.mark.parametrize("filenames, check", [
        ([], lambda r: r == []),
        (["Experiment_batch1_sample_A.mzML"],
         lambda r: r == ["Experiment_batch1_sample_A"]),
        (["alpha.mzML", "beta.mzML"],
         lambda r: r == ["alpha", "beta"]),
        (["same_name.mzML", "same_name.mzML"],
         lambda r: r == ["same_name", "same_name"]),
    ], ids=["empty", "single_file", "no_common_parts", "all_common_fallback"])
    def test_basic_cases(self, filenames, check):
        assert check(shorten_filenames(filenames))

    def test_common_prefix_suffix_produces_distinct_results(self):
        filenames = [
            "Experiment_batch1_sample_A.mzML",
            "Experiment_batch1_sample_B.mzML",
        ]
        result = shorten_filenames(filenames)
        assert len(result) == 2
        assert result[0] != result[1]


# ---------------------------------------------------------------------------
# Parametrized tests for identical/disjoint behavior (shared across 5 metrics)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("metric_fn, identical_val, disjoint_val", [
    (compute_r2_matching, 1.0, float("nan")),
    (compute_cosine_similarity, 1.0, float("nan")),
    (compute_shared_peak_dice, 1.0, 0.0),
    (compute_shared_peak_jaccard, 1.0, 0.0),
    (compute_shared_peak_ratio, 1.0, 0.0),
])
class TestMetricBehavior:
    def test_identical_returns_perfect(self, metric_fn, identical_val, disjoint_val):
        q, r = _identical_pair()
        val = metric_fn(q, r)
        if math.isnan(identical_val):
            assert math.isnan(val)
        else:
            assert val == pytest.approx(identical_val, abs=1e-6)

    def test_disjoint_returns_expected(self, metric_fn, identical_val, disjoint_val):
        q, r = _disjoint_pair()
        val = metric_fn(q, r)
        if math.isnan(disjoint_val):
            assert math.isnan(val)
        else:
            assert val == pytest.approx(disjoint_val, abs=1e-6)


# ---------------------------------------------------------------------------
# Edge cases unique to individual metrics
# ---------------------------------------------------------------------------

class TestR2Matching:
    def test_partial_overlap(self):
        q, r = _partial_overlap_pair()
        val = compute_r2_matching(q, r)
        # Only "b1" is shared -- need >=2 points for regression
        assert math.isnan(val)

    def test_good_correlation(self):
        q = _make_df([100, 200, 300], [0.2, 0.5, 1.0], ["b1", "b2", "b3"])
        r = _make_df([100, 200, 300], [0.25, 0.55, 0.95], ["b1", "b2", "b3"])
        val = compute_r2_matching(q, r)
        assert val > 0.9


class TestR2Binned:
    def test_identical_spectra(self):
        q, r = _identical_pair()
        val = compute_r2_binned(q, r)
        assert val == pytest.approx(1.0, abs=1e-6)

    def test_empty_spectrum(self):
        q = _make_df([], [], [])
        r = _make_df([100], [1.0], ["b1"])
        val = compute_r2_binned(q, r)
        assert math.isnan(val)


class TestSpectralEntropy:
    def test_identical(self):
        q, r = _identical_pair()
        val = compute_spectral_entropy_similarity(q, r)
        if not math.isnan(val):  # skip if ms_entropy not installed
            assert val > 0.9

    def test_empty(self):
        q = _make_df([], [], [])
        r = _make_df([100], [1.0], ["b1"])
        val = compute_spectral_entropy_similarity(q, r)
        assert math.isnan(val)


class TestSharedPeakDice:
    def test_no_annotations(self):
        q = _make_df([100], [1.0], ["unannotated"])
        r = _make_df([100], [1.0], ["unannotated"])
        val = compute_shared_peak_dice(q, r)
        assert math.isnan(val)


# ---------------------------------------------------------------------------
# compute_all_metrics
# ---------------------------------------------------------------------------

class TestComputeAllMetrics:
    def test_returns_all_keys(self):
        q, r = _identical_pair()
        m = compute_all_metrics(q, r)
        expected_keys = {
            "R2m", "R2b", "Cos", "Ent", "Dice", "Jacc", "ShR",
            "total_q", "total_r", "ann_q", "ann_r", "shared",
        }
        assert set(m.keys()) == expected_keys

    def test_counts_correct(self):
        q, r = _identical_pair()
        m = compute_all_metrics(q, r)
        assert m["total_q"] == 5
        assert m["ann_q"] == 5
        assert m["shared"] == 5

    def test_partial_counts(self):
        q, r = _partial_overlap_pair()
        m = compute_all_metrics(q, r)
        assert m["ann_q"] == 2   # b1, b2
        assert m["ann_r"] == 2   # b1, y1
        assert m["shared"] == 1  # b1


# ---------------------------------------------------------------------------
# format_metrics_text
# ---------------------------------------------------------------------------

class TestFormatMetricsText:
    @pytest.fixture
    def all_valid_metrics(self):
        return {"R2m": 0.85, "R2b": 0.72, "Cos": 0.91, "Ent": 0.88,
                "Dice": 0.65, "Jacc": 0.5, "ShR": 0.7}

    def test_format_with_valid_values(self, all_valid_metrics):
        text = format_metrics_text(all_valid_metrics)
        assert "Match" in text and "Binned" in text
        assert "0.850" in text and "0.720" in text
        # Should produce 7 lines (one per metric)
        assert len(text.strip().split("\n")) == 7

    def test_nan_shows_na(self):
        metrics = {"R2m": float("nan"), "R2b": 0.72, "Cos": float("nan"),
                   "Ent": 0.88, "Dice": float("nan"), "Jacc": 0.5, "ShR": 0.7}
        assert "N/A" in format_metrics_text(metrics)


# ---------------------------------------------------------------------------
# Integration tests with annotated spectra
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestMetricsWithRealAnnotation:
    """Compute metrics using annotated spectra from real data."""

    def test_self_comparison_perfect(self):
        """A spectrum compared against itself should have high similarity."""
        from mspecviz.peak_annotation import prepare_spectrum_for_plot, generate_theoretical_fragments

        # Use actual theoretical fragment m/z values so peaks get annotated
        frags = generate_theoretical_fragments("ASDDFDQASY", charge=1, ion_types=["b", "y"])
        mz = frags["mz"].values[:5].astype(np.float64)
        intensity = np.array([100, 200, 150, 300, 250], dtype=np.float64)[:len(mz)]

        df = prepare_spectrum_for_plot(
            mz, intensity, "ASDDFDQASY", charge=1,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
        )

        # Verify we got annotated peaks
        annotated_count = df["is_annotated"].sum()
        assert annotated_count > 0, "No peaks were annotated -- test data invalid"

        metrics = compute_all_metrics(df, df)
        # Self-comparison should yield perfect or near-perfect scores
        assert metrics["Cos"] == 1.0 or abs(metrics["Cos"] - 1.0) < 0.001
        assert metrics["Dice"] == 1.0 or abs(metrics["Dice"] - 1.0) < 0.001
