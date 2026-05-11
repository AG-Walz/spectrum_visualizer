"""
Unit tests for the seven comparison metrics in ``mspecviz.annotation``.

Pure-data tests built from hand-constructed DataFrames — no theoretical
fragment generation, no real spectra. End-to-end checks against annotated
spectra live in ``test_spectrum_metrics_integration.py``.
"""
import math

import pandas as pd
import pytest

from mspecviz.annotation import (
    compute_all_metrics,
    compute_cosine_similarity,
    compute_r2_binned,
    compute_r2_matching,
    compute_spectral_entropy_similarity,
    format_metrics_text,
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
    query = _make_df([100.0, 200.0], [1.0, 0.5], ["b1", "b2"])
    reference = _make_df([300.0, 400.0], [0.8, 0.4], ["y1", "y2"])
    return query, reference


def _partial_overlap_pair():
    """Two spectra sharing some annotations."""
    query = _make_df([100.0, 200.0, 300.0], [1.0, 0.5, 0.3], ["b1", "b2", "unannotated"])
    reference = _make_df([100.0, 200.0, 400.0], [0.9, 0.6, 0.4], ["b1", "y1", "unannotated"])
    return query, reference


@pytest.mark.parametrize("metric_fn, identical_value, disjoint_value", [
    (compute_r2_matching, 1.0, float("nan")),
    (compute_cosine_similarity, 1.0, float("nan")),
])
class TestMetricBehavior:
    def test_identical_returns_perfect(self, metric_fn, identical_value, disjoint_value):
        query, reference = _identical_pair()
        result = metric_fn(query, reference)
        if math.isnan(identical_value):
            assert math.isnan(result)
        else:
            assert result == pytest.approx(identical_value, abs=1e-6)

    def test_disjoint_returns_expected(self, metric_fn, identical_value, disjoint_value):
        query, reference = _disjoint_pair()
        result = metric_fn(query, reference)
        if math.isnan(disjoint_value):
            assert math.isnan(result)
        else:
            assert result == pytest.approx(disjoint_value, abs=1e-6)


class TestR2Matching:
    def test_partial_overlap(self):
        query, reference = _partial_overlap_pair()
        result = compute_r2_matching(query, reference)
        # Only "b1" is shared -- need >=2 points for regression
        assert math.isnan(result)

    def test_good_correlation(self):
        query = _make_df([100, 200, 300], [0.2, 0.5, 1.0], ["b1", "b2", "b3"])
        reference = _make_df([100, 200, 300], [0.25, 0.55, 0.95], ["b1", "b2", "b3"])
        result = compute_r2_matching(query, reference)
        assert result > 0.9


class TestR2Binned:
    def test_identical_spectra(self):
        query, reference = _identical_pair()
        result = compute_r2_binned(query, reference)
        assert result == pytest.approx(1.0, abs=1e-6)

    def test_empty_spectrum(self):
        query = _make_df([], [], [])
        reference = _make_df([100], [1.0], ["b1"])
        result = compute_r2_binned(query, reference)
        assert math.isnan(result)


class TestSpectralEntropy:
    def test_identical(self):
        query, reference = _identical_pair()
        result = compute_spectral_entropy_similarity(query, reference)
        if not math.isnan(result):  # skip if ms_entropy not installed
            assert result > 0.9

    def test_empty(self):
        query = _make_df([], [], [])
        reference = _make_df([100], [1.0], ["b1"])
        result = compute_spectral_entropy_similarity(query, reference)
        assert math.isnan(result)


class TestSetMetricsViaComputeAll:
    """Set metrics (Dice/Jaccard/ShR) are exposed via ``compute_all_metrics`` only."""

    def test_identical_returns_perfect(self):
        query, reference = _identical_pair()
        metrics = compute_all_metrics(query, reference)
        assert metrics["Dice"] == pytest.approx(1.0, abs=1e-6)
        assert metrics["Jacc"] == pytest.approx(1.0, abs=1e-6)
        assert metrics["ShR"] == pytest.approx(1.0, abs=1e-6)

    def test_disjoint_returns_zero(self):
        query, reference = _disjoint_pair()
        metrics = compute_all_metrics(query, reference)
        assert metrics["Dice"] == 0.0
        assert metrics["Jacc"] == 0.0
        assert metrics["ShR"] == 0.0

    def test_no_annotations(self):
        query = _make_df([100], [1.0], ["unannotated"])
        reference = _make_df([100], [1.0], ["unannotated"])
        metrics = compute_all_metrics(query, reference)
        assert math.isnan(metrics["Dice"])


class TestComputeAllMetrics:
    def test_returns_all_keys(self):
        query, reference = _identical_pair()
        metrics = compute_all_metrics(query, reference)
        expected_keys = {
            "R2m", "R2b", "Cos", "Ent", "Dice", "Jacc", "ShR",
            "query_peak_count", "reference_peak_count",
            "query_annotation_count", "reference_annotation_count", "shared",
        }
        assert set(metrics.keys()) == expected_keys

    def test_counts_correct(self):
        query, reference = _identical_pair()
        metrics = compute_all_metrics(query, reference)
        assert metrics["query_peak_count"] == 5
        assert metrics["query_annotation_count"] == 5
        assert metrics["shared"] == 5

    def test_partial_counts(self):
        query, reference = _partial_overlap_pair()
        metrics = compute_all_metrics(query, reference)
        assert metrics["query_annotation_count"] == 2     # b1, b2
        assert metrics["reference_annotation_count"] == 2 # b1, y1
        assert metrics["shared"] == 1                     # b1


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
