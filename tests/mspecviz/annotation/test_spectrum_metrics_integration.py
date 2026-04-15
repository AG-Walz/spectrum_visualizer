"""
End-to-end tests for the comparison metrics on annotated spectra.

Exercises ``compute_all_metrics`` on DataFrames produced by
``prepare_spectrum_for_plot`` rather than hand-constructed mini-frames.
Pure-data tests of individual metrics live in
``test_spectrum_metrics_unit.py``.
"""
import numpy as np
import pytest

from mspecviz.annotation import (
    compute_all_metrics,
    generate_theoretical_fragments,
    prepare_spectrum_for_plot,
)


@pytest.mark.integration
class TestMetricsWithRealAnnotation:
    """Compute metrics using annotated spectra from real data."""

    def test_self_comparison_perfect(self):
        """A spectrum compared against itself should have high similarity."""
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
