"""
End-to-end tests for the comparison metrics on annotated spectra.

Exercises ``compute_all_metrics`` on DataFrames produced by
``prepare_spectrum_for_plot`` rather than hand-constructed mini-frames.
Pure-data tests of individual metrics live in
``test_spectrum_metrics_unit.py``.
"""
import numpy as np
import pytest

from core.annotation import compute_all_metrics, prepare_spectrum_for_plot

# Singly-charged b1/b2/y1/y2/b3 m/z for ASDDFDQASY — pinned values from the
# spectrum_utils parity spike; saves the test from generating theoretical
# fragments just to learn where the peaks should land.
ASDDFDQASY_BY_MZ = np.array([72.0444, 159.0764, 182.0812, 269.1132, 274.1034])


@pytest.mark.integration
class TestMetricsWithRealAnnotation:
    """Compute metrics using annotated spectra from real data."""

    def test_self_comparison_perfect(self):
        """A spectrum compared against itself should have high similarity."""
        mz = ASDDFDQASY_BY_MZ.astype(np.float64)
        intensity = np.array([100, 200, 150, 300, 250], dtype=np.float64)

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
