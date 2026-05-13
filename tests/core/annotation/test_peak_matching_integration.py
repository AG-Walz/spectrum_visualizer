"""
End-to-end tests for ``core.annotation.prepare_spectrum_for_plot``.

``prepare_spectrum_for_plot`` composes annotation + ion-type filtering +
normalization. These tests exercise that composition on realistic peptide
inputs.
"""
import numpy as np
import pandas as pd
import pytest

from core.annotation import prepare_spectrum_for_plot

# Known singly-charged b2 / y2 m/z for ASDDFDQASY (computed via spectrum_utils
# during the parity spike — pinned here to avoid coupling tests to PyOpenMS).
ASDDFDQASY_B2 = 159.0764
ASDDFDQASY_Y2 = 269.1132


class TestPrepareSpectrumForPlot:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.mz = np.array([ASDDFDQASY_B2, ASDDFDQASY_Y2, 500.0, 600.0])
        self.intensity = np.array([1000.0, 800.0, 200.0, 150.0])

    def test_normalization_max_is_one(self):
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "ASDDFDQASY", charge=1, normalize=True
        )
        assert df["intensity"].max() == pytest.approx(1.0, abs=1e-9)

    def test_no_normalization_unchanged(self):
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "ASDDFDQASY", charge=1, normalize=False
        )
        assert df["intensity"].max() == pytest.approx(1000.0, abs=1e-6)

    def test_empty_mz_returns_empty_df(self):
        df = prepare_spectrum_for_plot(
            np.array([]), np.array([]), "ASDDFDQASY", charge=1
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_result_has_ion_type_column(self):
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "ASDDFDQASY", charge=1
        )
        assert "ion_type" in df.columns

    def test_result_length_matches_input(self):
        mz = np.array([600.0, 200.0, 400.0])
        intensity = np.array([100.0, 200.0, 150.0])
        df = prepare_spectrum_for_plot(mz, intensity, "ASDDFDQASY", charge=1)
        assert len(df) == 3
        assert "mz" in df.columns


@pytest.mark.integration
class TestAnnotateRealSpectrum:
    """Annotate spectra from real-like data."""

    def test_annotate_asddfdqasy(self):
        mz = np.array([175.12, 290.14, 405.17, 520.20, 667.27, 782.30], dtype=np.float64)
        intensity = np.array([100, 200, 150, 300, 250, 180], dtype=np.float64)
        result = prepare_spectrum_for_plot(
            mz, intensity, "ASDDFDQASY", charge=1,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
        )
        assert len(result) == len(mz)
        assert "annotation" in result.columns
        assert result["mz"].max() <= max(mz) + 0.01

    def test_annotate_with_losses(self):
        mz = np.array([175.12, 290.14, 405.17], dtype=np.float64)
        intensity = np.array([100, 200, 150], dtype=np.float64)
        result = prepare_spectrum_for_plot(
            mz, intensity, "RTEFNLNQY", charge=2,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
            add_neutral_losses=True,
        )
        assert len(result) == len(mz)
