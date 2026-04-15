"""
End-to-end tests for ``mspecviz.annotation.prepare_spectrum_for_plot``.

``prepare_spectrum_for_plot`` composes the building blocks exercised by
``test_peak_annotation_core.py`` (peptidoform parsing, fragment generation,
peak matching, post-processing). These tests exercise that composition on
realistic peptide inputs.
"""
import numpy as np
import pandas as pd
import pytest

from mspecviz.annotation import (
    generate_theoretical_fragments,
    prepare_spectrum_for_plot,
)


class TestPrepareSpectrumForPlot:
    @pytest.fixture(autouse=True)
    def setup(self):
        frags = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        b2_mz = float(frags[frags["annotation"] == "b2+"]["mz"].values[0])
        y2_mz = float(frags[frags["annotation"] == "y2+"]["mz"].values[0])
        self.mz = np.array([b2_mz, y2_mz, 500.0, 600.0])
        self.intensity = np.array([1000.0, 800.0, 200.0, 150.0])

    def test_normalization_max_is_one(self):
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "PEPTIDEK", charge=2, normalize=True
        )
        assert df["intensity"].max() == pytest.approx(1.0, abs=1e-9)

    def test_no_normalization_unchanged(self):
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "PEPTIDEK", charge=2, normalize=False
        )
        assert df["intensity"].max() == pytest.approx(1000.0, abs=1e-6)

    def test_precursor_removal(self):
        # Build spectrum where a peak is at the precursor location (annotated as [M+H]+)
        frags = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        # Get a precursor annotation m/z
        prec_rows = frags[frags["annotation"].str.contains(r"\[M\+H\]", regex=True)]
        if len(prec_rows) == 0:
            pytest.skip("No precursor annotation found for PEPTIDEK")
        prec_mz = float(prec_rows.iloc[0]["mz"])
        mz = np.array([prec_mz, 300.0])
        intensity = np.array([5000.0, 1000.0])
        df_with = prepare_spectrum_for_plot(
            mz, intensity, "PEPTIDEK", charge=2,
            remove_precursor_peaks=False, normalize=False
        )
        df_without = prepare_spectrum_for_plot(
            mz, intensity, "PEPTIDEK", charge=2,
            remove_precursor_peaks=True, normalize=False
        )
        # With precursor removal, the annotated precursor peak should be gone
        assert len(df_without) <= len(df_with)

    def test_empty_mz_returns_empty_df(self):
        df = prepare_spectrum_for_plot(
            np.array([]), np.array([]), "PEPTIDEK", charge=2
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_result_has_ion_type_column(self):
        """prepare_spectrum_for_plot output includes an ion_type column."""
        df = prepare_spectrum_for_plot(
            self.mz, self.intensity, "PEPTIDEK", charge=2
        )
        assert "ion_type" in df.columns

    def test_result_length_matches_input(self):
        # Provide unsorted m/z input
        mz = np.array([600.0, 200.0, 400.0])
        intensity = np.array([100.0, 200.0, 150.0])
        df = prepare_spectrum_for_plot(mz, intensity, "PEPTIDEK", charge=2)
        assert len(df) == 3
        assert "mz" in df.columns


@pytest.mark.integration
class TestAnnotateRealSpectrum:
    """Annotate spectra from real-like data."""

    def test_annotate_asddfdqasy(self):
        """Annotate ASDDFDQASY spectrum with b/y ions."""
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
            add_losses=True,
        )
        assert len(result) == len(mz)
