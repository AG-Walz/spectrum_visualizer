"""Unit tests for core.annotation.peak_matching (PyOpenMS + Koina annotators)."""
import numpy as np
import pandas as pd
import pytest

from core.annotation import (
    annotate_spectrum,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.annotation.peak_matching import _parse_koina_charge

# Singly-charged b2 / y2 m/z for ASDDFDQASY — pinned (computed via
# spectrum_utils during the parity spike) so the test isn't coupled to whichever
# theoretical-spectrum library we use today.
ASDDFDQASY_B2 = 159.0764
ASDDFDQASY_Y2 = 269.1132


class TestAnnotatePeaks:
    """annotate_spectrum — match peaks against theoretical fragments at a given tolerance."""

    def test_exact_match_annotated(self):
        mz = np.array([ASDDFDQASY_B2])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, "ASDDFDQASY", charge=1, mz_tolerance=0.02)
        assert bool(result.iloc[0]["is_annotated"]) is True
        assert result.iloc[0]["ion_type"] == "b"

    def test_unannotated_when_no_match(self):
        mz = np.array([999.0])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, "ASDDFDQASY", charge=1, mz_tolerance=0.02)
        assert result.iloc[0]["annotation"] == "unannotated"
        assert bool(result.iloc[0]["is_annotated"]) is False

    @pytest.mark.parametrize("offset, expected_annotated", [
        (0.019, True),
        (-0.019, True),
        (0.021, False),
        (-0.021, False),
        (0.0, True),
    ])
    def test_tolerance_boundary(self, offset, expected_annotated):
        mz = np.array([ASDDFDQASY_B2 + offset])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, "ASDDFDQASY", charge=1, mz_tolerance=0.02)
        assert bool(result.iloc[0]["is_annotated"]) is expected_annotated

    def test_empty_input_returns_empty_df(self):
        result = annotate_spectrum(
            np.array([]), np.array([]), "ASDDFDQASY", charge=1
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


class TestPrepareSpectrumForPlot:
    """prepare_spectrum_for_plot — annotate + filter + normalise composition."""

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


class TestAnnotateRealisticSpectrum:
    """End-to-end checks on realistic ASDDFDQASY / RTEFNLNQY peak lists."""

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


class TestParseKoinaCharge:
    @pytest.mark.parametrize("annotation, expected_charge", [
        ("y3+1", 1),
        ("b5+2", 2),
        ("y7+3", 3),
        ("y3++", 2),
        ("b5+++", 3),
        ("y3", 1),
        ("", 0),
        ("y3-H2O+1", 1),  # Koina neutral-loss annotation format
    ])
    def test_charge_parsing(self, annotation, expected_charge):
        assert _parse_koina_charge(annotation) == expected_charge


class TestPrepareKoinaSpectrumForPlot:
    @pytest.fixture
    def koina_frame(self):
        return pd.DataFrame({
            "mz": [175.12, 290.14, 405.17, 520.0],
            "intensity": [500.0, 1000.0, 300.0, 800.0],
            "annotation": ["y3+1", "b5+1", "y7+2", "b2+1"],
        })

    def test_annotations_ion_types_and_charges(self, koina_frame):
        """Koina native annotations are preserved with correct ion types and charges."""
        result = prepare_koina_spectrum_for_plot(koina_frame)
        for col in ("mz", "intensity", "annotation", "ion_type", "charge", "is_annotated"):
            assert col in result.columns
        assert result["annotation"].iloc[0] == "y3+1"
        assert result["annotation"].iloc[1] == "b5+1"
        assert result["ion_type"].iloc[0] == "y"
        assert result["ion_type"].iloc[1] == "b"
        assert result["charge"].iloc[0] == 1
        assert result["charge"].iloc[2] == 2

    @pytest.mark.parametrize("normalize, expected_max", [
        (True, 1.0),
        (False, 1000.0),
    ])
    def test_normalization(self, koina_frame, normalize, expected_max):
        result = prepare_koina_spectrum_for_plot(koina_frame, normalize=normalize)
        assert result["intensity"].max() == pytest.approx(expected_max)

    def test_ion_type_filtering(self, koina_frame):
        result = prepare_koina_spectrum_for_plot(koina_frame, ion_types=["y"])
        b_rows = result[result["annotation"].isin(["b5+1", "b2+1"])]
        assert (b_rows["ion_type"] == "unannotated").all()
        y_rows = result[result["annotation"].isin(["y3+1", "y7+2"])]
        assert (y_rows["ion_type"] == "y").all()
