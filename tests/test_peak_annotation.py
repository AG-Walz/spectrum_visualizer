"""
Unit tests for src/peak_annotation.py
"""
import numpy as np
import pandas as pd
import pytest
import pyopenms as oms

from src.peak_annotation import (
    annotate_spectrum,
    generate_theoretical_fragments,
    peptidoform_to_aasequence,
    prepare_spectrum_for_plot,
    prepare_koina_spectrum_for_plot,
    _generate_fragments_impl,
    _parse_koina_charge,
)


# ---------------------------------------------------------------------------
# TestPeptidoformToAASequence
# ---------------------------------------------------------------------------

class TestPeptidoformToAASequence:
    def test_simple_peptide(self):
        seq = peptidoform_to_aasequence("PEPTIDEK")
        assert seq is not None
        assert isinstance(seq, oms.AASequence)

    def test_bracket_notation(self):
        # [Oxidation] should be converted to (Oxidation) for PyOpenMS
        seq = peptidoform_to_aasequence("PEPTM[Oxidation]DEK")
        assert seq is not None
        assert isinstance(seq, oms.AASequence)

    def test_charge_suffix_stripped(self):
        seq_with_charge = peptidoform_to_aasequence("PEPTIDEK/2")
        seq_without = peptidoform_to_aasequence("PEPTIDEK")
        assert seq_with_charge.toString() == seq_without.toString()

    def test_invalid_peptide_raises(self):
        with pytest.raises((ValueError, Exception)):
            peptidoform_to_aasequence("!@#$%^&*")


# ---------------------------------------------------------------------------
# TestGenerateFragmentMzDf
# ---------------------------------------------------------------------------

class TestGenerateFragmentMzDf:
    def test_b_ions_present(self):
        df = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b"])
        b_ions = df[df["annotation"].str.startswith("b")]
        assert len(b_ions) > 0

    def test_y_ions_present(self):
        df = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["y"])
        y_ions = df[df["annotation"].str.startswith("y")]
        assert len(y_ions) > 0

    def test_fragment_count_reasonable(self):
        # PEPTIDEK (8 residues, charge=2) should produce >= 12 fragments
        df = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        assert len(df) >= 12

    def test_cache_returns_same_df(self):
        # lru_cache means the same DataFrame object is returned for identical args
        df1 = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        df2 = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        assert df1 is df2

    def test_different_ion_types_differ(self):
        df_by = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["b", "y"])
        df_aby = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=["a", "b", "y"])
        # Adding "a" ions should give more fragments
        assert len(df_aby) > len(df_by)


# ---------------------------------------------------------------------------
# TestAnnotatePeaks
# ---------------------------------------------------------------------------

class TestAnnotatePeaks:
    @pytest.fixture(autouse=True)
    def setup_frags(self):
        self.theo_frags = generate_theoretical_fragments(
            "PEPTIDEK", charge=2, ion_types=["b", "y"]
        )
        # Get a specific singly-charged b2 ion m/z
        b2_rows = self.theo_frags[self.theo_frags["annotation"] == "b2+"]
        self.b2_mz = float(b2_rows["mz"].values[0])

    def test_exact_match_annotated(self):
        mz = np.array([self.b2_mz])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, self.theo_frags, mz_tolerance=0.02)
        assert result.iloc[0]["annotation"] == "b2+"
        assert bool(result.iloc[0]["is_annotated"]) is True

    def test_unannotated_when_no_match(self):
        # Place a peak far from any fragment (m/z = 999.0, well outside any PEPTIDEK fragment)
        mz = np.array([999.0])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, self.theo_frags, mz_tolerance=0.02)
        assert result.iloc[0]["annotation"] == "unannotated"
        assert bool(result.iloc[0]["is_annotated"]) is False

    @pytest.mark.parametrize("offset, expected_annotated", [
        (0.019, True),    # just within tolerance (tol - 0.001)
        (-0.019, True),   # just within tolerance on the other side
        (0.021, False),   # just outside tolerance (tol + 0.001)
        (-0.021, False),  # just outside tolerance on the other side
        (0.0, True),      # exact match
    ])
    def test_tolerance_boundary(self, offset, expected_annotated):
        mz = np.array([self.b2_mz + offset])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, self.theo_frags, mz_tolerance=0.02)
        assert bool(result.iloc[0]["is_annotated"]) is expected_annotated

    def test_returns_required_columns(self):
        mz = np.array([self.b2_mz])
        intensity = np.array([1000.0])
        result = annotate_spectrum(mz, intensity, self.theo_frags, mz_tolerance=0.02)
        for col in ["mz", "intensity", "annotation", "is_annotated"]:
            assert col in result.columns


# ---------------------------------------------------------------------------
# TestPrepareSpectrumForPlot
# ---------------------------------------------------------------------------

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
        """MSInsight adds an ion_type column to prepare_spectrum_for_plot output."""
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


# ---------------------------------------------------------------------------
# Integration tests with annotation on real-like spectra
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Koina annotation helpers
# ---------------------------------------------------------------------------

class TestParseKoinaCharge:
    def test_plus_number(self):
        assert _parse_koina_charge("y3+1") == 1
        assert _parse_koina_charge("b5+2") == 2
        assert _parse_koina_charge("y7+3") == 3

    def test_double_plus(self):
        assert _parse_koina_charge("y3++") == 2
        assert _parse_koina_charge("b5+++") == 3

    def test_no_charge(self):
        assert _parse_koina_charge("y3") == 1

    def test_empty(self):
        assert _parse_koina_charge("") == 0

    def test_loss_annotation(self):
        # Koina sometimes includes losses like 'y3-H2O+1'
        assert _parse_koina_charge("y3-H2O+1") == 1


class TestPrepareKoinaSpectrumForPlot:
    def _make_koina_df(self):
        return pd.DataFrame({
            "mz": [175.12, 290.14, 405.17, 520.0],
            "intensity": [500.0, 1000.0, 300.0, 800.0],
            "annotation": ["y3+1", "b5+1", "y7+2", "b2+1"],
        })

    def test_preserves_annotations(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df)
        # Annotations should be Koina's native format, not PyOpenMS re-matched
        assert result["annotation"].iloc[0] == "y3+1"
        assert result["annotation"].iloc[1] == "b5+1"

    def test_correct_ion_types(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df)
        assert result["ion_type"].iloc[0] == "y"
        assert result["ion_type"].iloc[1] == "b"

    def test_correct_charges(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df)
        assert result["charge"].iloc[0] == 1
        assert result["charge"].iloc[2] == 2

    def test_normalization(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df, normalize=True)
        assert result["intensity"].max() == pytest.approx(1.0)

    def test_no_normalization(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df, normalize=False)
        assert result["intensity"].max() == 1000.0

    def test_required_columns(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df)
        for col in ("mz", "intensity", "annotation", "ion_type", "charge", "is_annotated"):
            assert col in result.columns

    def test_ion_type_filtering(self):
        koina_df = self._make_koina_df()
        result = prepare_koina_spectrum_for_plot(koina_df, ion_types=["y"])
        # b-ions should be marked unannotated
        b_rows = result[result["annotation"].isin(["b5+1", "b2+1"])]
        assert (b_rows["ion_type"] == "unannotated").all()
        # y-ions should keep their ion_type
        y_rows = result[result["annotation"].isin(["y3+1", "y7+2"])]
        assert (y_rows["ion_type"] == "y").all()
