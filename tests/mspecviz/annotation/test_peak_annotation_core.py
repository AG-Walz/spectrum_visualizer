"""
Unit tests for ``mspecviz.annotation`` — low-level building blocks.

Covers peptidoform parsing, theoretical fragment generation, the peak-to-
fragment matcher itself, and the Koina annotation helpers. End-to-end
pipeline composition (``prepare_spectrum_for_plot``) lives in the
``test_peak_annotation_integration.py`` sibling.
"""
import numpy as np
import pandas as pd
import pytest
import pyopenms as oms

from mspecviz.annotation import (
    annotate_spectrum,
    generate_theoretical_fragments,
    peptidoform_to_aasequence,
    prepare_koina_spectrum_for_plot,
)
from mspecviz.annotation.peak_matching import _parse_koina_charge


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


class TestGenerateFragmentMzDf:
    @pytest.mark.parametrize("ion_type", ["b", "y"])
    def test_ion_type_present(self, ion_type):
        df = generate_theoretical_fragments("PEPTIDEK", charge=2, ion_types=[ion_type])
        matching = df[df["annotation"].str.startswith(ion_type)]
        assert len(matching) > 0

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


class TestParseKoinaCharge:
    @pytest.mark.parametrize("annotation, expected_charge", [
        ("y3+1", 1),
        ("b5+2", 2),
        ("y7+3", 3),
        ("y3++", 2),
        ("b5+++", 3),
        ("y3", 1),
        ("", 0),
        ("y3-H2O+1", 1),  # Koina loss annotation format
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
        """Verify Koina native annotations are preserved with correct ion types and charges."""
        result = prepare_koina_spectrum_for_plot(koina_frame)
        # Required columns present
        for col in ("mz", "intensity", "annotation", "ion_type", "charge", "is_annotated"):
            assert col in result.columns
        # Annotations preserved in Koina format
        assert result["annotation"].iloc[0] == "y3+1"
        assert result["annotation"].iloc[1] == "b5+1"
        # Ion types correctly extracted
        assert result["ion_type"].iloc[0] == "y"
        assert result["ion_type"].iloc[1] == "b"
        # Charges correctly parsed
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
