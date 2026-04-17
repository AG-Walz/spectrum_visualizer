"""
Unit and integration tests for mspecviz/koina_prediction.py
"""
import numpy as np
import pandas as pd
import pytest

from mspecviz.annotation import (
    compute_all_metrics,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from mspecviz.spectra import (
    KOINA_AVAILABLE,
    KOINA_MODELS,
    INSTRUMENT_TYPES,
    DEFAULT_MODEL,
    DEFAULT_COLLISION_ENERGY,
    DEFAULT_INSTRUMENT_TYPE,
    convert_peptide_to_unimod,
    predict_koina_spectrum,
    get_koina_settings,
)


# ---------------------------------------------------------------------------
# TestConvertPeptideToUnimod
# ---------------------------------------------------------------------------

class TestConvertPeptideToUnimod:
    def test_simple_peptide(self):
        result = convert_peptide_to_unimod("PEPTIDE")
        assert isinstance(result, str)
        assert len(result) > 0
        assert "PEPTIDE" in result

    def test_charge_suffix_stripped(self):
        result_with = convert_peptide_to_unimod("PEPTIDE/2")
        result_without = convert_peptide_to_unimod("PEPTIDE")
        assert result_with == result_without

    def test_oxidation_modification(self):
        result = convert_peptide_to_unimod("PEPTM[Oxidation]DE")
        assert isinstance(result, str)
        assert "UNIMOD" in result.upper() or "UniMod" in result

    @pytest.mark.parametrize("peptide", ["ASDDFDQASY", "RTEFNLNQY"])
    def test_known_psm_peptides(self, peptide):
        result = convert_peptide_to_unimod(peptide)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_empty_string(self):
        try:
            result = convert_peptide_to_unimod("")
            assert isinstance(result, str)
        except (ValueError, Exception):
            pass

    def test_invalid_peptide_raises(self):
        with pytest.raises((ValueError, Exception)):
            convert_peptide_to_unimod("!@#$%^&*")


# ---------------------------------------------------------------------------
# TestGetKoinaSettings
# ---------------------------------------------------------------------------

class TestGetKoinaSettings:
    def test_defaults_from_settings(self, app_settings):
        result = get_koina_settings(app_settings)
        assert result["enabled"] is True
        assert result["model"] == "ms2pip_timsTOF2024"
        assert result["collision_energy"] == 35
        assert result["instrument_type"] == "TIMSTOF"

    def test_empty_settings_returns_defaults(self):
        result = get_koina_settings({})
        assert result["enabled"] is True
        assert result["model"] == DEFAULT_MODEL
        assert result["collision_energy"] == DEFAULT_COLLISION_ENERGY
        assert result["instrument_type"] == DEFAULT_INSTRUMENT_TYPE

    def test_partial_settings(self):
        settings = {"koina": {"model": "Prosit_2019_intensity"}}
        result = get_koina_settings(settings)
        assert result["model"] == "Prosit_2019_intensity"
        # Other keys should fall back to defaults
        assert result["enabled"] is True
        assert result["collision_energy"] == DEFAULT_COLLISION_ENERGY
        assert result["instrument_type"] == DEFAULT_INSTRUMENT_TYPE

    def test_all_settings_override(self):
        settings = {
            "koina": {
                "enabled": False,
                "model": "UniSpec",
                "collision_energy": 25,
                "instrument_type": "QE",
            }
        }
        result = get_koina_settings(settings)
        assert result["enabled"] is False
        assert result["model"] == "UniSpec"
        assert result["collision_energy"] == 25
        assert result["instrument_type"] == "QE"


# ---------------------------------------------------------------------------
# TestKoinaConstants
# ---------------------------------------------------------------------------

class TestKoinaConstants:
    def test_koina_available(self):
        assert KOINA_AVAILABLE is True, "koinapy should be installed in mspecviz env"

    def test_defaults_valid(self):
        assert len(KOINA_MODELS) > 0
        assert len(INSTRUMENT_TYPES) > 0
        assert DEFAULT_MODEL in KOINA_MODELS.values()
        assert DEFAULT_INSTRUMENT_TYPE in INSTRUMENT_TYPES
        assert DEFAULT_COLLISION_ENERGY > 0

    def test_model_keys_are_display_names_and_values_are_ids(self):
        for key, value in KOINA_MODELS.items():
            assert isinstance(key, str) and len(key) > 0
            assert isinstance(value, str) and " " not in value


# ---------------------------------------------------------------------------
# TestPredictKoinaSpectrum (network required)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestPredictKoinaSpectrum:
    @pytest.fixture(scope="class")
    def predicted_asddfdqasy(self):
        """Shared prediction result to avoid redundant network calls."""
        return predict_koina_spectrum("ASDDFDQASY", charge=1)

    def test_prediction_returns_valid_dataframe(self, predicted_asddfdqasy):
        """Verify prediction result is a non-empty DataFrame with required columns and valid values."""
        result = predicted_asddfdqasy
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        for col in ["mz", "intensity", "annotation"]:
            assert col in result.columns, f"Missing column: {col}"
        assert (result["mz"] > 0).all()
        assert (result["intensity"] > 0).all()

    @pytest.mark.parametrize("peptide,charge", [
        ("ASDDFDQASY", 1),
        ("RTEFNLNQY", 2),
    ])
    def test_predict_known_peptides(self, peptide, charge):
        result = predict_koina_spectrum(peptide, charge=charge)
        assert result is not None
        assert len(result) > 0

    def test_invalid_peptide_returns_none(self):
        result = predict_koina_spectrum("ZZZZZZZ", charge=1)
        assert result is None

    @pytest.mark.parametrize("model_id", ["Prosit_2019_intensity", "ms2pip_HCD2021"])
    def test_alternative_models(self, model_id):
        result = predict_koina_spectrum("ASDDFDQASY", charge=2, model_id=model_id)
        assert result is not None
        assert len(result) > 0

    def test_different_collision_energy(self):
        result = predict_koina_spectrum("RTEFNLNQY", charge=2, collision_energy=25)
        assert result is not None
        assert len(result) > 0


# ---------------------------------------------------------------------------
# TestKoinaPipelineIntegration (network required)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestKoinaPipelineIntegration:
    def test_predict_and_annotate_native(self):
        """Predict a spectrum via Koina and annotate using native Koina annotations."""
        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        annotated = prepare_koina_spectrum_for_plot(predicted, normalize=True)
        assert isinstance(annotated, pd.DataFrame)
        assert len(annotated) > 0
        assert "annotation" in annotated.columns
        assert "ion_type" in annotated.columns
        # Koina predictions should all be annotated (they are predicted fragments)
        annotated_count = (annotated["is_annotated"]).sum()
        assert annotated_count > 0, "Expected annotated peaks from Koina prediction"
        # Verify native annotations are preserved (not re-matched to PyOpenMS theory)
        # Koina annotations should be like 'y3+1', 'b5+1', not PyOpenMS format
        non_unannotated = annotated[annotated["annotation"] != "unannotated"]
        if len(non_unannotated) > 0:
            sample_ann = non_unannotated["annotation"].iloc[0]
            assert "+" in sample_ann, f"Expected Koina annotation format (e.g. y3+1), got '{sample_ann}'"

    def test_koina_as_reference_spectrum(self):
        """
        Full flow: predict reference via Koina using native annotations,
        annotate query with PyOpenMS, then compute comparison metrics.
        """
        # Simulate a query spectrum (synthetic peaks around known ASDDFDQASY fragments)
        query_mz = np.array([175.0, 290.0, 405.0, 520.0, 635.0, 750.0])
        query_intensity = np.array([500.0, 800.0, 1000.0, 600.0, 300.0, 200.0])

        query_frame = prepare_spectrum_for_plot(
            query_mz, query_intensity, "ASDDFDQASY", charge=1, normalize=True
        )

        # Predict reference via Koina using native annotations
        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        reference_frame = prepare_koina_spectrum_for_plot(predicted, normalize=True)

        # Compute metrics between query and reference
        metrics = compute_all_metrics(query_frame, reference_frame)
        assert isinstance(metrics, dict)
        assert len(metrics) > 0
        # Cosine similarity should be a number between 0 and 1
        if "cosine_similarity" in metrics:
            assert 0.0 <= metrics["cosine_similarity"] <= 1.0
