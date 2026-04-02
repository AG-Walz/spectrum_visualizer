"""
Unit and integration tests for src/koina_prediction.py
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.koina_prediction import (
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
        # Simple peptide without mods should come back essentially unchanged
        assert "PEPTIDE" in result

    def test_with_charge_suffix(self):
        result_with = convert_peptide_to_unimod("PEPTIDE/2")
        result_without = convert_peptide_to_unimod("PEPTIDE")
        # Charge suffix should be stripped; both should produce the same UniMod string
        assert result_with == result_without

    def test_oxidation_modification(self):
        result = convert_peptide_to_unimod("PEPTM[Oxidation]DE")
        assert isinstance(result, str)
        # Should contain UNIMOD notation with brackets
        assert "UNIMOD" in result.upper() or "UniMod" in result

    def test_known_psm_asddfdqasy(self):
        result = convert_peptide_to_unimod("ASDDFDQASY")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_known_psm_rtefnlnqy(self):
        result = convert_peptide_to_unimod("RTEFNLNQY")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_empty_string(self):
        # Empty string should either return empty or raise
        try:
            result = convert_peptide_to_unimod("")
            # If it doesn't raise, result should be a string
            assert isinstance(result, str)
        except (ValueError, Exception):
            pass  # Raising is acceptable

    def test_invalid_peptide(self):
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
    def test_koina_available_is_true(self):
        assert KOINA_AVAILABLE is True, "koinapy should be installed in mspecviz env"

    def test_koina_models_not_empty(self):
        assert len(KOINA_MODELS) > 0

    def test_instrument_types_not_empty(self):
        assert len(INSTRUMENT_TYPES) > 0

    def test_default_model_in_models(self):
        assert DEFAULT_MODEL in KOINA_MODELS.values()

    def test_default_instrument_in_types(self):
        assert DEFAULT_INSTRUMENT_TYPE in INSTRUMENT_TYPES

    def test_default_collision_energy_positive(self):
        assert DEFAULT_COLLISION_ENERGY > 0

    def test_model_keys_are_display_names(self):
        # Keys are human-readable display names (contain spaces or capitals)
        for key in KOINA_MODELS:
            assert isinstance(key, str)
            assert len(key) > 0

    def test_model_values_are_ids(self):
        # Values are machine-readable IDs
        for value in KOINA_MODELS.values():
            assert isinstance(value, str)
            assert " " not in value


# ---------------------------------------------------------------------------
# TestPredictKoinaSpectrum (network required)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestPredictKoinaSpectrum:
    def test_predict_simple_peptide(self):
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert result is not None

    def test_predict_second_peptide(self):
        result = predict_koina_spectrum("RTEFNLNQY", charge=2)
        assert result is not None

    def test_result_is_dataframe(self):
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert isinstance(result, pd.DataFrame)

    def test_result_has_required_columns(self):
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert result is not None
        for col in ["mz", "intensity", "annotation"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_result_not_empty(self):
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert result is not None
        assert len(result) > 0

    def test_mz_values_positive(self):
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert result is not None
        assert (result["mz"] > 0).all()

    def test_intensity_values_positive(self):
        # After filtering zero-intensity peaks, all remaining should be > 0
        result = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert result is not None
        assert (result["intensity"] > 0).all()

    def test_invalid_peptide_returns_none(self):
        result = predict_koina_spectrum("ZZZZZZZ", charge=1)
        assert result is None

    def test_model_prosit_2019(self):
        result = predict_koina_spectrum(
            "ASDDFDQASY", charge=2, model_id="Prosit_2019_intensity"
        )
        assert result is not None
        assert len(result) > 0

    def test_model_ms2pip_hcd2021(self):
        result = predict_koina_spectrum(
            "ASDDFDQASY", charge=2, model_id="ms2pip_HCD2021"
        )
        assert result is not None
        assert len(result) > 0

    def test_different_collision_energy(self):
        result = predict_koina_spectrum(
            "RTEFNLNQY", charge=2, collision_energy=25
        )
        assert result is not None
        assert len(result) > 0


# ---------------------------------------------------------------------------
# TestKoinaPipelineIntegration (network required)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestKoinaPipelineIntegration:
    def test_predict_and_annotate(self):
        """Predict a spectrum via Koina, then annotate it with prepare_spectrum_for_plot."""
        from src.peak_annotation import prepare_spectrum_for_plot

        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        annotated = prepare_spectrum_for_plot(
            predicted["mz"].values,
            predicted["intensity"].values,
            "ASDDFDQASY",
            charge=1,
            normalize=True,
        )
        assert isinstance(annotated, pd.DataFrame)
        assert len(annotated) > 0
        assert "annotation" in annotated.columns
        assert "ion_type" in annotated.columns
        # At least some peaks should be annotated (not all unannotated)
        annotated_count = (annotated["is_annotated"]).sum()
        assert annotated_count > 0, "Expected at least some annotated peaks from Koina prediction"

    def test_koina_as_reference_spectrum(self):
        """
        Full flow: predict reference via Koina, annotate both query and reference,
        then compute comparison metrics.
        """
        import numpy as np
        from src.peak_annotation import prepare_spectrum_for_plot
        from src.spectrum_metrics import compute_all_metrics

        # Simulate a query spectrum (synthetic peaks around known ASDDFDQASY fragments)
        query_mz = np.array([175.0, 290.0, 405.0, 520.0, 635.0, 750.0])
        query_intensity = np.array([500.0, 800.0, 1000.0, 600.0, 300.0, 200.0])

        query_df = prepare_spectrum_for_plot(
            query_mz, query_intensity, "ASDDFDQASY", charge=1, normalize=True
        )

        # Predict reference via Koina
        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        ref_df = prepare_spectrum_for_plot(
            predicted["mz"].values,
            predicted["intensity"].values,
            "ASDDFDQASY",
            charge=1,
            normalize=True,
        )

        # Compute metrics between query and reference
        metrics = compute_all_metrics(query_df, ref_df)
        assert isinstance(metrics, dict)
        assert len(metrics) > 0
        # Cosine similarity should be a number between 0 and 1
        if "cosine_similarity" in metrics:
            assert 0.0 <= metrics["cosine_similarity"] <= 1.0
