"""Unit + network-integration tests for ``core.spectra.koina_prediction`` (Koina AI fallback)."""
import numpy as np
import pandas as pd
import pytest

from core.annotation import (
    compute_all_metrics,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.common.defaults import DEFAULTS
from core.spectra import (
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
        assert "PEPTIDE" in result

    def test_charge_suffix_stripped(self):
        assert convert_peptide_to_unimod("PEPTIDE/2") == convert_peptide_to_unimod("PEPTIDE")

    def test_oxidation_modification(self):
        result = convert_peptide_to_unimod("PEPTM[Oxidation]DE")
        assert "UNIMOD" in result.upper() or "UniMod" in result


# ---------------------------------------------------------------------------
# TestGetKoinaSettings
# ---------------------------------------------------------------------------

class TestGetKoinaSettings:
    def test_defaults_from_app_settings(self, app_settings):
        result = get_koina_settings(app_settings)
        assert result["enabled"] is True
        assert result["model"] == "ms2pip_timsTOF2024"
        assert result["collision_energy"] == 35
        assert result["instrument_type"] == "TIMSTOF"

    def test_empty_settings_returns_defaults(self):
        result = get_koina_settings({})
        assert result["model"] == DEFAULTS.koina_model
        assert result["collision_energy"] == DEFAULTS.koina_collision_energy
        assert result["instrument_type"] == DEFAULTS.koina_instrument_type

    def test_full_override(self):
        result = get_koina_settings({
            "koina": {"enabled": False, "model": "UniSpec",
                      "collision_energy": 25, "instrument_type": "QE"},
        })
        assert result == {"enabled": False, "model": "UniSpec",
                          "collision_energy": 25, "instrument_type": "QE"}


# ---------------------------------------------------------------------------
# TestPredictKoinaSpectrum (network required)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestPredictKoinaSpectrum:
    @pytest.fixture(scope="class")
    def predicted_asddfdqasy(self):
        """Run one Koina prediction per class to keep network calls down."""
        return predict_koina_spectrum("ASDDFDQASY", charge=1)

    def test_prediction_returns_valid_dataframe(self, predicted_asddfdqasy):
        """Prediction is a non-empty DataFrame with mz/intensity/annotation columns and positive values."""
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
        """Predict via Koina then preserve native ``y3+1``-style annotations through ``prepare_koina_spectrum_for_plot``."""
        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        annotated = prepare_koina_spectrum_for_plot(predicted, normalize=True)
        assert isinstance(annotated, pd.DataFrame)
        assert len(annotated) > 0
        assert "annotation" in annotated.columns
        assert "ion_type" in annotated.columns
        # Every Koina-predicted peak is a theoretical fragment ⇒ all annotated
        annotated_count = (annotated["is_annotated"]).sum()
        assert annotated_count > 0, "Expected annotated peaks from Koina prediction"
        # Native annotations (y3+1, b5+1) must survive — we don't re-match to PyOpenMS
        non_unannotated = annotated[annotated["annotation"] != "unannotated"]
        if len(non_unannotated) > 0:
            sample_ann = non_unannotated["annotation"].iloc[0]
            assert "+" in sample_ann, f"Expected Koina annotation format (e.g. y3+1), got '{sample_ann}'"

    def test_koina_as_reference_spectrum(self):
        """End-to-end: predict reference via Koina, annotate query via PyOpenMS, compute metrics."""
        # Arrange: synthetic query peaks near the known ASDDFDQASY fragments
        query_mz = np.array([175.0, 290.0, 405.0, 520.0, 635.0, 750.0])
        query_intensity = np.array([500.0, 800.0, 1000.0, 600.0, 300.0, 200.0])

        query_frame = prepare_spectrum_for_plot(
            query_mz, query_intensity, "ASDDFDQASY", charge=1, normalize=True
        )

        # Act: predict reference via Koina (keeps native annotations)
        predicted = predict_koina_spectrum("ASDDFDQASY", charge=1)
        assert predicted is not None

        reference_frame = prepare_koina_spectrum_for_plot(predicted, normalize=True)

        # Assert: comparison metrics dict is populated and cosine is bounded
        metrics = compute_all_metrics(query_frame, reference_frame)
        assert isinstance(metrics, dict)
        assert len(metrics) > 0
        if "cosine_similarity" in metrics:
            assert 0.0 <= metrics["cosine_similarity"] <= 1.0
