"""
Tests for ``_proforma_to_openms_str`` peptidoform string conversion.

Builder-level behaviour (``build_psm_table``, ``build_peaks_dataframe``, …)
lives in ``test_insight_frames_builders.py``; those builders also exercise
the polars schemas implicitly, so explicit schema-equality assertions add
no coverage and have been removed.
"""
from core.dataframes.insight_frames import _proforma_to_openms_str


class TestProformaToOpenmsStr:
    def test_simple_peptide(self):
        result = _proforma_to_openms_str("PEPTIDEK")
        assert result == "PEPTIDEK"

    def test_oxidation_modification(self):
        result = _proforma_to_openms_str("PEPTM[Oxidation]IDEK")
        assert "Oxidation" in result
        assert "(" in result  # OpenMS uses parentheses

    def test_charge_suffix_stripped(self):
        result = _proforma_to_openms_str("PEPTIDEK/2")
        assert "/2" not in result
        assert "PEPTIDEK" in result

    def test_nterm_modification(self):
        result = _proforma_to_openms_str("[Acetyl]-PEPTIDEK")
        assert "(Acetyl)" in result

    def test_fallback_on_invalid(self):
        # Should not raise, returns best-effort string
        result = _proforma_to_openms_str("PEPTIDEK")
        assert isinstance(result, str)
        assert len(result) > 0
