"""Unit tests for core.dataframes._proforma_to_openms_str (peptidoform string conversion)."""
from core.dataframes import _proforma_to_openms_str


class TestProformaToOpenmsStr:
    def test_simple_peptide(self):
        result = _proforma_to_openms_str("PEPTIDEK")
        assert result == "PEPTIDEK"

    def test_oxidation_modification(self):
        result = _proforma_to_openms_str("PEPTM[Oxidation]IDEK")
        assert "Oxidation" in result
        # OpenMS notation uses parentheses (e.g. M(Oxidation)), not brackets
        assert "(" in result

    def test_charge_suffix_stripped(self):
        result = _proforma_to_openms_str("PEPTIDEK/2")
        assert "/2" not in result
        assert "PEPTIDEK" in result

    def test_nterm_modification(self):
        result = _proforma_to_openms_str("[Acetyl]-PEPTIDEK")
        assert "(Acetyl)" in result

    def test_fallback_on_invalid(self):
        # Unparseable ProForma must not raise — falls back to the bare sequence
        # (charge suffix dropped) instead of the parsed per-residue path.
        result = _proforma_to_openms_str("X1Z9/3")
        assert result == "X1Z9"
