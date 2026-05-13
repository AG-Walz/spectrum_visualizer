"""
Schema-level tests for ``core.dataframes``.

Covers schema-equality checks (polars-native via ``frame.collect_schema()``)
and the ``_proforma_to_openms_str`` peptidoform string converter. Builder-level
behaviour (``build_psm_table``, ``build_peaks_dataframe``, ...) lives in
``test_insight_frames_builders.py``.
"""
import polars as pl
import pytest

from core.dataframes import PEAKS_SCHEMA
from core.dataframes.insight_frames import _proforma_to_openms_str


class TestSchemaEquality:
    @pytest.mark.parametrize("columns, expected", [
        # Valid schema
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32)}, True),
        # Missing columns
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32)}, False),
        # Wrong type (Float64 instead of Float32)
        ({"x": pl.Series([1.0], dtype=pl.Float64),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32)}, False),
        # Extra columns
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32),
          "extra_col": [42]}, False),
    ], ids=["valid", "missing_columns", "wrong_type", "extra_columns"])
    def test_schema_validation(self, columns, expected):
        df = pl.DataFrame(columns)
        assert (dict(df.lazy().collect_schema()) == PEAKS_SCHEMA) is expected

    def test_empty_dataframe_valid(self):
        df = pl.DataFrame(schema=PEAKS_SCHEMA)
        assert dict(df.lazy().collect_schema()) == PEAKS_SCHEMA


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
