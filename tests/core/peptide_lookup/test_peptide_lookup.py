"""Unit + integration tests for core.peptide_lookup (summarize_lookup_results + idXML scan)."""
import pytest
import polars as pl
from psm_utils import PSMList

from core.peptide_lookup import (
    load_idxmls,
    lookup_peptides,
    summarize_lookup_results,
)
from tests.conftest import make_peptide_match


class TestSummarizeLookupResults:
    """Pure-data tests for summarize_lookup_results (LazyFrame columns, counts, empties)."""

    def test_returns_lazyframe(self):
        result = summarize_lookup_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
        })
        assert isinstance(result, pl.LazyFrame)

    def test_expected_columns(self):
        result = summarize_lookup_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
        })
        df = result.collect()
        expected_cols = {"id", "base_sequence", "modification", "display_label",
                         "psm_count", "idxml_files", "mzml_resolved", "unresolved"}
        assert expected_cols.issubset(set(df.columns))

    def test_psm_count(self):
        result = summarize_lookup_results({
            "PEPTIDEK": PSMList(psm_list=[
                make_peptide_match(spectrum_id="index=1"),
                make_peptide_match(spectrum_id="index=2"),
            ]),
        })
        df = result.collect()
        assert df["psm_count"][0] == 2

    def test_empty_peptidoform_psms(self):
        result = summarize_lookup_results({})
        assert isinstance(result, pl.LazyFrame)
        df = result.collect()
        assert len(df) == 0

    def test_multiple_peptidoforms(self):
        result = summarize_lookup_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
            "AAAVAAAK": PSMList(psm_list=[make_peptide_match(peptidoform="AAAVAAAK")]),
        })
        df = result.collect()
        assert len(df) == 2


# ---------------------------------------------------------------------------
# Integration tests using real test data
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestLookupWithRealData:
    """Look up the two expected peptides in the small_dataset idXML."""

    @pytest.fixture(scope="class")
    def lookup_results(self, idxml_path, target_peptides):
        """Load the idXML once, then scan the PSMList — every test below shares the result."""
        return lookup_peptides(load_idxmls([idxml_path]), target_peptides)

    def test_lookup_finds_expected_peptides(self, lookup_results):
        assert len(lookup_results) == 2
        assert "ASDDFDQASY" in lookup_results
        assert "RTEFNLNQY" in lookup_results

    def test_each_peptidoform_has_one_psm_with_valid_fields(self, lookup_results):
        for key, psm_list in lookup_results.items():
            assert len(psm_list) == 1
            psm = psm_list[0]
            assert psm.peptidoform is not None
            assert psm.spectrum_id is not None
            assert (psm.get_precursor_charge() or 0) > 0

    def test_summarize_returns_polars(self, lookup_results):
        summary = summarize_lookup_results(lookup_results)
        assert isinstance(summary, pl.LazyFrame)
        df = summary.collect()
        assert len(df) == 2
        assert "id" in df.columns
