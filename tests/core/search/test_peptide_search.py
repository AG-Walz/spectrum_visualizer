"""Tests for ``summarize_search_results`` and the search-pipeline integration."""
import pytest
import polars as pl
from psm_utils import PSMList

from core.search import (
    search_peptides_in_idxmls,
    summarize_search_results,
)
from tests.conftest import make_peptide_match


class TestSummarizeSearchResults:
    """Tests for ``summarize_search_results()`` which returns ``pl.LazyFrame``."""

    def test_returns_lazyframe(self):
        result = summarize_search_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
        })
        assert isinstance(result, pl.LazyFrame)

    def test_expected_columns(self):
        result = summarize_search_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
        })
        df = result.collect()
        expected_cols = {"id", "base_sequence", "modification", "display_label",
                         "psm_count", "idxml_files", "mzml_resolved", "unresolved"}
        assert expected_cols.issubset(set(df.columns))

    def test_psm_count(self):
        result = summarize_search_results({
            "PEPTIDEK": PSMList(psm_list=[
                make_peptide_match(spectrum_id="index=1"),
                make_peptide_match(spectrum_id="index=2"),
            ]),
        })
        df = result.collect()
        assert df["psm_count"][0] == 2

    def test_empty_variant_groups(self):
        result = summarize_search_results({})
        assert isinstance(result, pl.LazyFrame)
        df = result.collect()
        assert len(df) == 0

    def test_multiple_variants(self):
        result = summarize_search_results({
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match()]),
            "AAAVAAAK": PSMList(psm_list=[make_peptide_match(peptidoform="AAAVAAAK")]),
        })
        df = result.collect()
        assert len(df) == 2


# ---------------------------------------------------------------------------
# Integration tests using real test data
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestSearchWithRealData:
    """Tests using real idXML/mzML test data."""

    @pytest.fixture(scope="class")
    def search_results(self, idxml_path, target_peptides):
        """Shared fixture: search once for all tests in this class."""
        return search_peptides_in_idxmls([idxml_path], target_peptides)

    def test_search_finds_expected_peptides(self, search_results):
        assert len(search_results) == 2
        assert "ASDDFDQASY" in search_results
        assert "RTEFNLNQY" in search_results

    def test_each_variant_has_one_psm_with_valid_fields(self, search_results):
        for key, psm_list in search_results.items():
            assert len(psm_list) == 1
            psm = psm_list[0]
            assert psm.peptidoform is not None
            assert psm.spectrum_id is not None
            assert (psm.get_precursor_charge() or 0) > 0

    def test_summarize_returns_polars(self, search_results):
        summary = summarize_search_results(search_results)
        assert isinstance(summary, pl.LazyFrame)
        df = summary.collect()
        assert len(df) == 2
        assert "id" in df.columns
