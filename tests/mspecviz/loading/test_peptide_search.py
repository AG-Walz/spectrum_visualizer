"""Tests for peptide_search helper functions: abbreviate_mod, format_inline_sequence,
and summarize_search_results (returns pl.LazyFrame).
"""
import pytest
import polars as pl

from mspecviz.search import (
    abbreviate_mod,
    format_inline_sequence,
    get_unique_peptides_from_idxmls,
    MOD_ABBREVIATIONS,
    search_peptides_in_idxmls,
    VariantGroup,
    summarize_search_results,
)
from tests.conftest import make_peptide_match


class TestAbbreviateMod:
    """Tests for abbreviate_mod()."""

    def test_known_modifications(self):
        assert abbreviate_mod("Oxidation") == "Ox"
        assert abbreviate_mod("Phospho") == "Ph"
        assert abbreviate_mod("Carbamidomethyl") == "CAM"
        assert abbreviate_mod("Acetyl") == "Ac"
        assert abbreviate_mod("Deamidated") == "Deam"

    def test_unknown_modification_returns_full_name(self):
        assert abbreviate_mod("SomeUnknownMod") == "SomeUnknownMod"

    def test_all_abbreviations_in_map(self):
        for mod_name, abbrev in MOD_ABBREVIATIONS.items():
            assert abbreviate_mod(mod_name) == abbrev


class TestFormatInlineSequence:
    """Tests for format_inline_sequence()."""

    def test_no_modification(self):
        assert format_inline_sequence("PEPTIDEK", "") == "PEPTIDEK"

    def test_single_modification(self):
        result = format_inline_sequence("PEPTMDEK", "M5[Ox]")
        assert result == "PEPTM[Ox]DEK"

    def test_multiple_modifications(self):
        result = format_inline_sequence("PEPTMSEK", "M5[Ox]|S6[Ph]")
        assert result == "PEPTM[Ox]S[Ph]EK"

    def test_nterm_modification(self):
        result = format_inline_sequence("PEPTIDEK", "NTerm[Ac]")
        assert result == "[Ac]PEPTIDEK"

    def test_cterm_modification(self):
        result = format_inline_sequence("PEPTIDEK", "CTerm[Ac]")
        assert result == "PEPTIDEK[Ac]"

    def test_nterm_and_residue_modification(self):
        result = format_inline_sequence("PEPTMDEK", "NTerm[Ac]|M5[Ox]")
        assert result == "[Ac]PEPTM[Ox]DEK"

    def test_modification_at_first_position(self):
        result = format_inline_sequence("MPEPTIDEK", "M1[Ox]")
        assert result == "M[Ox]PEPTIDEK"

    def test_modification_at_last_position(self):
        result = format_inline_sequence("PEPTIDEK", "K8[Ac]")
        assert result == "PEPTIDEK[Ac]"

    def test_unknown_mod_keeps_full_name(self):
        result = format_inline_sequence("PEPTMDEK", "M5[SomeUnknownMod]")
        assert result == "PEPTM[SomeUnknownMod]DEK"


class TestSummarizeSearchResults:
    """Tests for summarize_search_results() which returns pl.LazyFrame."""

    def test_returns_lazyframe(self):
        variant_group = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        result = summarize_search_results({"PEPTIDEK": variant_group})
        assert isinstance(result, pl.LazyFrame)

    def test_expected_columns(self):
        variant_group = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        result = summarize_search_results({"PEPTIDEK": variant_group})
        df = result.collect()
        expected_cols = {"id", "base_sequence", "modification", "display_label",
                         "psm_count", "idxml_files", "mzml_resolved", "unresolved"}
        assert expected_cols.issubset(set(df.columns))

    def test_psm_count(self):
        variant_group = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[
                make_peptide_match(spectrum_id="index=1"),
                make_peptide_match(spectrum_id="index=2"),
            ],
        )
        result = summarize_search_results({"PEPTIDEK": variant_group})
        df = result.collect()
        assert df["psm_count"][0] == 2

    def test_empty_variant_groups(self):
        result = summarize_search_results({})
        assert isinstance(result, pl.LazyFrame)
        df = result.collect()
        assert len(df) == 0

    def test_multiple_variants(self):
        variant_group1 = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        variant_group2 = VariantGroup(
            base_sequence="AAAVAAAK",
            modification_key="",
            display_label="AAAVAAAK",
            matches=[make_peptide_match(peptidoform="AAAVAAAK", base_seq="AAAVAAAK")],
        )
        result = summarize_search_results({"PEPTIDEK": variant_group1, "AAAVAAAK": variant_group2})
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
        for key, variant_group in search_results.items():
            assert len(variant_group.matches) == 1
            match = variant_group.matches[0]
            assert match.peptidoform is not None
            assert match.spectrum_id is not None
            assert match.charge > 0

    def test_get_unique_peptides(self, idxml_path):
        peptides = get_unique_peptides_from_idxmls([idxml_path])
        assert len(peptides) > 0
        assert isinstance(peptides, (list, set))

    def test_summarize_returns_polars(self, search_results):
        summary = summarize_search_results(search_results)
        assert isinstance(summary, pl.LazyFrame)
        df = summary.collect()
        assert len(df) == 2
        assert "id" in df.columns
