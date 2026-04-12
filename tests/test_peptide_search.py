"""Tests for peptide_search helper functions: abbreviate_mod, format_inline_sequence,
and summarize_search_results (returns pl.LazyFrame in MSInsight).
"""
import pytest
import polars as pl

from src.peptide_search import (
    abbreviate_mod,
    format_inline_sequence,
    MOD_ABBREVIATIONS,
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
    """Tests for summarize_search_results() which returns pl.LazyFrame in MSInsight."""

    def test_returns_lazyframe(self):
        vg = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        result = summarize_search_results({"PEPTIDEK": vg})
        assert isinstance(result, pl.LazyFrame)

    def test_expected_columns(self):
        vg = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        result = summarize_search_results({"PEPTIDEK": vg})
        df = result.collect()
        expected_cols = {"id", "base_sequence", "modification", "display_label",
                         "psm_count", "idxml_files", "mzml_resolved", "unresolved"}
        assert expected_cols.issubset(set(df.columns))

    def test_psm_count(self):
        vg = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[
                make_peptide_match(spectrum_id="index=1"),
                make_peptide_match(spectrum_id="index=2"),
            ],
        )
        result = summarize_search_results({"PEPTIDEK": vg})
        df = result.collect()
        assert df["psm_count"][0] == 2

    def test_empty_variant_groups(self):
        result = summarize_search_results({})
        assert isinstance(result, pl.LazyFrame)
        df = result.collect()
        assert len(df) == 0

    def test_multiple_variants(self):
        vg1 = VariantGroup(
            base_sequence="PEPTIDEK",
            modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match()],
        )
        vg2 = VariantGroup(
            base_sequence="AAAVAAAK",
            modification_key="",
            display_label="AAAVAAAK",
            matches=[make_peptide_match(peptidoform="AAAVAAAK", base_seq="AAAVAAAK")],
        )
        result = summarize_search_results({"PEPTIDEK": vg1, "AAAVAAAK": vg2})
        df = result.collect()
        assert len(df) == 2


# ---------------------------------------------------------------------------
# Integration tests using real test data
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestSearchWithRealData:
    """Tests using real idXML/mzML test data."""

    def test_search_finds_two_psms(self, idxml_path, target_peptides):
        """Search known peptides in real idXML file."""
        from src.peptide_search import search_peptides_in_idxmls
        results = search_peptides_in_idxmls([idxml_path], target_peptides)
        assert len(results) == 2
        assert "ASDDFDQASY" in results
        assert "RTEFNLNQY" in results

    def test_each_variant_has_one_psm(self, idxml_path, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        results = search_peptides_in_idxmls([idxml_path], target_peptides)
        for key, vg in results.items():
            assert len(vg.matches) == 1

    def test_psm_has_expected_fields(self, idxml_path, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        results = search_peptides_in_idxmls([idxml_path], target_peptides)
        match = list(results.values())[0].matches[0]
        assert match.peptidoform is not None
        assert match.spectrum_id is not None
        assert match.charge > 0

    def test_get_unique_peptides(self, idxml_path):
        from src.peptide_search import get_unique_peptides_from_idxmls
        peptides = get_unique_peptides_from_idxmls([idxml_path])
        assert len(peptides) > 0
        assert isinstance(peptides, (list, set))

    def test_summarize_returns_polars(self, idxml_path, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        results = search_peptides_in_idxmls([idxml_path], target_peptides)
        summary = summarize_search_results(results)
        assert isinstance(summary, pl.LazyFrame)
        df = summary.collect()
        assert len(df) == 2
        assert "id" in df.columns
