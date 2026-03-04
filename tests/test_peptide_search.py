"""Tests for peptide_search helper functions: abbreviate_mod, format_inline_sequence."""
import pytest

from src.peptide_search import abbreviate_mod, format_inline_sequence, MOD_ABBREVIATIONS


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
