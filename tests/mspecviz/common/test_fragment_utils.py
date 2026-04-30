"""Unit tests for ``mspecviz.common.fragment_utils.peptidoform_to_aasequence``."""
import pyopenms as oms
import pytest

from mspecviz.common.fragment_utils import peptidoform_to_aasequence


class TestPeptidoformToAASequence:
    def test_simple_peptide(self):
        seq = peptidoform_to_aasequence("PEPTIDEK")
        assert seq is not None
        assert isinstance(seq, oms.AASequence)

    def test_bracket_notation(self):
        # [Oxidation] should be converted to (Oxidation) for PyOpenMS
        seq = peptidoform_to_aasequence("PEPTM[Oxidation]DEK")
        assert seq is not None
        assert isinstance(seq, oms.AASequence)

    def test_charge_suffix_stripped(self):
        seq_with_charge = peptidoform_to_aasequence("PEPTIDEK/2")
        seq_without = peptidoform_to_aasequence("PEPTIDEK")
        assert seq_with_charge.toString() == seq_without.toString()

    def test_invalid_peptide_raises(self):
        with pytest.raises((ValueError, Exception)):
            peptidoform_to_aasequence("!@#$%^&*")
