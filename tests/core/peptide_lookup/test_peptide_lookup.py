"""Unit + integration tests for core.peptide_lookup (summarize_lookup_results + idXML scan)."""
import pytest
import polars as pl
from psm_utils import PSM, Peptidoform, PSMList

from core.peptide_lookup import (
    load_idxmls,
    lookup_peptides,
    summarize_lookup_results,
)
from tests.conftest import make_peptide_match


def _psm_with_source(run, idxml_path, peptidoform="PEPTIDEK", spectrum_id="index=1"):
    """PSM carrying the idXML provenance the source-aware resolver reads (run + idxml path)."""
    from pathlib import Path
    psm = PSM(peptidoform=Peptidoform(f"{peptidoform}/2"), spectrum_id=spectrum_id, run=run)
    psm.metadata["idxml_resolved_path"] = idxml_path
    psm.metadata["idxml_file"] = Path(idxml_path).name
    return psm


class TestSourceAwareMzmlResolution:
    """lookup_peptides pairs each PSM with the mzML co-located with its source idXML."""

    def test_same_stem_mzml_resolves_per_source_no_shadowing(self):
        # illumina/sample and pac/sample share a stem; each PSM must resolve to the
        # mzML in its own folder rather than one shadowing the other.
        psms = PSMList(psm_list=[
            _psm_with_source("sample", "/illumina/sample.idXML", spectrum_id="index=1"),
            _psm_with_source("sample", "/pac/sample.idXML", spectrum_id="index=2"),
        ])
        grouped = lookup_peptides(
            psms, ["PEPTIDEK"], mzml_paths=["/illumina/sample.mzML", "/pac/sample.mzML"],
        )
        resolved = sorted(
            p.metadata["mzml_resolved_path"] for plist in grouped.values() for p in plist
        )
        assert resolved == ["/illumina/sample.mzML", "/pac/sample.mzML"]

    def test_falls_back_to_global_stem_when_not_colocated(self):
        # idXML and mzML in different folders but the stem is unique → still resolves.
        psms = PSMList(psm_list=[_psm_with_source("sample", "/idx/sample.idXML")])
        grouped = lookup_peptides(psms, ["PEPTIDEK"], mzml_paths=["/mz/sample.mzML"])
        psm = next(iter(grouped.values()))[0]
        assert psm.metadata["mzml_resolved_path"] == "/mz/sample.mzML"

    def test_ambiguous_stem_without_colocation_left_unresolved(self):
        # Stem appears in two folders and none is co-located → refuse to guess.
        psms = PSMList(psm_list=[_psm_with_source("sample", "/idx/sample.idXML")])
        grouped = lookup_peptides(
            psms, ["PEPTIDEK"], mzml_paths=["/illumina/sample.mzML", "/pac/sample.mzML"],
        )
        psm = next(iter(grouped.values()))[0]
        assert "mzml_resolved_path" not in psm.metadata


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

    def test_empty_grouped_psms(self):
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

    def test_counts_distinct_sources_for_same_basename(self):
        # Two PSMs, same peptidoform and same idXML basename but different source
        # folders, must count as 2 sources / 2 distinct files (not collapsed to 1).
        illumina = make_peptide_match(spectrum_id="index=1")
        illumina.metadata["idxml_source"] = "illumina"
        pac = make_peptide_match(spectrum_id="index=2")
        pac.metadata["idxml_source"] = "pac"
        df = summarize_lookup_results({
            "PEPTIDEK": PSMList(psm_list=[illumina, pac]),
        }).collect()
        assert df["sources"][0] == 2
        assert df["idxml_files"][0] == 2


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

    def test_idxml_source_stamped_with_parent_dir(self, lookup_results, idxml_path):
        from pathlib import Path
        expected_source = Path(idxml_path).parent.name
        for psm_list in lookup_results.values():
            for psm in psm_list:
                assert psm.metadata["idxml_source"] == expected_source

    def test_summarize_returns_polars(self, lookup_results):
        summary = summarize_lookup_results(lookup_results)
        assert isinstance(summary, pl.LazyFrame)
        df = summary.collect()
        assert len(df) == 2
        assert "id" in df.columns
