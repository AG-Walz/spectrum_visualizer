"""Integration tests for core.spectra.spectrum_loader against real mzML fixtures."""
import pytest
from pathlib import Path

from psm_utils import PSM, Peptidoform, PSMList

from core.peptide_lookup import load_idxmls, lookup_peptides
from core.spectra import load_spectra_streaming, psm_lookup_key
from core.spectra.spectrum_streaming import build_mzml_to_matches_map


def _resolved_key(path):
    return str(Path(path).resolve()).lower()


def _psm(run, mzml_resolved="", idxml_path="", spectrum_id="index=1"):
    psm = PSM(peptidoform=Peptidoform("PEPTIDEK/2"), spectrum_id=spectrum_id, run=run)
    if mzml_resolved:
        psm.metadata["mzml_resolved_path"] = mzml_resolved
    if idxml_path:
        psm.metadata["idxml_resolved_path"] = idxml_path
        psm.metadata["idxml_file"] = Path(idxml_path).name
    return psm


class TestBuildMzmlToMatchesMap:
    """Run-to-mzML bucketing must not cross-contaminate same-stem mzMLs from different folders."""

    def test_pre_resolved_paths_bucket_to_distinct_files(self):
        # Both runs are "sample" but each PSM carries its own resolved path.
        available = {
            _resolved_key("/illumina/sample.mzML"): Path("/illumina/sample.mzML"),
            _resolved_key("/pac/sample.mzML"): Path("/pac/sample.mzML"),
        }
        grouped = {"PEPTIDEK": PSMList(psm_list=[
            _psm("sample", mzml_resolved="/illumina/sample.mzML", spectrum_id="index=1"),
            _psm("sample", mzml_resolved="/pac/sample.mzML", spectrum_id="index=2"),
        ])}
        result = build_mzml_to_matches_map(grouped, available)
        assert set(result) == set(available)
        assert all(len(matches) == 1 for matches in result.values())

    def test_ambiguous_stem_fallback_dropped_without_colocation(self):
        # No pre-resolved path, stem "sample" maps to two files, idXML co-located with
        # neither → the run is dropped rather than silently mis-resolved.
        available = {
            _resolved_key("/illumina/sample.mzML"): Path("/illumina/sample.mzML"),
            _resolved_key("/pac/sample.mzML"): Path("/pac/sample.mzML"),
        }
        grouped = {"PEPTIDEK": PSMList(psm_list=[
            _psm("sample", idxml_path="/idx/sample.idXML"),
        ])}
        result = build_mzml_to_matches_map(grouped, available)
        assert result == {}

    def test_pre_resolved_path_mismatch_falls_back_to_unique_stem(self):
        # mzml_resolved_path is a relative/stale path that does NOT canonicalize into
        # `available`, so the pre-resolved branch misses; the unique stem match still
        # resolves the PSM. This is the load-bearing fallback that protects against
        # path-canonicalization drift between the lookup and streaming stages — the
        # invariant any future shared-resolver refactor must preserve.
        available = {_resolved_key("/data/sample.mzML"): Path("/data/sample.mzML")}
        grouped = {"PEPTIDEK": PSMList(psm_list=[
            _psm("sample", mzml_resolved="sample.mzML"),  # relative -> different key
        ])}
        result = build_mzml_to_matches_map(grouped, available)
        assert set(result) == set(available)
        assert sum(len(matches) for matches in result.values()) == 1

    def test_ambiguous_stem_resolves_to_colocated_mzml(self):
        # Stem "sample" maps to two files; the idXML is co-located with the pac/ one,
        # so the PSM resolves there rather than being dropped (the success path that
        # complements test_ambiguous_stem_fallback_dropped_without_colocation).
        available = {
            _resolved_key("/illumina/sample.mzML"): Path("/illumina/sample.mzML"),
            _resolved_key("/pac/sample.mzML"): Path("/pac/sample.mzML"),
        }
        grouped = {"PEPTIDEK": PSMList(psm_list=[
            _psm("sample", idxml_path="/pac/sample.idXML"),
        ])}
        result = build_mzml_to_matches_map(grouped, available)
        assert set(result) == {_resolved_key("/pac/sample.mzML")}


@pytest.mark.integration
class TestLoadRealSpectra:
    @pytest.fixture(scope="class")
    def loaded_spectra(self, idxml_path, mzml_paths, target_peptides):
        """Run lookup + streaming once per class — every test below shares the result."""
        peptidoform_psms = lookup_peptides(load_idxmls([idxml_path]), target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        spectra = []
        keys = []
        total_peaks = 0

        def on_query(result):
            nonlocal total_peaks
            spectra.append(result.query_psm)
            keys.append(psm_lookup_key(result.query_psm))
            total_peaks += len(result.query_psm.spectrum["mz"])

        load_spectra_streaming(
            peptidoform_psms=peptidoform_psms, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )
        return {"spectra": spectra, "keys": keys, "total_peaks": total_peaks}

    def test_streaming_loads_two_spectra(self, loaded_spectra):
        assert len(loaded_spectra["spectra"]) == 2

    def test_spectrum_has_expected_fields(self, loaded_spectra):
        for psm in loaded_spectra["spectra"]:
            assert len(psm.spectrum["mz"]) > 0
            assert len(psm.spectrum["intensity"]) > 0
            assert psm.peptidoform.modified_sequence
            assert psm.get_precursor_charge() > 0
            assert psm.retention_time > 0
            assert psm.run != ""

    def test_psm_lookup_keys_unique(self, loaded_spectra):
        keys = loaded_spectra["keys"]
        assert len(keys) == 2
        assert len(set(keys)) == 2

    def test_total_peaks_count(self, loaded_spectra):
        # Range check rather than exact total — spectrum-loader behaviour
        # shifts by a handful of peaks with library upgrades without
        # signalling a real regression.
        assert 350 < loaded_spectra["total_peaks"] < 400
