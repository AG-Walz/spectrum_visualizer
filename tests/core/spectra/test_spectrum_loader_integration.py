"""Integration tests for core.spectra.spectrum_loader against real mzML fixtures."""
import pytest
from pathlib import Path

from core.peptide_lookup import lookup_peptides_in_idxmls
from core.spectra import load_spectra_streaming, psm_lookup_key


@pytest.mark.integration
class TestLoadRealSpectra:
    @pytest.fixture(scope="class")
    def loaded_spectra(self, idxml_path, mzml_paths, target_peptides):
        """Run lookup + streaming once per class — every test below shares the result."""
        variant_groups = lookup_peptides_in_idxmls([idxml_path], target_peptides)
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
            variant_groups=variant_groups, mzml_paths=mzml_path_list,
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
