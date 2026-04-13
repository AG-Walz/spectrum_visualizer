"""Integration tests for spectrum_loader using real mzML data."""
import pytest
from pathlib import Path


@pytest.mark.integration
class TestLoadRealSpectra:
    @pytest.fixture(scope="class")
    def loaded_spectra(self, idxml_path, mzml_paths, target_peptides):
        """Shared fixture: search + stream spectra once for all tests in this class."""
        from src.peptide_search import search_peptides_in_idxmls
        from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key

        vgs = search_peptides_in_idxmls([idxml_path], target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        spectra = []
        keys = []
        total_peaks = 0

        def on_query(result):
            nonlocal total_peaks
            spectra.append(result.query_spectrum)
            keys.append(spectrum_lookup_key(result.query_spectrum))
            total_peaks += len(result.query_spectrum.mz)

        load_spectra_streaming(
            variant_groups=vgs, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )
        return {"spectra": spectra, "keys": keys, "total_peaks": total_peaks}

    def test_streaming_loads_two_spectra(self, loaded_spectra):
        assert len(loaded_spectra["spectra"]) == 2

    def test_spectrum_has_expected_fields(self, loaded_spectra):
        for spec in loaded_spectra["spectra"]:
            assert spec.mz is not None and len(spec.mz) > 0
            assert spec.intensity is not None and len(spec.intensity) > 0
            assert spec.peptidoform is not None
            assert spec.precursor_charge > 0
            assert spec.rt > 0
            assert spec.source_file != ""

    def test_spectrum_lookup_keys_unique(self, loaded_spectra):
        keys = loaded_spectra["keys"]
        assert len(keys) == 2
        assert len(set(keys)) == 2

    def test_total_peaks_count(self, loaded_spectra):
        assert loaded_spectra["total_peaks"] == 373
