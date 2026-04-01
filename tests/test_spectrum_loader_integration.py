"""Integration tests for spectrum_loader using real mzML data."""
import pytest
from pathlib import Path


@pytest.mark.integration
class TestLoadRealSpectra:
    def test_streaming_loads_two_spectra(self, idxml_path, mzml_paths, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        from src.spectrum_loader import load_spectra_streaming

        vgs = search_peptides_in_idxmls([idxml_path], target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        loaded = []
        def on_query(result):
            loaded.append(result.query_spectrum)

        load_spectra_streaming(
            variant_groups=vgs, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )
        assert len(loaded) == 2

    def test_spectrum_has_expected_fields(self, idxml_path, mzml_paths, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        from src.spectrum_loader import load_spectra_streaming

        vgs = search_peptides_in_idxmls([idxml_path], target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        spectra = []
        def on_query(result):
            spectra.append(result.query_spectrum)

        load_spectra_streaming(
            variant_groups=vgs, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )

        for spec in spectra:
            assert spec.mz is not None and len(spec.mz) > 0
            assert spec.intensity is not None and len(spec.intensity) > 0
            assert spec.peptidoform is not None
            assert spec.precursor_charge > 0
            assert spec.rt > 0
            assert spec.source_file != ""

    def test_spectrum_lookup_key_unique(self, idxml_path, mzml_paths, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key

        vgs = search_peptides_in_idxmls([idxml_path], target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        keys = []
        def on_query(result):
            keys.append(spectrum_lookup_key(result.query_spectrum))

        load_spectra_streaming(
            variant_groups=vgs, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )
        assert len(keys) == 2
        assert len(set(keys)) == 2  # all unique

    def test_total_peaks_count(self, idxml_path, mzml_paths, target_peptides):
        """Verify total peak count matches expected 373."""
        from src.peptide_search import search_peptides_in_idxmls
        from src.spectrum_loader import load_spectra_streaming

        vgs = search_peptides_in_idxmls([idxml_path], target_peptides)
        mzml_path_list = [Path(p) for p in mzml_paths]

        total_peaks = 0
        def on_query(result):
            nonlocal total_peaks
            total_peaks += len(result.query_spectrum.mz)

        load_spectra_streaming(
            variant_groups=vgs, mzml_paths=mzml_path_list,
            reference_substring="", on_query_ready=on_query,
        )
        assert total_peaks == 373
