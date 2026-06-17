"""Unit tests for native-ID -> spectrum resolution in core.spectra.spectrum_streaming.

Builds synthetic mzML files in tmp_path (no external fixtures), so these run in any
environment. They lock in the two-tier resolver inside _load_spectra_from_mzml:
the pyopenms getSpectrumByNativeId exact path and the SpectrumLookup.findByScanNumber
fallback for references that carry only the scan number.
"""
import pyopenms as oms
import pytest
from psm_utils import PSM

from core.spectra.spectrum_streaming import _load_spectra_from_mzml


def _write_mzml(path, native_ids):
    """Write an mzML whose spectra carry the given native IDs; peak mz0 encodes the index (100+i)."""
    experiment = oms.MSExperiment()
    for i, native_id in enumerate(native_ids):
        spectrum = oms.MSSpectrum()
        spectrum.setNativeID(native_id)
        spectrum.setMSLevel(2)
        spectrum.setRT(float(i))
        spectrum.set_peaks(([100.0 + i, 200.0 + i, 300.0 + i], [1.0, 2.0, 3.0]))
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)
    return str(path)


def _matches(*spectrum_ids):
    return [(PSM(peptidoform="PEPTIDEK/2", spectrum_id=sid, run="r"), "PEPTIDEK/2") for sid in spectrum_ids]


class TestExactNativeIdResolution:
    """getSpectrumByNativeId resolves the full native-ID string directly, format-agnostic."""

    def test_bruker_one_based_index(self, tmp_path):
        # index=1 sits at array position 0 — the timsTOF 1-based gotcha; exact string match
        # resolves it with no offset arithmetic.
        path = _write_mzml(tmp_path / "bruker.mzML", [f"index={i}" for i in range(1, 6)])
        result = _load_spectra_from_mzml(path, "bruker.mzML", _matches("index=3"), is_reference=False)
        assert result.loaded_count == 1
        assert result.not_found_count == 0
        # index=3 -> array position 2 -> mz0 == 100 + 2
        assert result.peptidoform_psm_pairs[0][1].spectrum["mz"][0] == pytest.approx(102.0)

    def test_thermo_controller_scan(self, tmp_path):
        ids = [f"controllerType=0 controllerNumber=1 scan={i}" for i in range(1, 6)]
        path = _write_mzml(tmp_path / "thermo.mzML", ids)
        result = _load_spectra_from_mzml(path, "thermo.mzML", _matches(ids[3]), is_reference=False)
        assert result.loaded_count == 1
        assert result.peptidoform_psm_pairs[0][1].spectrum["mz"][0] == pytest.approx(103.0)


class TestScanNumberFallback:
    """A reference carrying only the scan number falls back to SpectrumLookup.findByReference."""

    def test_bare_scan_reference_against_full_native_ids(self, tmp_path):
        # mzML stores full Thermo IDs; the PSM reference is the bare 'scan=3' — exact match misses,
        # the fallback resolves it by scan number.
        ids = [f"controllerType=0 controllerNumber=1 scan={i}" for i in range(1, 6)]
        path = _write_mzml(tmp_path / "partial.mzML", ids)
        result = _load_spectra_from_mzml(path, "partial.mzML", _matches("scan=3"), is_reference=True)
        assert result.loaded_count == 1
        assert result.not_found_count == 0
        assert result.peptidoform_psm_pairs[0][1].spectrum["mz"][0] == pytest.approx(102.0)

    def test_fallback_resolves_same_spectrum_as_exact_index(self, tmp_path):
        # Off-by-one guard for the timsTOF 1-based index= gotcha: the findByReference fallback must
        # land on the SAME array position the exact native-ID path would. The file stores 1-based
        # index= IDs; a bare '3' reference misses the exact path and resolves by scan number to index=3.
        path = _write_mzml(tmp_path / "idx.mzML", [f"index={i}" for i in range(1, 6)])
        exact = _load_spectra_from_mzml(path, "idx.mzML", _matches("index=3"), is_reference=False)
        fallback = _load_spectra_from_mzml(path, "idx.mzML", _matches("3"), is_reference=True)
        assert exact.loaded_count == 1 and fallback.loaded_count == 1
        # index=3 -> array position 2 -> mz0 == 102.0, identical peaks via both resolution paths.
        assert fallback.peptidoform_psm_pairs[0][1].spectrum["mz"][0] == pytest.approx(102.0)
        assert list(fallback.peptidoform_psm_pairs[0][1].spectrum["mz"]) == \
            list(exact.peptidoform_psm_pairs[0][1].spectrum["mz"])


class TestMissingSpectrum:
    """An absent native ID is recorded, not raised — one miss never aborts the batch."""

    def test_missing_id_recorded_as_not_found(self, tmp_path):
        path = _write_mzml(tmp_path / "present.mzML", [f"index={i}" for i in range(1, 4)])
        result = _load_spectra_from_mzml(
            path, "present.mzML", _matches("index=2", "index=999"), is_reference=False,
        )
        assert result.loaded_count == 1
        assert result.not_found_count == 1
        assert any("index=999" in error for error in result.errors)
