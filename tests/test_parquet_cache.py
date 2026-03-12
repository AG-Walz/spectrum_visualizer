"""
Unit tests for src/parquet_cache.py — automatic Parquet caching for mzML spectra.
"""
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.parquet_cache import (
    compute_cache_fingerprint,
    find_existing_cache,
    get_parquet_row_counts,
    save_spectra_to_parquet,
)
from src.spectrum_loader import LoadedSpectrum


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def tmp_dir(tmp_path):
    """Provide a temporary directory."""
    return tmp_path


@pytest.fixture
def sample_files(tmp_path):
    """Create temporary files that can be used as idXML/mzML stand-ins."""
    idxml = tmp_path / "test.idXML"
    idxml.write_text("<idXML/>")

    mzml = tmp_path / "test.mzML"
    mzml.write_text("<mzML/>")

    return str(idxml), str(mzml)


@pytest.fixture
def sample_spectra():
    """Create sample LoadedSpectrum objects for testing."""
    ref = LoadedSpectrum(
        spectrum_id="scan=100",
        mz=np.array([100.0, 200.0, 300.0]),
        intensity=np.array([1000.0, 2000.0, 500.0]),
        precursor_mz=450.0,
        precursor_charge=2,
        peptidoform="PEPTIDE",
        base_sequence="PEPTIDE",
        source_file="ref.mzML",
        is_reference=True,
        rt=10.5,
    )

    query = LoadedSpectrum(
        spectrum_id="scan=200",
        mz=np.array([101.0, 201.0, 301.0, 401.0]),
        intensity=np.array([800.0, 1500.0, 600.0, 200.0]),
        precursor_mz=450.0,
        precursor_charge=2,
        peptidoform="PEPTIDE",
        base_sequence="PEPTIDE",
        source_file="query.mzML",
        is_reference=False,
        rt=12.3,
    )

    return ref, query


# ============================================================================
# Fingerprint tests
# ============================================================================

class TestFingerprint:
    def test_deterministic(self, sample_files):
        """Same inputs produce same fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "Synpep")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "Synpep")
        assert fp1 == fp2
        assert len(fp1) == 12

    def test_changes_on_different_peptides(self, sample_files):
        """Different peptides produce different fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "Synpep")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["OTHERPEP"], "Synpep")
        assert fp1 != fp2

    def test_changes_on_different_ref_substring(self, sample_files):
        """Different reference substring produces different fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "Synpep")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "Standard")
        assert fp1 != fp2

    def test_order_independent_peptides(self, sample_files):
        """Shuffled peptide order produces same fingerprint (sorted internally)."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["AAA", "BBB", "CCC"], "")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["CCC", "AAA", "BBB"], "")
        assert fp1 == fp2

    def test_case_insensitive_peptides(self, sample_files):
        """Peptide case is normalized (uppercased)."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["peptide"], "")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "")
        assert fp1 == fp2

    def test_changes_on_idxml_filter(self, sample_files):
        """Different idxml_filter produces different fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", idxml_filter="")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", idxml_filter="filtered")
        assert fp1 != fp2

    def test_changes_on_mzml_filter(self, sample_files):
        """Different mzml_filter produces different fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", mzml_filter="")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", mzml_filter="standard")
        assert fp1 != fp2

    def test_filter_case_insensitive(self, sample_files):
        """Filter string is case-normalised in the fingerprint."""
        idxml, mzml = sample_files
        fp1 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", idxml_filter="Filtered")
        fp2 = compute_cache_fingerprint([idxml], [mzml], ["PEPTIDE"], "", idxml_filter="filtered")
        assert fp1 == fp2

    def test_order_independent_files(self, tmp_path):
        """Shuffled file order produces same fingerprint (sorted internally)."""
        f1 = tmp_path / "a.mzML"
        f2 = tmp_path / "b.mzML"
        f1.write_text("data1")
        f2.write_text("data2")

        fp1 = compute_cache_fingerprint([], [str(f1), str(f2)], ["PEP"], "")
        fp2 = compute_cache_fingerprint([], [str(f2), str(f1)], ["PEP"], "")
        assert fp1 == fp2


# ============================================================================
# Cache lookup tests
# ============================================================================

class TestFindCache:
    def test_found(self, tmp_dir):
        """Returns path when cache file exists."""
        cache_file = tmp_dir / "abc123def456.parquet"
        cache_file.write_bytes(b"parquet-data")

        result = find_existing_cache(str(tmp_dir), "abc123def456")
        assert result == str(cache_file)

    def test_not_found(self, tmp_dir):
        """Returns None when cache file doesn't exist."""
        result = find_existing_cache(str(tmp_dir), "nonexistent123")
        assert result is None

    def test_empty_file_not_found(self, tmp_dir):
        """Returns None when cache file exists but is empty."""
        cache_file = tmp_dir / "abc123def456.parquet"
        cache_file.write_bytes(b"")

        result = find_existing_cache(str(tmp_dir), "abc123def456")
        assert result is None


# ============================================================================
# Save and load round-trip tests
# ============================================================================

class TestSaveAndLoad:
    def test_round_trip(self, tmp_dir, sample_spectra):
        """Save spectra then load and verify data matches."""
        ref, query = sample_spectra
        output_path = str(tmp_dir / "cache.parquet")

        reference_spectra = {"PEP|mod": [ref]}
        query_spectra = [("PEP|mod", query)]

        result = save_spectra_to_parquet(
            reference_spectra=reference_spectra,
            query_spectra=query_spectra,
            output_path=output_path,
            cache_metadata={"test": True},
        )

        assert result is not None
        assert Path(result).exists()

        # Reload and verify
        import pyarrow.parquet as pq
        table = pq.read_table(result)
        df = table.to_pandas()

        assert len(df) == 2  # 1 ref + 1 query
        assert df["is_reference"].sum() == 1

        # Verify reference row
        ref_row = df[df["is_reference"]].iloc[0]
        assert ref_row["spectrum_id"] == "scan=100"
        assert ref_row["peptidoform"] == "PEPTIDE"
        assert ref_row["precursor_mz"] == 450.0
        assert ref_row["charge"] == 2
        np.testing.assert_array_almost_equal(ref_row["mz"], [100.0, 200.0, 300.0])
        np.testing.assert_array_almost_equal(ref_row["intensity"], [1000.0, 2000.0, 500.0])

        # Verify query row
        query_row = df[~df["is_reference"]].iloc[0]
        assert query_row["spectrum_id"] == "scan=200"
        assert len(query_row["mz"]) == 4

    def test_empty_spectra_returns_none(self, tmp_dir):
        """Empty input returns None."""
        output_path = str(tmp_dir / "empty.parquet")
        result = save_spectra_to_parquet(
            reference_spectra={},
            query_spectra=[],
            output_path=output_path,
        )
        assert result is None

    def test_sidecar_json_written(self, tmp_dir, sample_spectra):
        """Metadata JSON sidecar is created alongside parquet file."""
        ref, query = sample_spectra
        output_path = str(tmp_dir / "cache.parquet")

        save_spectra_to_parquet(
            reference_spectra={"PEP|mod": [ref]},
            query_spectra=[("PEP|mod", query)],
            output_path=output_path,
            cache_metadata={"fingerprint": "test123", "custom_key": "value"},
        )

        sidecar_path = Path(output_path).with_suffix(".json")
        assert sidecar_path.exists()

        with open(sidecar_path) as f:
            sidecar = json.load(f)

        assert sidecar["total_spectra"] == 2
        assert sidecar["reference_spectra"] == 1
        assert sidecar["query_spectra"] == 1
        assert sidecar["fingerprint"] == "test123"
        assert sidecar["custom_key"] == "value"
        assert "created_at" in sidecar

    def test_parquet_schema_compatibility(self, tmp_dir, sample_spectra):
        """Verify saved Parquet uses the same schema as convert_mzml_to_parquet.py."""
        from src.spectrum_loader import PARQUET_COLUMNS

        ref, query = sample_spectra
        output_path = str(tmp_dir / "cache.parquet")

        save_spectra_to_parquet(
            reference_spectra={"PEP|mod": [ref]},
            query_spectra=[("PEP|mod", query)],
            output_path=output_path,
        )

        import pyarrow.parquet as pq
        table = pq.read_table(output_path)
        parquet_cols = table.column_names

        for col in PARQUET_COLUMNS:
            assert col in parquet_cols, f"Missing column: {col}"


# ============================================================================
# Row count tests
# ============================================================================

class TestGetRowCounts:
    def test_correct_split(self, tmp_dir, sample_spectra):
        """get_parquet_row_counts returns correct ref/query counts."""
        ref, query = sample_spectra
        output_path = str(tmp_dir / "counts.parquet")

        # Save 2 refs and 3 queries
        ref2 = LoadedSpectrum(
            spectrum_id="scan=101", mz=np.array([150.0]), intensity=np.array([500.0]),
            precursor_mz=400.0, precursor_charge=2, peptidoform="PEP",
            base_sequence="PEP", source_file="ref.mzML", is_reference=True,
        )

        query2 = LoadedSpectrum(
            spectrum_id="scan=201", mz=np.array([151.0]), intensity=np.array([501.0]),
            precursor_mz=400.0, precursor_charge=2, peptidoform="PEP",
            base_sequence="PEP", source_file="query.mzML", is_reference=False,
        )

        query3 = LoadedSpectrum(
            spectrum_id="scan=202", mz=np.array([152.0]), intensity=np.array([502.0]),
            precursor_mz=400.0, precursor_charge=2, peptidoform="PEP",
            base_sequence="PEP", source_file="query2.mzML", is_reference=False,
        )

        save_spectra_to_parquet(
            reference_spectra={"PEP|": [ref, ref2]},
            query_spectra=[("PEP|", query), ("PEP|", query2), ("PEP|", query3)],
            output_path=output_path,
        )

        refs, queries = get_parquet_row_counts(output_path)
        assert refs == 2
        assert queries == 3

    def test_nonexistent_file(self):
        """Returns (0, 0) for nonexistent file."""
        refs, queries = get_parquet_row_counts("/nonexistent/path.parquet")
        assert refs == 0
        assert queries == 0
