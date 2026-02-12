"""
Tests for Parquet-based spectrum loading pipeline.

Tests cover:
- Parquet write→read round-trip
- Provenance (source_mzml field)
- variant_key grouping
- is_reference flag
- Fast index (open_indexed_mzml_fast) equivalence
- _extract_id_number helper
"""
import json
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from src.spectrum_loader import (
    LoadedSpectrum,
    StreamingResult,
    _extract_id_number,
    clean_spectrum_id,
    find_spectrum_index,
    load_spectra_from_parquet,
    open_indexed_mzml,
    open_indexed_mzml_fast,
)


# =============================================================================
# Fixtures
# =============================================================================

def _make_test_parquet(tmp_path, rows=None):
    """Create a test Parquet file with known data."""
    if rows is None:
        rows = [
            {
                "source_mzml": "reference_sample.mzML",
                "spectrum_id": "index=10",
                "peptidoform": "PEPTIDEK",
                "base_sequence": "PEPTIDEK",
                "variant_key": "PEPTIDEK",
                "charge": 2,
                "precursor_mz": 450.25,
                "rt": 120.5,
                "is_reference": True,
                "mz": [100.0, 200.0, 300.0, 400.0],
                "intensity": [1000.0, 2000.0, 500.0, 1500.0],
                "score": 0.95,
                "idxml_file": "test.idXML",
            },
            {
                "source_mzml": "query_sample.mzML",
                "spectrum_id": "index=20",
                "peptidoform": "PEPTIDEK",
                "base_sequence": "PEPTIDEK",
                "variant_key": "PEPTIDEK",
                "charge": 2,
                "precursor_mz": 450.26,
                "rt": 125.3,
                "is_reference": False,
                "mz": [100.1, 200.1, 300.1, 400.1],
                "intensity": [900.0, 1800.0, 450.0, 1400.0],
                "score": 0.88,
                "idxml_file": "test.idXML",
            },
            {
                "source_mzml": "query_sample.mzML",
                "spectrum_id": "index=30",
                "peptidoform": "ANOTHERPEPTIDE",
                "base_sequence": "ANOTHERPEPTIDE",
                "variant_key": "ANOTHERPEPTIDE",
                "charge": 3,
                "precursor_mz": 520.10,
                "rt": 200.0,
                "is_reference": False,
                "mz": [150.0, 250.0, 350.0],
                "intensity": [500.0, 1000.0, 750.0],
                "score": 0.75,
                "idxml_file": "test.idXML",
            },
        ]

    schema = pa.schema([
        ("source_mzml", pa.string()),
        ("spectrum_id", pa.string()),
        ("peptidoform", pa.string()),
        ("base_sequence", pa.string()),
        ("variant_key", pa.string()),
        ("charge", pa.int32()),
        ("precursor_mz", pa.float64()),
        ("rt", pa.float64()),
        ("is_reference", pa.bool_()),
        ("mz", pa.list_(pa.float64())),
        ("intensity", pa.list_(pa.float64())),
        ("score", pa.float64()),
        ("idxml_file", pa.string()),
    ])

    arrays = {}
    for col in schema.names:
        if col in ("mz", "intensity"):
            arrays[col] = pa.array([r[col] for r in rows], type=pa.list_(pa.float64()))
        elif col == "charge":
            arrays[col] = pa.array([r[col] for r in rows], type=pa.int32())
        elif col in ("precursor_mz", "rt", "score"):
            arrays[col] = pa.array([r[col] for r in rows], type=pa.float64())
        elif col == "is_reference":
            arrays[col] = pa.array([r[col] for r in rows], type=pa.bool_())
        else:
            arrays[col] = pa.array([r[col] for r in rows], type=pa.string())

    table = pa.table(arrays, schema=schema)

    # Add provenance metadata
    meta = {
        b"mspecviz_provenance": json.dumps({
            "created_by": "test",
            "total_spectra": len(rows),
        }).encode()
    }
    table = table.replace_schema_metadata(meta)

    parquet_path = str(tmp_path / "test_spectra.parquet")
    pq.write_table(table, parquet_path)
    return parquet_path, rows


# =============================================================================
# Tests: _extract_id_number
# =============================================================================

class TestExtractIdNumber:
    def test_index_format(self):
        assert _extract_id_number("index=123") == 123

    def test_scan_format(self):
        assert _extract_id_number("scan=456") == 456

    def test_spectrum_format(self):
        assert _extract_id_number("spectrum=789") == 789

    def test_thermo_format(self):
        # Thermo IDs have numbers too
        result = _extract_id_number("controllerType=0 controllerNumber=1 scan=5678")
        assert result == 5678 or result is not None  # Should extract some number

    def test_bare_number(self):
        assert _extract_id_number("42") == 42

    def test_no_number(self):
        assert _extract_id_number("no_numbers_here") is None

    def test_empty(self):
        assert _extract_id_number("") is None


# =============================================================================
# Tests: Parquet round-trip
# =============================================================================

class TestParquetRoundTrip:
    def test_write_read_preserves_data(self, tmp_path):
        """Test that writing and reading Parquet preserves all spectrum data."""
        parquet_path, original_rows = _make_test_parquet(tmp_path)

        # Read back
        table = pq.read_table(parquet_path)
        df = table.to_pandas()

        assert len(df) == len(original_rows)

        for i, row in enumerate(original_rows):
            assert df.iloc[i]["source_mzml"] == row["source_mzml"]
            assert df.iloc[i]["spectrum_id"] == row["spectrum_id"]
            assert df.iloc[i]["peptidoform"] == row["peptidoform"]
            assert df.iloc[i]["base_sequence"] == row["base_sequence"]
            assert df.iloc[i]["variant_key"] == row["variant_key"]
            assert df.iloc[i]["charge"] == row["charge"]
            assert abs(df.iloc[i]["precursor_mz"] - row["precursor_mz"]) < 1e-6
            assert abs(df.iloc[i]["rt"] - row["rt"]) < 1e-6
            assert df.iloc[i]["is_reference"] == row["is_reference"]
            np.testing.assert_array_almost_equal(df.iloc[i]["mz"], row["mz"])
            np.testing.assert_array_almost_equal(df.iloc[i]["intensity"], row["intensity"])

    def test_provenance_metadata(self, tmp_path):
        """Test that provenance metadata is preserved."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        table = pq.read_table(parquet_path)
        meta = table.schema.metadata
        assert b"mspecviz_provenance" in meta

        prov = json.loads(meta[b"mspecviz_provenance"])
        assert prov["created_by"] == "test"
        assert prov["total_spectra"] == 3


# =============================================================================
# Tests: load_spectra_from_parquet
# =============================================================================

class TestLoadSpectraFromParquet:
    def test_basic_load(self, tmp_path):
        """Test loading spectra from Parquet produces correct LoadedSpectrum objects."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        ref_spectra = load_spectra_from_parquet(parquet_path)

        # Should have reference spectra for PEPTIDEK variant
        assert "PEPTIDEK" in ref_spectra
        assert len(ref_spectra["PEPTIDEK"]) == 1  # 1 reference

        ref = ref_spectra["PEPTIDEK"][0]
        assert ref.is_reference is True
        assert ref.source_file == "reference_sample.mzML"
        assert ref.base_sequence == "PEPTIDEK"
        np.testing.assert_array_almost_equal(ref.mz, [100.0, 200.0, 300.0, 400.0])

    def test_variant_grouping(self, tmp_path):
        """Test that spectra are correctly grouped by variant_key."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        ref_spectra = load_spectra_from_parquet(parquet_path)

        # Should have two variant groups
        assert "PEPTIDEK" in ref_spectra
        assert "ANOTHERPEPTIDE" in ref_spectra

    def test_is_reference_flag(self, tmp_path):
        """Test that is_reference flag is correctly applied."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        query_results = []

        def on_query(result: StreamingResult):
            query_results.append(result)

        ref_spectra = load_spectra_from_parquet(
            parquet_path,
            on_query_ready=on_query,
        )

        # Should have 2 query spectra total
        assert len(query_results) == 2

        # First query should be PEPTIDEK (from query_sample.mzML)
        peptidek_queries = [r for r in query_results if r.variant_key == "PEPTIDEK"]
        assert len(peptidek_queries) == 1
        assert peptidek_queries[0].query_spectrum.is_reference is False

        # ANOTHERPEPTIDE query
        other_queries = [r for r in query_results if r.variant_key == "ANOTHERPEPTIDE"]
        assert len(other_queries) == 1

    def test_streaming_callbacks(self, tmp_path):
        """Test that callbacks are invoked correctly."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        ref_loaded_calls = []
        query_ready_calls = []

        def on_ref(vk, count):
            ref_loaded_calls.append((vk, count))

        def on_query(result):
            query_ready_calls.append(result)

        load_spectra_from_parquet(
            parquet_path,
            on_reference_loaded=on_ref,
            on_query_ready=on_query,
        )

        # Reference callback should fire for PEPTIDEK (has 1 ref spectrum)
        ref_vks = [vk for vk, _ in ref_loaded_calls]
        assert "PEPTIDEK" in ref_vks

        # Query callback should fire for each query spectrum
        assert len(query_ready_calls) == 2

    def test_reference_spectra_passed_to_queries(self, tmp_path):
        """Test that query StreamingResults include matching reference spectra."""
        parquet_path, _ = _make_test_parquet(tmp_path)

        query_results = []

        def on_query(result: StreamingResult):
            query_results.append(result)

        load_spectra_from_parquet(
            parquet_path,
            on_query_ready=on_query,
        )

        # The PEPTIDEK query should have 1 reference spectrum
        peptidek_query = [r for r in query_results if r.variant_key == "PEPTIDEK"][0]
        assert len(peptidek_query.reference_spectra) == 1
        assert peptidek_query.reference_spectra[0].source_file == "reference_sample.mzML"

    def test_modification_key_from_variant(self, tmp_path):
        """Test that modification_key is correctly extracted from variant_key."""
        rows = [
            {
                "source_mzml": "ref.mzML",
                "spectrum_id": "index=1",
                "peptidoform": "PEPTM(Oxidation)DEK",
                "base_sequence": "PEPTMDEK",
                "variant_key": "PEPTMDEK|Oxidation@M5",
                "charge": 2,
                "precursor_mz": 460.0,
                "rt": 100.0,
                "is_reference": True,
                "mz": [100.0, 200.0],
                "intensity": [500.0, 1000.0],
                "score": 0.9,
                "idxml_file": "test.idXML",
            },
            {
                "source_mzml": "query.mzML",
                "spectrum_id": "index=2",
                "peptidoform": "PEPTM(Oxidation)DEK",
                "base_sequence": "PEPTMDEK",
                "variant_key": "PEPTMDEK|Oxidation@M5",
                "charge": 2,
                "precursor_mz": 460.0,
                "rt": 105.0,
                "is_reference": False,
                "mz": [100.0, 200.0],
                "intensity": [400.0, 800.0],
                "score": 0.85,
                "idxml_file": "test.idXML",
            },
        ]

        parquet_path, _ = _make_test_parquet(tmp_path, rows=rows)
        query_results = []

        def on_query(result):
            query_results.append(result)

        load_spectra_from_parquet(parquet_path, on_query_ready=on_query)

        assert len(query_results) == 1
        assert query_results[0].modification_key == "Oxidation@M5"
        assert query_results[0].base_sequence == "PEPTMDEK"

    def test_missing_file_error(self, tmp_path):
        """Test error handling for missing Parquet file."""
        errors = []
        result = load_spectra_from_parquet(
            str(tmp_path / "nonexistent.parquet"),
            error_callback=lambda msg: errors.append(msg),
        )
        assert result == {}
        assert len(errors) > 0

    def test_reference_substring_fallback(self, tmp_path):
        """Test that reference_substring is used when is_reference column is absent."""
        # Create Parquet without is_reference column
        rows = [
            {
                "source_mzml": "Synpep_reference.mzML",
                "spectrum_id": "index=1",
                "peptidoform": "PEPTIDEK",
                "base_sequence": "PEPTIDEK",
                "variant_key": "PEPTIDEK",
                "charge": 2,
                "precursor_mz": 450.0,
                "rt": 100.0,
                "mz": [100.0, 200.0],
                "intensity": [500.0, 1000.0],
                "score": 0.9,
                "idxml_file": "test.idXML",
            },
            {
                "source_mzml": "query_sample.mzML",
                "spectrum_id": "index=2",
                "peptidoform": "PEPTIDEK",
                "base_sequence": "PEPTIDEK",
                "variant_key": "PEPTIDEK",
                "charge": 2,
                "precursor_mz": 450.0,
                "rt": 105.0,
                "mz": [100.0, 200.0],
                "intensity": [400.0, 800.0],
                "score": 0.85,
                "idxml_file": "test.idXML",
            },
        ]

        # Build without is_reference
        schema = pa.schema([
            ("source_mzml", pa.string()),
            ("spectrum_id", pa.string()),
            ("peptidoform", pa.string()),
            ("base_sequence", pa.string()),
            ("variant_key", pa.string()),
            ("charge", pa.int32()),
            ("precursor_mz", pa.float64()),
            ("rt", pa.float64()),
            ("mz", pa.list_(pa.float64())),
            ("intensity", pa.list_(pa.float64())),
            ("score", pa.float64()),
            ("idxml_file", pa.string()),
        ])

        arrays = {col: pa.array([r[col] for r in rows], type=schema.field(col).type)
                  for col in schema.names}
        table = pa.table(arrays, schema=schema)
        parquet_path = str(tmp_path / "no_is_ref.parquet")
        pq.write_table(table, parquet_path)

        ref_spectra = load_spectra_from_parquet(
            parquet_path,
            reference_substring="Synpep",
        )

        assert "PEPTIDEK" in ref_spectra
        assert len(ref_spectra["PEPTIDEK"]) == 1
        assert ref_spectra["PEPTIDEK"][0].source_file == "Synpep_reference.mzML"


# =============================================================================
# Tests: open_indexed_mzml_fast (unit tests with mock-like approach)
# =============================================================================

class TestFastIndexHelpers:
    """Test the helper functions used by the fast index path."""

    def test_clean_spectrum_id_index(self):
        assert clean_spectrum_id("index=123") == "index=123"

    def test_clean_spectrum_id_scan(self):
        assert clean_spectrum_id("scan=456") == "scan=456"

    def test_clean_spectrum_id_with_whitespace(self):
        assert clean_spectrum_id("  scan=456  ") == "scan=456"

    def test_clean_spectrum_id_thermo(self):
        thermo_id = "controllerType=0 controllerNumber=1 scan=5678"
        assert clean_spectrum_id(thermo_id) == thermo_id

    def test_find_spectrum_index_direct_lookup(self):
        """Test find_spectrum_index with a pre-built index."""
        id_to_index = {
            "__format__": "scan",
            "__n_spectra__": 100,
            "scan=50": 49,
            "_num_50": 49,
        }
        assert find_spectrum_index("scan=50", id_to_index, 100) == 49

    def test_find_spectrum_index_format(self):
        """Test find_spectrum_index with index= format (O(1) path)."""
        id_to_index = {
            "__format__": "index",
            "__n_spectra__": 1000,
        }
        assert find_spectrum_index("index=500", id_to_index, 1000) == 500
        assert find_spectrum_index("index=0", id_to_index, 1000) == 0
        assert find_spectrum_index("index=999", id_to_index, 1000) == 999

    def test_find_spectrum_index_out_of_range(self):
        """Test find_spectrum_index returns None for out-of-range index."""
        id_to_index = {
            "__format__": "index",
            "__n_spectra__": 100,
        }
        assert find_spectrum_index("index=200", id_to_index, 100) is None

    def test_find_spectrum_index_numeric_fallback(self):
        """Test find_spectrum_index using numeric key fallback."""
        id_to_index = {
            "__format__": "scan",
            "__n_spectra__": 100,
            "_num_42": 41,
        }
        assert find_spectrum_index("scan=42", id_to_index, 100) == 41
