"""
Unit tests for src/insight_dataframes.py.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pytest

from src.insight_dataframes import (
    PSM_TABLE_SCHEMA,
    PEAKS_SCHEMA,
    SEQUENCE_SCHEMA,
    HEATMAP_SCHEMA,
    validate_dataframe,
    _proforma_to_openms_str,
    build_psm_table,
    build_peaks_dataframe,
    build_sequence_dataframe,
    build_heatmap_dataframe,
    save_to_parquet,
)
from src.peptide_search import VariantGroup
from src.spectrum_loader import LoadedSpectrum


from tests.conftest import make_peptide_match


def _make_variant_groups():
    """Build a simple dict of variant groups for testing."""
    vg = VariantGroup(
        base_sequence="PEPTIDEK",
        modification_key="",
        display_label="PEPTIDEK",
        matches=[
            make_peptide_match(spectrum_id="index=1"),
            make_peptide_match(spectrum_id="index=2"),
        ],
    )
    return {"PEPTIDEK": vg}


def _make_annotated_df(n_peaks=5):
    """Build a mock annotated DataFrame as returned by prepare_spectrum_for_plot."""
    return pd.DataFrame({
        "mz": np.linspace(100.0, 600.0, n_peaks).astype(np.float64),
        "intensity": np.linspace(0.2, 1.0, n_peaks).astype(np.float64),
        "annotation": ["b1", "y2", "b3", "y4", "unannotated"][:n_peaks],
        "is_annotated": [True, True, True, True, False][:n_peaks],
        "ion_type": ["b", "y", "b", "y", "unannotated"][:n_peaks],
    })


def _make_loaded_spectrum(spectrum_id="index=1", peptidoform="PEPTIDEK",
                          charge=2, rt=120.0, n_peaks=5):
    """Build a mock LoadedSpectrum."""
    return LoadedSpectrum(
        spectrum_id=spectrum_id,
        mz=np.linspace(100.0, 600.0, n_peaks).astype(np.float64),
        intensity=np.linspace(100.0, 1000.0, n_peaks).astype(np.float64),
        precursor_mz=500.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence="PEPTIDEK",
        source_file="test.mzML",
        rt=rt,
    )


# ---------------------------------------------------------------------------
# TestValidateDataframe
# ---------------------------------------------------------------------------

class TestValidateDataframe:
    @pytest.mark.parametrize("columns, expected", [
        # Valid schema
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32)}, True),
        # Missing columns
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32)}, False),
        # Wrong type (Float64 instead of Float32)
        ({"x": pl.Series([1.0], dtype=pl.Float64),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32)}, False),
        # Extra columns
        ({"x": pl.Series([1.0], dtype=pl.Float32),
          "y": pl.Series([2.0], dtype=pl.Float32),
          "highlight": pl.Series([1], dtype=pl.Int32),
          "annotation": ["b1"],
          "ion_type": ["b"],
          "psm_id": pl.Series([0], dtype=pl.Int32),
          "extra_col": [42]}, False),
    ], ids=["valid", "missing_columns", "wrong_type", "extra_columns"])
    def test_schema_validation(self, columns, expected):
        df = pl.DataFrame(columns)
        assert validate_dataframe(df.lazy(), PEAKS_SCHEMA, "test") is expected

    def test_empty_dataframe_valid(self):
        df = pl.DataFrame(schema=PEAKS_SCHEMA)
        assert validate_dataframe(df.lazy(), PEAKS_SCHEMA, "test") is True


# ---------------------------------------------------------------------------
# TestProformaToOpenmsStr
# ---------------------------------------------------------------------------

class TestProformaToOpenmsStr:
    def test_simple_peptide(self):
        result = _proforma_to_openms_str("PEPTIDEK")
        assert result == "PEPTIDEK"

    def test_oxidation_modification(self):
        result = _proforma_to_openms_str("PEPTM[Oxidation]IDEK")
        assert "Oxidation" in result
        assert "(" in result  # OpenMS uses parentheses

    def test_charge_suffix_stripped(self):
        result = _proforma_to_openms_str("PEPTIDEK/2")
        assert "/2" not in result
        assert "PEPTIDEK" in result

    def test_nterm_modification(self):
        result = _proforma_to_openms_str("[Acetyl]-PEPTIDEK")
        assert "(Acetyl)" in result

    def test_fallback_on_invalid(self):
        # Should not raise, returns best-effort string
        result = _proforma_to_openms_str("PEPTIDEK")
        assert isinstance(result, str)
        assert len(result) > 0


# ---------------------------------------------------------------------------
# TestBuildPsmTable
# ---------------------------------------------------------------------------

class TestBuildPsmTable:
    def test_basic_build(self):
        """Build PSM table from standard variant groups: type, schema, row count, lookup keys."""
        vgs = _make_variant_groups()
        lf, lookup = build_psm_table(vgs)
        assert isinstance(lf, pl.LazyFrame)
        assert isinstance(lookup, dict)
        assert validate_dataframe(lf, PSM_TABLE_SCHEMA, "PSM table")
        df = lf.collect()
        assert len(df) == 2
        assert set(lookup.keys()) == {"index=1|PEPTIDEK", "index=2|PEPTIDEK"}

    def test_empty_variant_groups(self):
        lf, lookup = build_psm_table({})
        assert len(lf.collect()) == 0
        assert len(lookup) == 0

    def test_multiple_variants(self):
        vg1 = VariantGroup(
            base_sequence="PEPTIDEK", modification_key="",
            display_label="PEPTIDEK",
            matches=[make_peptide_match(spectrum_id="index=1")],
        )
        vg2 = VariantGroup(
            base_sequence="AAAVAAAK", modification_key="",
            display_label="AAAVAAAK",
            matches=[make_peptide_match(peptidoform="AAAVAAAK", base_seq="AAAVAAAK",
                                 spectrum_id="index=3")],
        )
        lf, lookup = build_psm_table({"PEPTIDEK": vg1, "AAAVAAAK": vg2})
        assert len(lf.collect()) == 2
        assert len(lookup) == 2


# ---------------------------------------------------------------------------
# TestBuildPeaksDataframe
# ---------------------------------------------------------------------------

class TestBuildPeaksDataframe:
    def test_basic_build(self):
        """Build peaks DF: type, schema, peak count, highlight column."""
        annotated = [(0, _make_annotated_df(n_peaks=5))]
        lf = build_peaks_dataframe(annotated)
        assert isinstance(lf, pl.LazyFrame)
        assert validate_dataframe(lf, PEAKS_SCHEMA, "Peaks")
        df = lf.collect()
        assert len(df) == 5
        # First 4 peaks are annotated (highlight=1), last is not (highlight=0)
        assert df["highlight"][0] == 1
        assert df["highlight"][4] == 0

    def test_empty_input(self):
        assert len(build_peaks_dataframe([]).collect()) == 0

    def test_multiple_spectra(self):
        annotated = [
            (0, _make_annotated_df(n_peaks=3)),
            (1, _make_annotated_df(n_peaks=4)),
        ]
        df = build_peaks_dataframe(annotated).collect()
        assert len(df) == 7
        assert set(df["psm_id"].to_list()) == {0, 1}


# ---------------------------------------------------------------------------
# TestBuildSequenceDataframe
# ---------------------------------------------------------------------------

class TestBuildSequenceDataframe:
    def test_basic_build(self):
        """Build sequence DF: type, schema, row count."""
        spectra = [
            (0, _make_loaded_spectrum(spectrum_id="index=1")),
            (1, _make_loaded_spectrum(spectrum_id="index=2")),
        ]
        lf = build_sequence_dataframe(spectra)
        assert isinstance(lf, pl.LazyFrame)
        assert validate_dataframe(lf, SEQUENCE_SCHEMA, "Sequence")
        assert len(lf.collect()) == 2

    def test_empty_input(self):
        assert len(build_sequence_dataframe([]).collect()) == 0

    def test_sequence_is_openms_notation(self):
        spectra = [(0, _make_loaded_spectrum(peptidoform="PEPTM[Oxidation]IDEK"))]
        seq = build_sequence_dataframe(spectra).collect()["sequence"][0]
        assert "(" in seq or "PEPTM" in seq


# ---------------------------------------------------------------------------
# TestBuildHeatmapDataframe
# ---------------------------------------------------------------------------

class TestBuildHeatmapDataframe:
    def test_basic_build(self):
        """Build heatmap DF: type, schema, point count."""
        spectra = [_make_loaded_spectrum(n_peaks=10)]
        lf = build_heatmap_dataframe(spectra)
        assert isinstance(lf, pl.LazyFrame)
        assert validate_dataframe(lf, HEATMAP_SCHEMA, "Heatmap")
        assert len(lf.collect()) == 10

    def test_empty_input(self):
        assert len(build_heatmap_dataframe([]).collect()) == 0

    def test_min_rt_filtering(self):
        spec_low_rt = _make_loaded_spectrum(rt=10.0, n_peaks=5)
        spec_high_rt = _make_loaded_spectrum(rt=200.0, n_peaks=5)
        df = build_heatmap_dataframe([spec_low_rt, spec_high_rt], min_rt=100.0).collect()
        assert len(df) == 5
        assert all(df["x"].to_list()[i] == pytest.approx(200.0, abs=0.1) for i in range(5))

    def test_multiple_spectra(self):
        spectra = [
            _make_loaded_spectrum(rt=100.0, n_peaks=3),
            _make_loaded_spectrum(rt=200.0, n_peaks=4),
        ]
        assert len(build_heatmap_dataframe(spectra).collect()) == 7


# ---------------------------------------------------------------------------
# TestSaveToParquet (round-trip)
# ---------------------------------------------------------------------------

class TestSaveToParquet:
    @pytest.mark.parametrize("build_fn", ["peaks", "psm"], ids=["peaks", "psm_table"])
    def test_round_trip(self, tmp_path, build_fn):
        if build_fn == "peaks":
            lf = build_peaks_dataframe([(0, _make_annotated_df())])
        else:
            lf, _ = build_psm_table(_make_variant_groups())
        parquet_path = str(tmp_path / f"{build_fn}.parquet")
        save_to_parquet(lf, parquet_path)
        assert Path(parquet_path).exists()
        read_back = pl.read_parquet(parquet_path)
        original = lf.collect()
        assert read_back.shape == original.shape
        assert read_back.columns == original.columns

    def test_creates_parent_directories(self, tmp_path):
        lf = pl.DataFrame(schema=PEAKS_SCHEMA).lazy()
        parquet_path = str(tmp_path / "subdir" / "nested" / "test.parquet")
        save_to_parquet(lf, parquet_path)
        assert Path(parquet_path).exists()


# ---------------------------------------------------------------------------
# Integration tests using real test data
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestBuildFromRealData:
    """Build DataFrames from real search results."""

    @pytest.fixture
    def variant_groups(self, idxml_path, target_peptides):
        from src.peptide_search import search_peptides_in_idxmls
        return search_peptides_in_idxmls([idxml_path], target_peptides)

    def test_build_psm_table_real(self, variant_groups):
        psm_lf, lookup = build_psm_table(variant_groups)
        assert validate_dataframe(psm_lf, PSM_TABLE_SCHEMA, "psm")
        df = psm_lf.collect()
        assert len(df) == 2
        assert len(lookup) == 2

    def test_build_all_dataframes_real(self, variant_groups, mzml_paths):
        """Full pipeline: search -> stream -> annotate -> build DataFrames."""
        from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key
        from src.peak_annotation import prepare_spectrum_for_plot

        psm_lf, psm_lookup = build_psm_table(variant_groups)
        annotated = []

        def on_query(result):
            key = spectrum_lookup_key(result.query_spectrum)
            pid = psm_lookup.get(key)
            if pid is None:
                return
            df = prepare_spectrum_for_plot(
                result.query_spectrum.mz, result.query_spectrum.intensity,
                result.query_spectrum.peptidoform, result.query_spectrum.precursor_charge,
                ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
            )
            annotated.append((pid, result.query_spectrum, df))

        mzml_path_list = [Path(p) for p in mzml_paths]
        load_spectra_streaming(
            variant_groups=variant_groups,
            mzml_paths=mzml_path_list,
            reference_substring="",
            on_query_ready=on_query,
        )

        assert len(annotated) == 2

        peaks_lf = build_peaks_dataframe([(pid, df) for pid, _, df in annotated])
        seq_lf = build_sequence_dataframe([(pid, s) for pid, s, _ in annotated])
        heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated])

        assert validate_dataframe(peaks_lf, PEAKS_SCHEMA, "peaks")
        assert validate_dataframe(seq_lf, SEQUENCE_SCHEMA, "seq")
        assert validate_dataframe(heatmap_lf, HEATMAP_SCHEMA, "heatmap")

        peaks_df = peaks_lf.collect()
        assert len(peaks_df) == 373
        assert peaks_df["psm_id"].n_unique() == 2
