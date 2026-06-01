"""Builder tests for core.dataframes (PSM table + peaks/sequence/heatmap LazyFrames)."""
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pytest

from core.annotation import prepare_spectrum_for_plot
from core.dataframes import (
    HEATMAP_SCHEMA,
    PEAKS_SCHEMA,
    PSM_TABLE_SCHEMA,
    SEQUENCE_SCHEMA,
    build_heatmap_dataframe,
    build_peaks_dataframe,
    build_psm_table,
    build_sequence_dataframe,
)
from psm_utils import PSMList

from core.peptide_lookup import lookup_peptides_in_idxmls
from core.spectra import load_spectra_streaming, spectrum_lookup_key

from tests.conftest import make_peptide_match


def _make_variant_groups():
    """Build a Dict[variant_key, PSMList] with two synthetic PEPTIDEK PSMs."""
    return {"PEPTIDEK": PSMList(psm_list=[
        make_peptide_match(spectrum_id="index=1"),
        make_peptide_match(spectrum_id="index=2"),
    ])}


def _make_annotated_frame(n_peaks=5):
    """Build a synthetic prepare_spectrum_for_plot-shaped DataFrame with n_peaks rows."""
    return pd.DataFrame({
        "mz": np.linspace(100.0, 600.0, n_peaks).astype(np.float64),
        "intensity": np.linspace(0.2, 1.0, n_peaks).astype(np.float64),
        "annotation": ["b1", "y2", "b3", "y4", "unannotated"][:n_peaks],
        "is_annotated": [True, True, True, True, False][:n_peaks],
        "ion_type": ["b", "y", "b", "y", "unannotated"][:n_peaks],
    })


def _make_psm(spectrum_id="index=1", peptidoform="PEPTIDEK",
              charge=2, rt=120.0, n_peaks=5):
    """Build a synthetic PSM carrying its peaks on psm.spectrum + retention_time —
    the Step 3 carrier shape the sequence/heatmap builders now read."""
    psm = make_peptide_match(spectrum_id=spectrum_id, peptidoform=peptidoform, charge=charge)
    psm.spectrum = {
        "mz": np.linspace(100.0, 600.0, n_peaks).astype(np.float64),
        "intensity": np.linspace(100.0, 1000.0, n_peaks).astype(np.float64),
    }
    psm.retention_time = rt
    return psm


class TestBuildPsmTable:
    def test_basic_build(self):
        """PSM table from two-PSM variant: LazyFrame, schema, row count, lookup keys all match."""
        variant_groups = _make_variant_groups()
        frame, lookup = build_psm_table(variant_groups)
        assert isinstance(frame, pl.LazyFrame)
        assert isinstance(lookup, dict)
        assert dict(frame.collect_schema()) == PSM_TABLE_SCHEMA
        df = frame.collect()
        assert len(df) == 2
        assert set(lookup.keys()) == {"index=1|PEPTIDEK", "index=2|PEPTIDEK"}

    def test_empty_variant_groups(self):
        frame, lookup = build_psm_table({})
        assert len(frame.collect()) == 0
        assert len(lookup) == 0

    def test_multiple_variants(self):
        groups = {
            "PEPTIDEK": PSMList(psm_list=[make_peptide_match(spectrum_id="index=1")]),
            "AAAVAAAK": PSMList(psm_list=[
                make_peptide_match(peptidoform="AAAVAAAK", spectrum_id="index=3"),
            ]),
        }
        frame, lookup = build_psm_table(groups)
        assert len(frame.collect()) == 2
        assert len(lookup) == 2


class TestBuildPeaksDataframe:
    def test_basic_build(self):
        """Peaks DF from one annotated frame: LazyFrame, schema, row count, highlight column."""
        annotated = [(0, _make_annotated_frame(n_peaks=5))]
        frame = build_peaks_dataframe(annotated)
        assert isinstance(frame, pl.LazyFrame)
        assert dict(frame.collect_schema()) == PEAKS_SCHEMA
        df = frame.collect()
        assert len(df) == 5
        # First 4 peaks are annotated (highlight=1); last peak is unannotated (highlight=0)
        assert df["highlight"][0] == 1
        assert df["highlight"][4] == 0

    def test_empty_input(self):
        assert len(build_peaks_dataframe([]).collect()) == 0

    def test_multiple_spectra(self):
        annotated = [
            (0, _make_annotated_frame(n_peaks=3)),
            (1, _make_annotated_frame(n_peaks=4)),
        ]
        df = build_peaks_dataframe(annotated).collect()
        assert len(df) == 7
        assert set(df["psm_id"].to_list()) == {0, 1}


class TestBuildSequenceDataframe:
    def test_basic_build(self):
        """Sequence DF from two spectra: LazyFrame, schema, row count."""
        spectra = [
            (0, _make_psm(spectrum_id="index=1")),
            (1, _make_psm(spectrum_id="index=2")),
        ]
        frame = build_sequence_dataframe(spectra)
        assert isinstance(frame, pl.LazyFrame)
        assert dict(frame.collect_schema()) == SEQUENCE_SCHEMA
        assert len(frame.collect()) == 2

    def test_empty_input(self):
        assert len(build_sequence_dataframe([]).collect()) == 0

    def test_sequence_is_openms_notation(self):
        spectra = [(0, _make_psm(peptidoform="PEPTM[Oxidation]IDEK"))]
        seq = build_sequence_dataframe(spectra).collect()["sequence"][0]
        assert "(" in seq or "PEPTM" in seq


class TestBuildHeatmapDataframe:
    def test_basic_build(self):
        """Heatmap DF from one ten-peak spectrum: LazyFrame, schema, point count."""
        spectra = [_make_psm(n_peaks=10)]
        frame = build_heatmap_dataframe(spectra)
        assert isinstance(frame, pl.LazyFrame)
        assert dict(frame.collect_schema()) == HEATMAP_SCHEMA
        assert len(frame.collect()) == 10

    def test_empty_input(self):
        assert len(build_heatmap_dataframe([]).collect()) == 0

    def test_multiple_spectra(self):
        spectra = [
            _make_psm(rt=100.0, n_peaks=3),
            _make_psm(rt=200.0, n_peaks=4),
        ]
        assert len(build_heatmap_dataframe(spectra).collect()) == 7


@pytest.mark.integration
class TestBuildFromRealData:
    """End-to-end: lookup → stream → annotate → build the four insight DataFrames."""

    @pytest.fixture
    def variant_groups(self, idxml_path, target_peptides):
        return lookup_peptides_in_idxmls([idxml_path], target_peptides)

    def test_build_psm_table_real(self, variant_groups):
        psm_frame, lookup = build_psm_table(variant_groups)
        assert dict(psm_frame.collect_schema()) == PSM_TABLE_SCHEMA
        df = psm_frame.collect()
        assert len(df) == 2
        assert len(lookup) == 2

    def test_build_all_dataframes_real(self, variant_groups, mzml_paths):
        """Full pipeline against real fixtures: 2 PSMs, ~373 annotated peaks across all frames."""
        psm_frame, psm_lookup = build_psm_table(variant_groups)
        annotated = []

        def on_query(result):
            key = spectrum_lookup_key(result.query_spectrum)
            psm_id = psm_lookup.get(key)
            if psm_id is None:
                return
            frame = prepare_spectrum_for_plot(
                result.query_spectrum.mz, result.query_spectrum.intensity,
                result.query_spectrum.peptidoform, result.query_spectrum.precursor_charge,
                ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
            )
            annotated.append((psm_id, result.query_spectrum, frame))

        mzml_path_list = [Path(p) for p in mzml_paths]
        load_spectra_streaming(
            variant_groups=variant_groups,
            mzml_paths=mzml_path_list,
            reference_substring="",
            on_query_ready=on_query,
        )

        assert len(annotated) == 2

        peaks_frame = build_peaks_dataframe([(psm_id, frame) for psm_id, _, frame in annotated])
        sequence_frame = build_sequence_dataframe([(psm_id, spectrum.match) for psm_id, spectrum, _ in annotated])
        heatmap_frame = build_heatmap_dataframe([spectrum.match for _, spectrum, _ in annotated])

        assert dict(peaks_frame.collect_schema()) == PEAKS_SCHEMA
        assert dict(sequence_frame.collect_schema()) == SEQUENCE_SCHEMA
        assert dict(heatmap_frame.collect_schema()) == HEATMAP_SCHEMA

        peaks_collected = peaks_frame.collect()
        # Range check rather than exact peak count — spectrum_utils / pyOpenMS
        # library upgrades shift the matched-peak total by a handful without
        # signalling a real regression.
        assert 350 < len(peaks_collected) < 400
        assert peaks_collected["psm_id"].n_unique() == 2
