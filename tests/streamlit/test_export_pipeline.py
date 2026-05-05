"""Unit tests for ``streamlit_app.views.export.pipeline._run_export``.

Drives ``_run_export`` with a mocked ``st`` module + session_state to verify
each format flag creates the expected artifacts inside the versioned
``export_*/`` directory. ``export_controls_fragment`` bug-fix checks live in
``test_export_bugfixes.py``.
"""
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import polars as pl
import pytest

from mspecviz.annotation import METRIC_SHORT_TO_DISPLAY
from mspecviz.search import PeptideMatch, VariantGroup
from mspecviz.spectra import LoadedSpectrum
from streamlit_app.views.export.pipeline import _run_export


def _make_loaded_spectrum(peptidoform="PEPTIDEK/2", base_seq="PEPTIDEK",
                          source_file="query.mzML", is_ref=False):
    return LoadedSpectrum(
        spectrum_id="index=1",
        mz=np.array([100.0, 200.0, 300.0, 400.0]),
        intensity=np.array([1000.0, 2000.0, 3000.0, 500.0]),
        precursor_mz=500.0,
        precursor_charge=2,
        peptidoform=peptidoform,
        base_sequence=base_seq,
        source_file=source_file,
        is_reference=is_ref,
        rt=120.0,
    )


def _make_annotated_frame():
    return pd.DataFrame({
        "mz": [100.0, 200.0, 300.0, 400.0],
        "intensity": [0.33, 0.67, 1.0, 0.17],
        "annotation": ["b1", "y2", "b3", "unannotated"],
        "is_annotated": [True, True, True, False],
        "ion_type": ["b", "y", "b", "unannotated"],
    })


def _make_comparison_frame():
    """Build a minimal comparison DataFrame matching the real schema."""
    rows = []
    for i, pep in enumerate(["PEPTIDEK", "ANOTHERPEP"]):
        row = {
            "psm_id": i,
            "Peptide": pep,
            "Modification": "(unmodified)",
            "Charge": 2,
            "Query": "query_file",
            "Reference": "ref_file",
            "query_spectrum_id": f"index={i}",
            "ref_spectrum_id": f"index={i+10}",
            "query_rt": 120.0 + i * 10,
            "ref_rt": 121.0 + i * 10,
        }
        for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
            row[display_name] = 0.85 + i * 0.05
        for count_key in ("total_q", "total_r", "ann_q", "ann_r", "shared"):
            row[count_key] = 5
        rows.append(row)
    return pd.DataFrame(rows)


def _make_variant_groups():
    match = PeptideMatch(
        peptidoform="PEPTIDEK/2", base_sequence="PEPTIDEK",
        spectrum_id="index=1", charge=2, score=0.95,
        idxml_file="test.idXML", mzml_filename="query.mzML",
        id_merge_index=-1,
    )
    variant_group = VariantGroup(
        base_sequence="PEPTIDEK",
        modification_key=None,
        display_label="PEPTIDEK (unmodified)",
        matches=[match],
    )
    return {"PEPTIDEK|(unmodified)": variant_group}


def _make_psm_lazyframe():
    return pl.DataFrame({
        "id": [0],
        "variant_key": ["PEPTIDEK|(unmodified)"],
        "base_sequence": ["PEPTIDEK"],
        "modification": ["(unmodified)"],
        "display_label": ["PEPTIDEK (unmodified)"],
        "peptidoform": ["PEPTIDEK/2"],
        "charge": [2],
        "score": [0.95],
        "source_file": ["query.mzML"],
        "spectrum_id": ["index=1"],
        "idxml_file": ["test.idXML"],
    }).lazy()


def _build_session_state(tmp_path):
    """Build a mock session_state dict with all keys _run_export needs."""
    query_spectrum = _make_loaded_spectrum()
    query_frame = _make_annotated_frame()
    reference_spectrum = _make_loaded_spectrum(source_file="ref.mzML", is_ref=True)
    reference_frame = _make_annotated_frame()

    return {
        "peptide_mirror.annotated_queries": {0: (query_spectrum, query_frame)},
        "peptide_mirror.reference_store": {"PEPTIDEK|(unmodified)": [(reference_spectrum, reference_frame)]},
        "peptide_mirror.variant_groups": _make_variant_groups(),
        "peptide_mirror.comparison_frame": _make_comparison_frame(),
        "peptide_mirror.psm_frame": _make_psm_lazyframe(),
        "peptide_mirror.peaks_frame": pl.DataFrame({"x": [1.0], "y": [1.0]}).lazy(),
        "peptide_mirror.sequence_frame": pl.DataFrame({"seq": ["PEP"]}).lazy(),
        "peptide_mirror.heatmap_frame": pl.DataFrame({"x": [1.0], "y": [1.0]}).lazy(),
        "peptide_mirror.data_ready": True,
        "peptide_mirror.cache_version": 1,
        "workspace": tmp_path,
        "settings": {"enable_workspaces": True},
        "export_directory": str(tmp_path / "exports"),
        "intensity_mode": "Relative",
    }


class TestRunExport:
    """Test _run_export produces expected output files."""

    @pytest.fixture
    def mock_st(self, tmp_path):
        """Patch Streamlit + helpers to let _run_export run outside the app."""
        state = _build_session_state(tmp_path)

        with patch("streamlit_app.views.export.pipeline.st") as mock_st, \
             patch("streamlit_app.views.export.pipeline.build_run_settings") as mock_brs, \
             patch("streamlit_app.views.export.pipeline.get_state_key",
                   side_effect=lambda k: f"peptide_mirror.{k}"), \
             patch("streamlit_app.views.export.pipeline.get_param") as mock_gp, \
             patch("streamlit_app.views.export.gallery.st", mock_st):

            mock_st.session_state = state
            mock_st.status.return_value.__enter__ = MagicMock(return_value=MagicMock())
            mock_st.status.return_value.__exit__ = MagicMock(return_value=False)
            mock_st.caption = MagicMock()
            mock_st.expander.return_value.__enter__ = MagicMock()
            mock_st.expander.return_value.__exit__ = MagicMock(return_value=False)
            mock_st.warning = MagicMock()

            mock_brs.return_value = {
                "peptides": ["PEPTIDEK"],
                "ion_types": ["b", "y"],
                "_created_at": "2026-01-01",
                "_last_run": "2026-01-01",
            }

            def _get_param(params, key, default=None):
                return state.get(key, default)
            mock_gp.side_effect = _get_param

            yield mock_st, tmp_path

    def _invoke_export(self, tmp_path, **flags):
        export_dir = str(tmp_path / "exports")
        defaults = dict(do_parquet=False, do_csv=False, do_mirror=False,
                        do_summary=False, do_html=False)
        defaults.update(flags)
        _run_export({}, export_dir, **defaults)
        return export_dir

    def test_parquet_export(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_parquet=True)

        export_dirs = list(Path(export_dir).iterdir())
        assert len(export_dirs) == 1
        export_root = export_dirs[0]

        assert (export_root / "run_settings.json").exists()
        parquet_files = list(export_root.glob("*.parquet"))
        assert len(parquet_files) >= 1

    def test_csv_export(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_csv=True)

        export_root = list(Path(export_dir).iterdir())[0]
        assert (export_root / "comparison_summary.csv").exists()
        frame = pd.read_csv(export_root / "comparison_summary.csv")
        assert len(frame) == 2
        assert "Peptide" in frame.columns

    def test_mirror_plot_export(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_mirror=True)

        export_root = list(Path(export_dir).iterdir())[0]
        png_files = list(export_root.rglob("*.png"))
        svg_files = list(export_root.rglob("*.svg"))
        assert len(png_files) >= 1, "Expected at least one PNG mirror plot"
        assert len(svg_files) >= 1, "Expected at least one SVG mirror plot"

    def test_summary_plot_export(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_summary=True)

        export_root = list(Path(export_dir).iterdir())[0]
        summary_dir = export_root / "summary"
        assert summary_dir.exists(), "Summary directory should be created"
        summary_files = list(summary_dir.rglob("*.*"))
        assert len(summary_files) >= 2, "Expected dual_axis + ridgeline plots"

    def test_full_export(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_parquet=True, do_csv=True,
                               do_mirror=True, do_summary=True)

        export_root = list(Path(export_dir).iterdir())[0]
        assert (export_root / "run_settings.json").exists()
        assert (export_root / "comparison_summary.csv").exists()
        assert list(export_root.rglob("*.parquet"))
        assert list(export_root.rglob("*.png"))
        assert (export_root / "summary").exists()

    def test_run_settings_json_valid(self, mock_st):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path)

        export_root = list(Path(export_dir).iterdir())[0]
        with open(export_root / "run_settings.json") as f:
            settings = json.load(f)
        assert "peptides" in settings

    def test_export_creates_versioned_directory(self, mock_st):
        """Streamlit now uses the CLI's fingerprinted ``{project}_{NNN}_{YYYYMMDD}`` layout."""
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, do_parquet=True)

        export_root = list(Path(export_dir).iterdir())[0]
        # Shared orchestrator: ``MSpecViz_analysis_001_20260420`` (or whatever
        # project name / next number / today is).
        assert export_root.name.startswith("MSpecViz_analysis_")
        parts = export_root.name.rsplit("_", 2)
        assert parts[1].isdigit() and len(parts[1]) == 3
        assert parts[2].isdigit() and len(parts[2]) == 8
