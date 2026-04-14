"""
Unit tests for streamlit_app/views/export.py — export pipeline.

Tests _run_export with mocked Streamlit state to verify file creation,
and verifies the export_controls_fragment bug fix (export_dir NameError).
"""
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import polars as pl
import pytest

from mspecviz.spectrum_loader import LoadedSpectrum
from mspecviz.peptide_search import PeptideMatch, VariantGroup
from mspecviz.spectrum_metrics import METRIC_SHORT_TO_DISPLAY


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def _make_annotated_df():
    return pd.DataFrame({
        "mz": [100.0, 200.0, 300.0, 400.0],
        "intensity": [0.33, 0.67, 1.0, 0.17],
        "annotation": ["b1", "y2", "b3", "unannotated"],
        "is_annotated": [True, True, True, False],
        "ion_type": ["b", "y", "b", "unannotated"],
    })


def _make_comparison_df():
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
    vg = VariantGroup(
        base_sequence="PEPTIDEK",
        modification_key=None,
        display_label="PEPTIDEK (unmodified)",
        matches=[match],
    )
    return {"PEPTIDEK|(unmodified)": vg}


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
    query_spec = _make_loaded_spectrum()
    query_df = _make_annotated_df()
    ref_spec = _make_loaded_spectrum(source_file="ref.mzML", is_ref=True)
    ref_df = _make_annotated_df()

    return {
        "peptide_mirror.annotated_queries": {0: (query_spec, query_df)},
        "peptide_mirror.reference_store": {"PEPTIDEK|(unmodified)": [(ref_spec, ref_df)]},
        "peptide_mirror.variant_groups": _make_variant_groups(),
        "peptide_mirror.comparison_df": _make_comparison_df(),
        "peptide_mirror.psm_lf": _make_psm_lazyframe(),
        "peptide_mirror.peaks_lf": pl.DataFrame({"x": [1.0], "y": [1.0]}).lazy(),
        "peptide_mirror.seq_lf": pl.DataFrame({"seq": ["PEP"]}).lazy(),
        "peptide_mirror.heatmap_lf": pl.DataFrame({"x": [1.0], "y": [1.0]}).lazy(),
        "peptide_mirror.data_ready": True,
        "peptide_mirror.cache_version": 1,
        "workspace": tmp_path,
        "settings": {"enable_workspaces": True},
        "export-dir": str(tmp_path / "exports"),
        "intensity-mode": "Relative",
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRunExport:
    """Test _run_export produces expected output files."""

    @pytest.fixture
    def mock_st(self, tmp_path):
        """Patch Streamlit to allow _run_export to run outside of the app."""
        state = _build_session_state(tmp_path)

        with patch("streamlit_app.views.export.st") as mock_st, \
             patch("streamlit_app.views.export.build_run_settings") as mock_brs, \
             patch("streamlit_app.views.export.get_state_key", side_effect=lambda k: f"peptide_mirror.{k}"), \
             patch("streamlit_app.views.export.get_param") as mock_gp:

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

    def test_parquet_export(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=True, do_csv=False,
                    do_mirror=False, do_summary=False, do_html=False)

        export_dirs = list(Path(export_dir).iterdir())
        assert len(export_dirs) == 1
        export_root = export_dirs[0]

        assert (export_root / "run_settings.json").exists()
        parquet_files = list(export_root.glob("*.parquet"))
        assert len(parquet_files) >= 1

    def test_csv_export(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=False, do_csv=True,
                    do_mirror=False, do_summary=False, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        assert (export_root / "comparison_summary.csv").exists()
        df = pd.read_csv(export_root / "comparison_summary.csv")
        assert len(df) == 2
        assert "Peptide" in df.columns

    def test_mirror_plot_export(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=False, do_csv=False,
                    do_mirror=True, do_summary=False, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        png_files = list(export_root.rglob("*.png"))
        svg_files = list(export_root.rglob("*.svg"))
        assert len(png_files) >= 1, "Expected at least one PNG mirror plot"
        assert len(svg_files) >= 1, "Expected at least one SVG mirror plot"

    def test_summary_plot_export(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=False, do_csv=False,
                    do_mirror=False, do_summary=True, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        summary_dir = export_root / "summary"
        assert summary_dir.exists(), "Summary directory should be created"
        summary_files = list(summary_dir.rglob("*.*"))
        assert len(summary_files) >= 2, "Expected dual_axis + ridgeline plots"

    def test_full_export(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=True, do_csv=True,
                    do_mirror=True, do_summary=True, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        assert (export_root / "run_settings.json").exists()
        assert (export_root / "comparison_summary.csv").exists()
        assert list(export_root.rglob("*.parquet"))
        assert list(export_root.rglob("*.png"))
        assert (export_root / "summary").exists()

    def test_run_settings_json_valid(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=False, do_csv=False,
                    do_mirror=False, do_summary=False, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        with open(export_root / "run_settings.json") as f:
            settings = json.load(f)
        assert "peptides" in settings

    def test_export_creates_versioned_directory(self, mock_st):
        from streamlit_app.views.export import _run_export
        _, tmp_path = mock_st
        export_dir = str(tmp_path / "exports")

        _run_export(export_dir, do_parquet=True, do_csv=False,
                    do_mirror=False, do_summary=False, do_html=False)

        export_root = list(Path(export_dir).iterdir())[0]
        assert export_root.name.startswith("export_")


class TestExportControlsBugFix:
    """Verify the NameError fix: export_dir is no longer referenced as undefined."""

    def test_export_controls_no_name_error(self):
        """The function source should reference get_param for export-dir, not bare export_dir."""
        import inspect
        from streamlit_app.views.export import export_controls_fragment
        source = inspect.getsource(export_controls_fragment)
        # The bug was: _run_export(export_dir, ...) where export_dir was undefined.
        # The fix uses get_param(params, 'export-dir', '') or export_fallback.
        assert "get_param" in source or "export_fallback" in source
        # Make sure the bare "export_dir," pattern on _run_export call is gone
        lines = source.split("\n")
        for line in lines:
            if "_run_export(" in line and "export_dir," in line:
                # Should NOT have bare export_dir as first arg
                assert "get_param" in line or "export_fallback" in line, \
                    f"export_dir should not be referenced directly: {line.strip()}"
