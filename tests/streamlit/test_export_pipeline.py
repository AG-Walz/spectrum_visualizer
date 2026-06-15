"""Integration tests for streamlit_app.views.export.export_actions._run_export (mocked Streamlit)."""
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import polars as pl
import pytest

from core.comparison_metrics import METRIC_SHORT_TO_DISPLAY
from psm_utils import PSMList

from tests.conftest import make_peptide_match
from streamlit_app.views.export.export_actions import _run_export


def _make_psm(peptidoform="PEPTIDEK", source_file="query.mzML"):
    """Build a carrier PSM (peaks on psm.spectrum) — the unit the pipeline streams."""
    mz = np.array([100.0, 200.0, 300.0, 400.0])
    intensity = np.array([1000.0, 2000.0, 3000.0, 500.0])
    psm = make_peptide_match(peptidoform=peptidoform, spectrum_id="index=1",
                             charge=2, mzml_filename=source_file)
    psm.spectrum = {"mz": mz, "intensity": intensity}
    psm.retention_time = 120.0
    psm.precursor_mz = 500.0
    return psm


def _make_annotated_frame():
    return pd.DataFrame({
        "mz": [100.0, 200.0, 300.0, 400.0],
        "intensity": [0.33, 0.67, 1.0, 0.17],
        "annotation": ["b1", "y2", "b3", "unannotated"],
        "is_annotated": [True, True, True, False],
        "ion_type": ["b", "y", "b", "unannotated"],
    })


def _make_comparison_frame():
    """Build a two-row comparison_frame with every metric column the export expects."""
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
        for count_key in (
            "query_peak_count", "reference_peak_count",
            "query_annotation_count", "reference_annotation_count", "shared",
        ):
            row[count_key] = 5
        rows.append(row)
    return pd.DataFrame(rows)


def _make_peptidoform_psms():
    psm = make_peptide_match(
        peptidoform="PEPTIDEK/2", spectrum_id="index=1", charge=2, score=0.95,
        idxml_file="test.idXML", mzml_filename="query.mzML",
    )
    return {"PEPTIDEK": PSMList(psm_list=[psm])}


def _make_psm_lazyframe():
    return pl.DataFrame({
        "id": [0],
        "base_sequence": ["PEPTIDEK"],
        "modification": ["(unmodified)"],
        "display_label": ["PEPTIDEK"],
        "peptidoform": ["PEPTIDEK"],
        "charge": [2],
        "score": [0.95],
        "source_file": ["query.mzML"],
        "spectrum_id": ["index=1"],
        "idxml_file": ["test.idXML"],
    }).lazy()


def _build_session_state(tmp_path):
    """Assemble the peptide_mirror.* session_state _run_export reads end-to-end."""
    query_spectrum = _make_psm()
    query_frame = _make_annotated_frame()
    reference_spectrum = _make_psm(source_file="ref.mzML")
    reference_frame = _make_annotated_frame()

    return {
        "peptide_mirror.annotated_queries": {0: (query_spectrum, query_frame)},
        "peptide_mirror.reference_store": {"PEPTIDEK": [(reference_spectrum, reference_frame)]},
        "peptide_mirror.peptidoform_psms": _make_peptidoform_psms(),
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
    """Drive _run_export outside the live app and check each format produces its artefacts."""

    @pytest.fixture
    def mock_st(self, tmp_path):
        """Patch st + helper modules so _run_export runs end-to-end without a Streamlit runtime."""
        state = _build_session_state(tmp_path)

        with patch("streamlit_app.views.export.export_actions.st") as mock_st, \
             patch("streamlit_app.views.export.export_actions.build_run_settings") as mock_brs, \
             patch("streamlit_app.views.export.export_actions.get_state_key",
                   side_effect=lambda k: f"peptide_mirror.{k}"), \
             patch("streamlit_app.views.export.export_actions.get_param") as mock_gp, \
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
        defaults = dict(export_parquet=False, export_csv=False, export_mirror=False,
                        export_summary=False, export_html=False)
        defaults.update(flags)
        _run_export({}, export_dir, **defaults)
        return export_dir

    @pytest.mark.parametrize("flag, check", [
        # (selection flag, callable(export_root) -> asserts artefact exists)
        ("export_parquet", lambda root: list(root.glob("*.parquet"))),
        ("export_csv",     lambda root: [root / "comparison_summary.csv"]
                                    if (root / "comparison_summary.csv").exists() else []),
        ("export_mirror",  lambda root: list(root.rglob("*.png")) + list(root.rglob("*.svg"))),
        ("export_summary", lambda root: list((root / "summary").rglob("*.*"))
                                    if (root / "summary").exists() else []),
    ], ids=["parquet", "csv", "mirror", "summary"])
    def test_each_format_writes_expected_artefacts(self, mock_st, flag, check):
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, **{flag: True})
        export_root = list(Path(export_dir).iterdir())[0]
        artefacts = check(export_root)
        assert len(artefacts) >= 1, f"{flag} produced no artefacts under {export_root}"
        # run_settings.json is always written regardless of which flags are on
        assert (export_root / "run_settings.json").exists()

    def test_csv_columns_round_trip(self, mock_st):
        """CSV export round-trip: Peptide column survives and row count matches comparison_frame."""
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, export_csv=True)
        export_root = list(Path(export_dir).iterdir())[0]
        frame = pd.read_csv(export_root / "comparison_summary.csv")
        assert len(frame) == 2
        assert "Peptide" in frame.columns

    def test_export_creates_versioned_directory(self, mock_st):
        """Streamlit reuses the CLI's {project}_{NNN}_{YYYYMMDD} versioned layout."""
        _, tmp_path = mock_st
        export_dir = self._invoke_export(tmp_path, export_parquet=True)

        export_root = list(Path(export_dir).iterdir())[0]
        # Expected layout: MSpecViz_analysis_001_20260420 (NNN counter + YYYYMMDD)
        assert export_root.name.startswith("MSpecViz_analysis_")
        parts = export_root.name.rsplit("_", 2)
        assert parts[1].isdigit() and len(parts[1]) == 3
        assert parts[2].isdigit() and len(parts[2]) == 8
