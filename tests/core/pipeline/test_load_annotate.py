"""Unit tests for core.pipeline.load_annotate (network-free via monkeypatch)."""
import math

import numpy as np
import pandas as pd

from psm_utils import PSM

from core.comparison_metrics import rt_delta
from core.pipeline import load_annotate


def test_get_koina_prediction_attaches_synthetic_psm(monkeypatch):
    """Step 3: the Koina reference is backed by a PSM carrying the predicted peaks."""
    predicted = pd.DataFrame({"mz": [100.0, 200.0, 300.0], "intensity": [0.5, 1.0, 0.25]})
    # Stub the network call and the annotation so the test stays offline and focused.
    monkeypatch.setattr(load_annotate, "predict_koina_spectrum", lambda *_, **__: predicted)
    monkeypatch.setattr(load_annotate, "prepare_koina_spectrum_for_plot", lambda frame, **_: frame)

    psm, frame = load_annotate.get_koina_prediction("PEPTIDEK", 2, {})

    assert isinstance(psm, PSM)
    assert psm.peptidoform.modified_sequence == "PEPTIDEK"
    assert psm.get_precursor_charge() == 2
    assert psm.run == "Koina:ms2pip_timsTOF2024"
    assert np.array_equal(psm.spectrum["mz"], predicted["mz"].to_numpy())
    assert np.array_equal(psm.spectrum["intensity"], predicted["intensity"].to_numpy())
    # Step 5b: precursor m/z is the theoretical m/z (no measured precursor).
    assert psm.precursor_mz == psm.peptidoform.theoretical_mz


def test_get_koina_prediction_returns_none_on_empty(monkeypatch):
    monkeypatch.setattr(load_annotate, "predict_koina_spectrum", lambda *_, **__: None)
    assert load_annotate.get_koina_prediction("PEPTIDEK", 2, {}) is None


class TestRetentionTimeColumns:
    """RT of each match rides along in the comparison summary (issue #113)."""

    def test_delta_is_query_minus_reference(self):
        assert rt_delta(180.0, 150.0) == 30.0
        assert rt_delta(150.0, 180.0) == -30.0

    def test_delta_is_nan_for_predicted_reference(self):
        # Koina references carry retention_time 0.0 — no real elution time to compare against.
        assert math.isnan(rt_delta(180.0, 0.0))

    def test_delta_is_nan_when_either_rt_is_missing(self):
        assert math.isnan(rt_delta(None, 150.0))
        assert math.isnan(rt_delta(180.0, None))

    def test_rt_columns_are_exported_in_order(self):
        # The RT block sits in the identity columns, so it reaches the CSV/XLSX export.
        for column in load_annotate.COMPARISON_RT_COLUMNS:
            assert column in load_annotate.COMPARISON_SUMMARY_COLUMNS

    def test_build_comparison_frame_keeps_rt_columns(self):
        row = {column: 0 for column in load_annotate.COMPARISON_SUMMARY_COLUMNS}
        row.update({"query_rt": 180.0, "ref_rt": 150.0, "rt_delta": 30.0})
        frame = load_annotate.build_comparison_frame([row])
        assert frame.loc[0, "rt_delta"] == 30.0
        assert list(frame.columns) == list(load_annotate.COMPARISON_SUMMARY_COLUMNS)
