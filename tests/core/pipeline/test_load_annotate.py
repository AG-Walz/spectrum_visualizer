"""Unit tests for core.pipeline.load_annotate (network-free via monkeypatch)."""
import numpy as np
import pandas as pd

from psm_utils import PSM

from core.pipeline import load_annotate


def test_get_koina_prediction_attaches_synthetic_psm(monkeypatch):
    """Step 3: the Koina reference is backed by a PSM carrying the predicted peaks."""
    predicted = pd.DataFrame({"mz": [100.0, 200.0, 300.0], "intensity": [0.5, 1.0, 0.25]})
    # Stub the network call and the annotation so the test stays offline and focused.
    monkeypatch.setattr(load_annotate, "predict_koina_spectrum", lambda **_: predicted)
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
    monkeypatch.setattr(load_annotate, "predict_koina_spectrum", lambda **_: None)
    assert load_annotate.get_koina_prediction("PEPTIDEK", 2, {}) is None
