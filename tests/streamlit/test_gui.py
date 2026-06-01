"""AppTest GUI smoke tests for streamlit_app/ (page launches + Mirror-Plots picker render)."""
from streamlit.testing.v1 import AppTest
import numpy as np
import pandas as pd
import pytest
import json
from pathlib import Path

from psm_utils import PSMList

from core.spectra.loaded_spectrum import LoadedSpectrum
from tests.conftest import make_peptide_match

SETTINGS_PATH = Path(__file__).resolve().parents[2] / "streamlit_app" / "config" / "settings.json"
STATE_PREFIX = "peptide_mirror"


@pytest.fixture
def launch(request):
    """Build an AppTest for the file in request.param and seed session_state.settings."""
    test = AppTest.from_file(request.param)

    # Seed session_state.settings from the real config so page_setup() succeeds.
    with open(SETTINGS_PATH, "r") as f:
        test.session_state.settings = json.load(f)
    return test


# Smoke-test every page boots without raising.
@pytest.mark.slow
@pytest.mark.parametrize(
    "launch",
    (
        "streamlit_app/views/setup_page.py",
        "streamlit_app/views/load_annotate/load_annotate_page.py",
        "streamlit_app/views/explore/explore_page.py",
        "streamlit_app/views/export/export_page.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Each top-level page runs to completion without an uncaught exception."""
    launch.run(timeout=30)
    assert not launch.exception


def _seed_mirror_picker_state(at: AppTest) -> None:
    """Inject one variant + one query + one experimental ref so the picker fragment renders end-to-end.

    Annotation frames carry the columns compute_all_metrics requires so the auto-rendered first-pair plot doesn't blow up."""
    match = make_peptide_match(
        peptidoform="PEPTIDEK", base_seq="PEPTIDEK",
        spectrum_id="scan=42", charge=2, score=0.9, mzml_filename="q.mzML",
    )
    match.retention_time = 12.5
    # Experimental references are PSM-backed too; carry a match so the consumers
    # that read spec.match (labels, picker) work.
    ref_match = make_peptide_match(
        peptidoform="PEPTIDEK", base_seq="PEPTIDEK",
        spectrum_id="scan=10", charge=2, score=0.8, mzml_filename="ref.mzML",
    )
    ref_match.retention_time = 11.0
    spec_kwargs = dict(
        mz=np.array([100.0, 200.0]),
        intensity=np.array([1000.0, 2000.0]),
        precursor_mz=400.0, precursor_charge=2,
        peptidoform="PEPTIDEK", base_sequence="PEPTIDEK",
    )
    annotated_frame = pd.DataFrame({
        "mz": [100.0, 200.0],
        "intensity": [1000.0, 2000.0],
        "annotation": ["b1", "y2"],
        "is_annotated": [True, True],
    })
    query_spec = LoadedSpectrum(spectrum_id="scan=42", source_file="q.mzML",
                                rt=12.5, match=match, **spec_kwargs)
    ref_spec = LoadedSpectrum(spectrum_id="scan=10", source_file="ref.mzML",
                              rt=11.0, is_reference=True, match=ref_match, **spec_kwargs)
    variant_key = "PEPTIDEK"
    variant_groups = {variant_key: PSMList(psm_list=[match])}
    at.session_state[f"{STATE_PREFIX}.variant_groups"] = variant_groups
    at.session_state[f"{STATE_PREFIX}.annotated_queries"] = {0: (query_spec, annotated_frame)}
    at.session_state[f"{STATE_PREFIX}.reference_store"] = {variant_key: [(ref_spec, annotated_frame)]}
    at.session_state[f"{STATE_PREFIX}.data_ready"] = True


@pytest.mark.slow
def test_mirror_picker_renders_when_state_seeded():
    """Mirror-Plots picker fragment runs without exception and shows the 'pick a row' hint."""
    at = AppTest.from_file("streamlit_app/views/explore/explore_page.py")
    with open(SETTINGS_PATH, "r") as f:
        at.session_state.settings = json.load(f)
    _seed_mirror_picker_state(at)
    # Default radio is "Spectrum Explorer" — force "Mirror Plots" so the
    # picker fragment is the one that renders.
    at.session_state[f"{STATE_PREFIX}_explore_mode"] = "Mirror Plots"

    at.run(timeout=30)
    assert not at.exception
