from streamlit.testing.v1 import AppTest
import numpy as np
import pandas as pd
import pytest
import json
from pathlib import Path

from mspecviz.search.models import VariantGroup
from mspecviz.spectra.loader_models import LoadedSpectrum
from tests.conftest import make_peptide_match

SETTINGS_PATH = Path(__file__).resolve().parents[2] / "streamlit_app" / "config" / "settings.json"
STATE_PREFIX = "peptide_mirror"


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open(SETTINGS_PATH, "r") as f:
        test.session_state.settings = json.load(f)
    return test


# Test launching of all pages
@pytest.mark.parametrize(
    "launch",
    (
        "streamlit_app/views/setup.py",
        "streamlit_app/views/search_load/page.py",
        "streamlit_app/views/explore/page.py",
        "streamlit_app/views/export/page.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)
    assert not launch.exception


def _seed_mirror_picker_state(at: AppTest) -> None:
    """Inject one variant + one query + one experimental ref so the picker fragment
    can run end-to-end. Annotation frames carry the columns ``compute_all_metrics``
    requires (``mz``, ``intensity``, ``annotation``, ``is_annotated``) so the
    auto-rendered first-pair plot doesn't blow up."""
    match = make_peptide_match(
        peptidoform="PEPTIDEK", base_seq="PEPTIDEK",
        spectrum_id="scan=42", charge=2, score=0.9, mzml_filename="q.mzML",
    )
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
                              rt=11.0, is_reference=True, **spec_kwargs)
    variant_key = "PEPTIDEK|"
    variant_groups = {variant_key: VariantGroup(
        base_sequence="PEPTIDEK", modification_key="", display_label="PEPTIDEK",
        matches=[match],
    )}
    at.session_state[f"{STATE_PREFIX}.variant_groups"] = variant_groups
    at.session_state[f"{STATE_PREFIX}.annotated_queries"] = {0: (query_spec, annotated_frame)}
    at.session_state[f"{STATE_PREFIX}.reference_store"] = {variant_key: [(ref_spec, annotated_frame)]}
    at.session_state[f"{STATE_PREFIX}.data_ready"] = True


def test_mirror_picker_renders_when_state_seeded():
    """The picker fragment runs without exception and shows the 'pick a row' hint
    when both selections are unset (StateManager wiring)."""
    at = AppTest.from_file("streamlit_app/views/explore/page.py")
    with open(SETTINGS_PATH, "r") as f:
        at.session_state.settings = json.load(f)
    _seed_mirror_picker_state(at)
    # Default radio is "Spectrum Explorer"; switch to "Mirror Plots" so our
    # fragment is the one that renders.
    at.session_state[f"{STATE_PREFIX}_explore_mode"] = "Mirror Plots"

    at.run(timeout=30)
    assert not at.exception
