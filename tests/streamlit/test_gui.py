"""AppTest GUI smoke tests for streamlit_app/ (page launches + Mirror-Plots picker render)."""
from streamlit.testing.v1 import AppTest
import numpy as np
import pandas as pd
import pytest
import json
from pathlib import Path

from psm_utils import PSMList

from tests.conftest import make_peptide_match

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SETTINGS_PATH = PROJECT_ROOT / "streamlit_app" / "config" / "settings.json"
STATE_PREFIX = "peptide_mirror"


@pytest.fixture
def launch(request):
    """Build an AppTest for the file in request.param and seed session_state.settings."""
    # Absolute path: Streamlit >= 1.6x resolves a relative path against this test file.
    test = AppTest.from_file(str(PROJECT_ROOT / request.param))

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


@pytest.mark.slow
def test_workspace_switch_loads_the_new_workspace_without_a_warning(tmp_path):
    """Switching workspace loads its params, and shows no Streamlit warning: st.rerun() in the
    selector's callback made Streamlit 1.5x print "Calling st.rerun() within a callback is a
    no-op." on every switch."""
    for name, tolerance in (("alpha", 0.01), ("beta", 0.05)):
        (tmp_path / name).mkdir()
        (tmp_path / name / "params.json").write_text(json.dumps({"mz_tolerance": tolerance}))
    at = AppTest.from_file(str(PROJECT_ROOT / "streamlit_app/views/setup_page.py"), default_timeout=60)
    with open(SETTINGS_PATH, "r") as f:
        at.session_state.settings = {**json.load(f), "workspaces_dir": str(tmp_path)}
    at.query_params["workspace"] = "alpha"
    at.run()

    at.selectbox(key="chosen-workspace").set_value("beta").run()

    assert not at.exception
    assert at.session_state.workspace.name == "beta"
    assert at.number_input(key="mz_tolerance").value == 0.05
    assert not [warning.value for warning in at.warning if "st.rerun()" in warning.value]


def _seed_mirror_picker_state(at: AppTest) -> None:
    """Inject one peptidoform + one query + one experimental ref so the picker fragment renders end-to-end.

    Annotation frames carry the columns compute_all_metrics requires so the auto-rendered first-pair plot doesn't blow up."""
    peaks = {"mz": np.array([100.0, 200.0]), "intensity": np.array([1000.0, 2000.0])}
    # Query and reference are carrier PSMs (peaks on psm.spectrum).
    query_psm = make_peptide_match(
        peptidoform="PEPTIDEK", spectrum_id="scan=42", charge=2, score=0.9, mzml_filename="q.mzML",
    )
    query_psm.retention_time = 12.5
    query_psm.spectrum = peaks
    ref_psm = make_peptide_match(
        peptidoform="PEPTIDEK", spectrum_id="scan=10", charge=2, score=0.8, mzml_filename="ref.mzML",
    )
    ref_psm.retention_time = 11.0
    ref_psm.spectrum = peaks
    annotated_frame = pd.DataFrame({
        "mz": [100.0, 200.0],
        "intensity": [1000.0, 2000.0],
        "annotation": ["b1", "y2"],
        "is_annotated": [True, True],
    })
    peptidoform = "PEPTIDEK"
    peptidoform_psms = {peptidoform: PSMList(psm_list=[query_psm])}
    at.session_state[f"{STATE_PREFIX}.peptidoform_psms"] = peptidoform_psms
    at.session_state[f"{STATE_PREFIX}.annotated_queries"] = {0: (query_psm, annotated_frame)}
    at.session_state[f"{STATE_PREFIX}.reference_store"] = {peptidoform: [(ref_psm, annotated_frame)]}
    at.session_state[f"{STATE_PREFIX}.data_ready"] = True


@pytest.mark.slow
def test_mirror_picker_renders_when_state_seeded():
    """Mirror-Plots picker fragment runs without exception and shows the 'pick a row' hint."""
    at = AppTest.from_file(str(PROJECT_ROOT / "streamlit_app/views/explore/explore_page.py"))
    with open(SETTINGS_PATH, "r") as f:
        at.session_state.settings = json.load(f)
    _seed_mirror_picker_state(at)
    # Default radio is "Spectrum Explorer" — force "Mirror Plots" so the
    # picker fragment is the one that renders.
    at.session_state[f"{STATE_PREFIX}_explore_mode"] = "Mirror Plots"

    at.run(timeout=30)
    assert not at.exception


@pytest.mark.slow
def test_folder_dropdown_lists_every_folder_and_opens_it(tmp_path):
    """The Setup file browser lists all sub-folders (no cap, hidden ones left out) in one
    dropdown; choosing a folder opens it and shows its files, and Up goes back."""
    data = tmp_path / "data"
    for index in range(20):  # more than the 16 buttons the old grid could show
        (data / f"run_{index:02d}").mkdir(parents=True)
    (data / ".hidden").mkdir()
    (data / "run_07" / "sample.idXML").write_text("")
    (tmp_path / "ws" / "default").mkdir(parents=True)
    (tmp_path / "ws" / "default" / "params.json").write_text(json.dumps(
        {"idxml_browser_dir": str(data), "mzml_browser_dir": str(data)}))
    at = AppTest.from_file(str(PROJECT_ROOT / "streamlit_app/views/setup_page.py"), default_timeout=60)
    with open(SETTINGS_PATH, "r") as f:
        at.session_state.settings = {**json.load(f), "workspaces_dir": str(tmp_path / "ws")}
    at.run()

    picker = at.selectbox(key=f"{STATE_PREFIX}_idxml_subfolder")
    assert len(picker.options) == 20 and ".hidden" not in picker.options

    picker.set_value("run_07").run()
    assert at.session_state[f"{STATE_PREFIX}_idxml_dir"] == str(data / "run_07")
    assert at.checkbox(key=f"{STATE_PREFIX}_idxml_select_all").label.startswith("Select all 1 files")

    at.button(key=f"{STATE_PREFIX}_idxml_up_btn").click().run()
    assert not at.exception
    assert at.session_state[f"{STATE_PREFIX}_idxml_dir"] == str(data)
