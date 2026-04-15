from streamlit.testing.v1 import AppTest
import pytest
import json
from pathlib import Path

SETTINGS_PATH = Path(__file__).resolve().parents[2] / "streamlit_app" / "config" / "settings.json"


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
