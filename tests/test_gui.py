from streamlit.testing.v1 import AppTest
import pytest
import json


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    return test


# Test launching of all pages
@pytest.mark.parametrize(
    "launch",
    (
        "content/peptide_mirror.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)
    assert not launch.exception
