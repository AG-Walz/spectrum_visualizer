"""
Shared fixtures and markers for MSInsight tests.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

SMALL_DATASET = Path("/mnt/volume/nacho/data/ProCa_runs/small_dataset")
IDXML_PATH = SMALL_DATASET / "ProT001_1_pout_filtered.idXML"
MZML_DIR = SMALL_DATASET
PEPTIDES_PATH = Path("/mnt/volume/nacho/data/peptide_lists/peptides.txt")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: marks tests that require network access (deselect with -m 'not slow')")
    config.addinivalue_line("markers", "integration: marks tests requiring real test data files")


@pytest.fixture
def idxml_path():
    return str(IDXML_PATH)


@pytest.fixture
def mzml_dir():
    return Path("/mnt/volume/nacho/data/ProCa_runs/small_dataset/")


@pytest.fixture
def mzml_paths():
    return sorted(str(p) for p in MZML_DIR.glob("231204_DF_ProCa_ProT001_*.mzML"))


@pytest.fixture
def target_peptides():
    return ["ASDDFDQASY", "RTEFNLNQY"]


@pytest.fixture
def peptides_path():
    return Path("/mnt/volume/nacho/data/peptide_lists/peptides.txt")


@pytest.fixture
def app_settings():
    settings_path = PROJECT_ROOT / "settings.json"
    with open(settings_path) as f:
        return json.load(f)


@pytest.fixture
def spectrum_df():
    """Standard annotated spectrum DataFrame for testing."""
    return pd.DataFrame({
        "mz": [100.0, 200.0, 300.0, 400.0, 500.0],
        "intensity": [1000.0, 2000.0, 3000.0, 1500.0, 500.0],
        "annotation": ["b1", "y2", "b3", "unannotated", ""],
        "is_annotated": [True, True, True, False, False],
    })


@pytest.fixture
def empty_spectrum_df():
    """Empty spectrum DataFrame."""
    return pd.DataFrame({"mz": [], "intensity": [], "annotation": [], "is_annotated": []})
