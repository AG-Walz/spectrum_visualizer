"""
Shared fixtures and markers for MSpecViz tests.
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


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: marks tests that require network access (deselect with -m 'not slow')")
    config.addinivalue_line("markers", "integration: marks tests requiring real test data files")


def pytest_collection_modifyitems(config, items):
    """Auto-skip @pytest.mark.integration tests when test data is unavailable."""
    if IDXML_PATH.exists():
        return
    skip_marker = pytest.mark.skip(reason=f"Test data not found: {SMALL_DATASET}")
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(skip_marker)


@pytest.fixture(scope="session")
def idxml_path():
    return str(IDXML_PATH)


@pytest.fixture(scope="session")
def mzml_paths():
    return sorted(str(p) for p in MZML_DIR.glob("231204_DF_ProCa_ProT001_*.mzML"))


@pytest.fixture(scope="session")
def target_peptides():
    return ["ASDDFDQASY", "RTEFNLNQY"]


@pytest.fixture
def app_settings():
    settings_path = PROJECT_ROOT / "streamlit_app" / "config" / "settings.json"
    with open(settings_path) as f:
        return json.load(f)


def make_peptide_match(peptidoform="PEPTIDEK", base_seq="PEPTIDEK",
                       spectrum_id="index=1", charge=2, score=0.95,
                       idxml_file="test.idXML", mzml_filename="test.mzML"):
    """Create a PeptideMatch for testing. Importable via ``from tests.conftest import ...``."""
    from mspecviz.search import PeptideMatch
    return PeptideMatch(
        peptidoform=peptidoform, base_sequence=base_seq,
        spectrum_id=spectrum_id, charge=charge, score=score,
        idxml_file=idxml_file, mzml_filename=mzml_filename,
        id_merge_index=-1,
    )


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
