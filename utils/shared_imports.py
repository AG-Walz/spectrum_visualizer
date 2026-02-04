"""
Centralized imports for heavy packages.

Import from this module instead of importing directly to avoid redundant loading.
Heavy packages (PyOpenMS, pandas, numpy, etc.) are imported once and cached.
"""
import streamlit as st


@st.cache_resource
def _load_heavy_imports():
    """Load heavy packages once and cache them."""
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import logging

    # PyOpenMS imports
    import pyopenms as oms

    # psm_utils imports
    from psm_utils.io import read_file
    from psm_utils import Peptidoform, PSMList

    return {
        'pd': pd,
        'np': np,
        'Path': Path,
        'logging': logging,
        'oms': oms,
        'read_file': read_file,
        'Peptidoform': Peptidoform,
        'PSMList': PSMList,
    }


# Module-level access for convenience
_imports = None

def _get_imports():
    global _imports
    if _imports is None:
        _imports = _load_heavy_imports()
    return _imports


# Lazy loading accessors
def get_pandas():
    return _get_imports()['pd']

def get_numpy():
    return _get_imports()['np']

def get_path():
    return _get_imports()['Path']

def get_logging():
    return _get_imports()['logging']

def get_pyopenms():
    return _get_imports()['oms']

def get_psm_utils():
    imports = _get_imports()
    return imports['read_file'], imports['Peptidoform'], imports['PSMList']
