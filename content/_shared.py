"""
Shared utilities for MSpecViz page modules.

Contains session state helpers, analysis cache management, settings access,
filesystem helpers, and run settings serialization used across all pages.
"""
import io
import json
import logging
import os
import shutil
import zipfile

import streamlit as st
import polars as pl

from pathlib import Path
from datetime import datetime
from typing import Optional

from utils.helpers import format_time, get_nested as get_settings_value, parse_peptide_list

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Session state key management
# ---------------------------------------------------------------------------

STATE_PREFIX = "peptide_mirror"


def get_state_key(key: str) -> str:
    """Get namespaced session state key."""
    return f"{STATE_PREFIX}.{key}"


def _get_settings() -> dict:
    """Get settings from session state (loaded by page_setup)."""
    return st.session_state.get("settings", {})


# ---------------------------------------------------------------------------
# Analysis cache constants
# ---------------------------------------------------------------------------

ANALYSIS_CACHE_DIR = "analysis_cache"

# Keys whose values are Polars LazyFrames (saved as Parquet)
_LAZYFRAME_KEYS = ("psm_lf", "peaks_lf", "seq_lf", "heatmap_lf")

# Keys whose values are pickle-safe Python objects
_PICKLE_KEYS = ("annotated_queries", "reference_store", "variant_groups")


# ---------------------------------------------------------------------------
# Analysis cache management
# ---------------------------------------------------------------------------

def _analysis_cache_path() -> Path:
    """Return the analysis cache directory inside the current workspace."""
    return Path(st.session_state.workspace) / ANALYSIS_CACHE_DIR


def save_analysis_cache() -> None:
    """Persist all analysis results to disk so they survive page reloads."""
    import pickle
    import pandas as pd

    cache_dir = _analysis_cache_path()
    cache_dir.mkdir(parents=True, exist_ok=True)

    # 1. Polars LazyFrames -> Parquet
    for key in _LAZYFRAME_KEYS:
        lf = st.session_state.get(get_state_key(key))
        if lf is not None:
            lf.collect().write_parquet(cache_dir / f"{key}.parquet")

    # 2. comparison_df (pandas) -> Parquet
    comp_df = st.session_state.get(get_state_key("comparison_df"))
    if comp_df is not None and not comp_df.empty:
        comp_df.to_parquet(cache_dir / "comparison_df.parquet")
    elif comp_df is not None:
        pd.DataFrame().to_parquet(cache_dir / "comparison_df.parquet")

    # 3. psm_lookup (dict) -> JSON
    psm_lookup = st.session_state.get(get_state_key("psm_lookup"))
    if psm_lookup is not None:
        with open(cache_dir / "psm_lookup.json", "w") as f:
            json.dump(psm_lookup, f)

    # 4. Complex objects -> pickle
    for key in _PICKLE_KEYS:
        obj = st.session_state.get(get_state_key(key))
        if obj is not None:
            with open(cache_dir / f"{key}.pkl", "wb") as f:
                pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)

    # 5. Metadata marker
    meta = {
        "saved_at": datetime.now().isoformat(),
        "cache_version": st.session_state.get(get_state_key("cache_version"), 0),
    }
    with open(cache_dir / "meta.json", "w") as f:
        json.dump(meta, f)

    logger.info("Analysis cache saved to %s", cache_dir)


@st.cache_resource(show_spinner="Restoring analysis from cache...")
def _read_analysis_cache(cache_dir_str: str, _cache_mtime: float):
    """Read all analysis cache files from disk (cached across sessions/users)."""
    import pickle
    import pandas as pd

    cache_dir = Path(cache_dir_str)
    data = {}

    for key in _LAZYFRAME_KEYS:
        p = cache_dir / f"{key}.parquet"
        if p.exists():
            data[key] = pl.scan_parquet(p)

    comp_path = cache_dir / "comparison_df.parquet"
    if comp_path.exists():
        data["comparison_df"] = pd.read_parquet(comp_path)

    lookup_path = cache_dir / "psm_lookup.json"
    if lookup_path.exists():
        with open(lookup_path, "r") as f:
            data["psm_lookup"] = json.load(f)

    for key in _PICKLE_KEYS:
        pkl_path = cache_dir / f"{key}.pkl"
        if pkl_path.exists():
            with open(pkl_path, "rb") as f:
                data[key] = pickle.load(f)

    with open(cache_dir / "meta.json", "r") as f:
        data["_meta"] = json.load(f)

    return data


def load_analysis_cache() -> bool:
    """Restore analysis results from the disk cache into session state."""
    cache_dir = _analysis_cache_path()
    meta_path = cache_dir / "meta.json"

    if not meta_path.exists():
        return False

    try:
        mtime = meta_path.stat().st_mtime
        data = _read_analysis_cache(str(cache_dir), mtime)

        for key in list(_LAZYFRAME_KEYS) + ["comparison_df", "psm_lookup"] + list(_PICKLE_KEYS):
            if key in data:
                st.session_state[get_state_key(key)] = data[key]

        meta = data.get("_meta", {})
        st.session_state[get_state_key("cache_version")] = meta.get("cache_version", 1)
        st.session_state[get_state_key("data_ready")] = True

        logger.info("Analysis cache restored from %s", cache_dir)
        return True

    except Exception as e:
        logger.warning("Failed to restore analysis cache: %s", e)
        return False


def clear_analysis_cache() -> None:
    """Remove the analysis cache directory and invalidate the in-memory cache."""
    cache_dir = _analysis_cache_path()
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
        _read_analysis_cache.clear()
        logger.info("Analysis cache cleared: %s", cache_dir)


# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------

def init_session_state(params: dict = None):
    """Initialize session state with defaults from settings.json."""
    if params:
        for key, value in params.get("peptide_mirror", {}).items():
            if key.startswith("peptide_mirror.") and key not in st.session_state:
                st.session_state[key] = value

    _settings = _get_settings()
    ion_types_default = get_settings_value(_settings, "ion_types", "default", default=["b", "y"])
    mz_tolerance_default = get_settings_value(_settings, "mz_tolerance", "default", default=0.02)
    workspace = st.session_state.get('workspace', None)
    export_dir_default = (
        str(Path(workspace) / "exports")
        if workspace
        else get_settings_value(_settings, "export", "default_directory", default="./exports")
    )
    add_losses_default = get_settings_value(_settings, "annotation", "add_losses", default=False)
    add_first_prefix_default = get_settings_value(_settings, "annotation", "add_first_prefix_ion", default=True)
    remove_precursor_default = get_settings_value(_settings, "annotation", "remove_precursor_peaks", default=False)
    workers_default = get_settings_value(_settings, "performance", "workers", default=4)

    defaults = {
        'selected_idxmls': [],
        'selected_mzmls': [],
        'target_peptides': '',
        'reference_substring': '',
        'ion_types': ion_types_default,
        'mz_tolerance': mz_tolerance_default,
        'add_losses': add_losses_default,
        'add_first_prefix_ion': add_first_prefix_default,
        'remove_precursor_peaks': remove_precursor_default,
        'export_dir': export_dir_default,
        'workers': workers_default,
        'intensity_mode': 'Relative',
        'variant_groups': None,
        'data_ready': False,
        'koina_enabled': False,
    }

    for key, default in defaults.items():
        state_key = get_state_key(key)
        if state_key not in st.session_state:
            st.session_state[state_key] = default

    # Auto-restore analysis results from disk cache on reload
    if not st.session_state.get(get_state_key('data_ready'), False):
        load_analysis_cache()


# ---------------------------------------------------------------------------
# Cached filesystem helpers
# ---------------------------------------------------------------------------

@st.cache_data(ttl=10)
def list_subdirectories(path: str) -> list:
    try:
        return sorted(str(d) for d in Path(path).iterdir() if d.is_dir())
    except (PermissionError, OSError):
        return []


@st.cache_data(ttl=10)
def list_files(path: str, extensions: tuple, recursive: bool = False, _max_files: int = 5000) -> list:
    try:
        files = []
        for ext in extensions:
            gen = Path(path).rglob(f"*{ext}") if recursive else Path(path).glob(f"*{ext}")
            for f in gen:
                files.append(f)
                if len(files) >= _max_files:
                    break
            if len(files) >= _max_files:
                break
        return sorted(str(f) for f in files)
    except (PermissionError, OSError):
        return []


# ---------------------------------------------------------------------------
# Run settings serialization
# ---------------------------------------------------------------------------

def build_run_settings() -> dict:
    """Capture current settings for run_settings.json."""
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    peptides = parse_peptide_list(peptide_text)
    return {
        "idxml_files": st.session_state.get(get_state_key('selected_idxmls'), []),
        "mzml_files": st.session_state.get(get_state_key('selected_mzmls'), []),
        "peptides": peptides,
        "reference_substring": st.session_state.get(get_state_key('reference_substring'), ''),
        "ion_types": st.session_state.get(get_state_key('ion_types'), ['b', 'y']),
        "mz_tolerance": st.session_state.get(get_state_key('mz_tolerance'), 0.02),
        "add_losses": st.session_state.get(get_state_key('add_losses'), False),
        "add_first_prefix_ion": st.session_state.get(get_state_key('add_first_prefix_ion'), True),
        "remove_precursor_peaks": st.session_state.get(get_state_key('remove_precursor_peaks'), False),
        "intensity_mode": st.session_state.get(get_state_key('intensity_mode'), 'Relative'),
        "koina_enabled": st.session_state.get(get_state_key('koina_enabled'), False),
        "_created_at": datetime.now().isoformat(),
        "_last_run": datetime.now().isoformat(),
    }


def apply_run_settings(run_cfg: dict):
    """Load a run_settings.json dict back into session state."""
    mapping = {
        'selected_idxmls':        run_cfg.get('idxml_files', []),
        'selected_mzmls':         run_cfg.get('mzml_files', []),
        'target_peptides':        '\n'.join(run_cfg.get('peptides', [])),
        'reference_substring':    run_cfg.get('reference_substring', ''),
        'ion_types':              run_cfg.get('ion_types', ['b', 'y']),
        'mz_tolerance':           run_cfg.get('mz_tolerance', 0.02),
        'add_losses':             run_cfg.get('add_losses', False),
        'add_first_prefix_ion':   run_cfg.get('add_first_prefix_ion', True),
        'remove_precursor_peaks': run_cfg.get('remove_precursor_peaks', False),
        'intensity_mode':         run_cfg.get('intensity_mode', 'Relative'),
        'koina_enabled':          run_cfg.get('koina_enabled', False),
    }
    for key, value in mapping.items():
        st.session_state[get_state_key(key)] = value
