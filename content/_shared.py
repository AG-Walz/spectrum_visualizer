"""
Shared utilities for MSpecViz page modules.

Contains session state helpers, analysis cache management, settings access,
filesystem helpers, and run settings serialization used across all pages.
"""
import io
import json
import logging
import os
import pickle
import shutil
import zipfile

import pandas as pd
import streamlit as st
import polars as pl

from pathlib import Path
from datetime import datetime
from typing import Optional

from src.common.common import load_params
from utils.helpers import format_time, parse_peptide_list

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

def restore_analysis_if_needed():
    """Auto-restore analysis results from disk cache on reload."""
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
    _params = load_params()

    def _p(key, default=None):
        return st.session_state.get(key, _params.get(key, default))

    peptide_text = _p('target-peptides', '')
    peptides = parse_peptide_list(peptide_text)
    return {
        "idxml_files": _p('selected-idxmls', []),
        "mzml_files": _p('selected-mzmls', []),
        "peptides": peptides,
        "reference_substring": _p('reference-substring', ''),
        "ion_types": _p('ion-types', ['b', 'y']),
        "mz_tolerance": _p('mz-tolerance', 0.02),
        "add_losses": _p('add-losses', False),
        "add_first_prefix_ion": _p('add-first-prefix-ion', True),
        "remove_precursor_peaks": _p('remove-precursor-peaks', False),
        "intensity_mode": _p('intensity-mode', 'Relative'),
        "koina_enabled": _p('koina-enabled', False),
        "_created_at": datetime.now().isoformat(),
        "_last_run": datetime.now().isoformat(),
    }


def apply_run_settings(run_cfg: dict):
    """Load a run_settings.json dict back into session state."""
    mapping = {
        'selected-idxmls':        run_cfg.get('idxml_files', []),
        'selected-mzmls':         run_cfg.get('mzml_files', []),
        'target-peptides':        '\n'.join(run_cfg.get('peptides', [])),
        'reference-substring':    run_cfg.get('reference_substring', ''),
        'ion-types':              run_cfg.get('ion_types', ['b', 'y']),
        'mz-tolerance':           run_cfg.get('mz_tolerance', 0.02),
        'add-losses':             run_cfg.get('add_losses', False),
        'add-first-prefix-ion':   run_cfg.get('add_first_prefix_ion', True),
        'remove-precursor-peaks': run_cfg.get('remove_precursor_peaks', False),
        'intensity-mode':         run_cfg.get('intensity_mode', 'Relative'),
        'koina-enabled':          run_cfg.get('koina_enabled', False),
    }
    for key, value in mapping.items():
        st.session_state[key] = value
