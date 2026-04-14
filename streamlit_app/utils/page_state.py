"""
Analysis-state helpers for MSpecViz view modules.

Handles the on-disk analysis cache (Polars LazyFrames as Parquet, pandas as
Parquet, pickled Python objects, JSON metadata), filesystem listings cached
for the directory browser, and ``run_settings.json`` serialization used by
the Export tab.

Session state primitives (``STATE_PREFIX``, ``get_state_key``, ``get_param``)
live in ``streamlit_app.utils.state`` alongside the rest of the Streamlit
framework plumbing.
"""
import json
import logging
import pickle
import shutil

import pandas as pd
import streamlit as st
import polars as pl

from pathlib import Path
from datetime import datetime

from streamlit_app.utils.state import (
    _persisted_key,
    get_param,
    get_state_key,
    load_params,
)
from mspecviz.common.text_utils import parse_peptide_list

logger = logging.getLogger(__name__)


ANALYSIS_CACHE_DIR = "analysis_cache"

# Keys whose values are Polars LazyFrames (saved as Parquet)
_LAZYFRAME_KEYS = ("psm_lf", "peaks_lf", "seq_lf", "heatmap_lf")

# Keys whose values are pickle-safe Python objects
_PICKLE_KEYS = ("annotated_queries", "reference_store", "variant_groups")


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


def restore_analysis_if_needed():
    """Auto-restore analysis results from disk cache on reload."""
    if not st.session_state.get(get_state_key('data_ready'), False):
        load_analysis_cache()


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
# run_settings.json serialization
# ---------------------------------------------------------------------------
# Single source of truth for which params land in run_settings.json.
# Each entry: (run_settings_key, param_key, default_value).
# ``target-peptides`` needs a list<->newline-joined-string conversion, so it
# is handled separately below.
RUN_SETTINGS_FIELDS: list[tuple[str, str, object]] = [
    ("idxml_files",            "selected-idxmls",        []),
    ("mzml_files",             "selected-mzmls",         []),
    ("reference_substring",    "reference-substring",    ""),
    ("ion_types",              "ion-types",              ["b", "y"]),
    ("mz_tolerance",           "mz-tolerance",           0.02),
    ("add_losses",             "add-losses",             False),
    ("add_first_prefix_ion",   "add-first-prefix-ion",   True),
    ("remove_precursor_peaks", "remove-precursor-peaks", False),
    ("intensity_mode",         "intensity-mode",         "Relative"),
    ("koina_enabled",          "koina-enabled",          False),
]


def build_run_settings() -> dict:
    """Capture current settings for run_settings.json."""
    p = load_params()
    cfg: dict = {
        rs_key: get_param(p, param_key, default)
        for rs_key, param_key, default in RUN_SETTINGS_FIELDS
    }
    cfg["peptides"] = parse_peptide_list(get_param(p, "target-peptides", ""))
    now = datetime.now().isoformat()
    cfg["_created_at"] = now
    cfg["_last_run"] = now
    return cfg


def apply_run_settings(run_cfg: dict):
    """Load a run_settings.json dict back into session state.

    Writes both the widget key and its ``_persisted`` backup so values survive
    page navigation even after Streamlit reaps the owning widget.
    """
    updates: list[tuple[str, object]] = [
        (param_key, run_cfg.get(rs_key, default))
        for rs_key, param_key, default in RUN_SETTINGS_FIELDS
    ]
    updates.append(
        ("target-peptides", "\n".join(run_cfg.get("peptides", [])))
    )

    for param_key, value in updates:
        st.session_state[param_key] = value
        st.session_state[_persisted_key(param_key)] = value
