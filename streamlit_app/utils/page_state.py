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
from mspecviz.common.helpers import parse_peptide_list
from mspecviz.export import (
    RUN_SETTINGS_FIELDS,
    build_run_settings as _build_run_settings_core,
)

logger = logging.getLogger(__name__)


ANALYSIS_CACHE_DIR = "analysis_cache"

# Keys whose values are Polars LazyFrames (saved as Parquet)
_LAZYFRAME_KEYS = ("psm_frame", "peaks_frame", "sequence_frame", "heatmap_frame")

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
        frame = st.session_state.get(get_state_key(key))
        if frame is not None:
            frame.collect().write_parquet(cache_dir / f"{key}.parquet")

    # 2. comparison_frame (pandas) -> Parquet
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
    if comparison_frame is not None and not comparison_frame.empty:
        comparison_frame.to_parquet(cache_dir / "comparison_frame.parquet")
    elif comparison_frame is not None:
        pd.DataFrame().to_parquet(cache_dir / "comparison_frame.parquet")

    # 3. psm_lookup (dict) -> JSON
    psm_lookup = st.session_state.get(get_state_key("psm_lookup"))
    if psm_lookup is not None:
        with open(cache_dir / "psm_lookup.json", "w") as handle:
            json.dump(psm_lookup, handle)

    # 4. Complex objects -> pickle
    for key in _PICKLE_KEYS:
        cached_value = st.session_state.get(get_state_key(key))
        if cached_value is not None:
            with open(cache_dir / f"{key}.pkl", "wb") as handle:
                pickle.dump(cached_value, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # 5. Metadata marker
    metadata = {
        "saved_at": datetime.now().isoformat(),
        "cache_version": st.session_state.get(get_state_key("cache_version"), 0),
    }
    with open(cache_dir / "meta.json", "w") as handle:
        json.dump(metadata, handle)

    logger.info("Analysis cache saved to %s", cache_dir)


@st.cache_resource(show_spinner="Restoring analysis from cache...")
def _read_analysis_cache(cache_dir_str: str, _cache_modification_time: float):
    """Read all analysis cache files from disk (cached across sessions/users)."""
    cache_dir = Path(cache_dir_str)
    cached_data: dict = {}

    for key in _LAZYFRAME_KEYS:
        parquet_path = cache_dir / f"{key}.parquet"
        if parquet_path.exists():
            cached_data[key] = pl.scan_parquet(parquet_path)

    comparison_path = cache_dir / "comparison_frame.parquet"
    if comparison_path.exists():
        cached_data["comparison_frame"] = pd.read_parquet(comparison_path)

    lookup_path = cache_dir / "psm_lookup.json"
    if lookup_path.exists():
        with open(lookup_path, "r") as handle:
            cached_data["psm_lookup"] = json.load(handle)

    for key in _PICKLE_KEYS:
        pickle_path = cache_dir / f"{key}.pkl"
        if pickle_path.exists():
            with open(pickle_path, "rb") as handle:
                cached_data[key] = pickle.load(handle)

    with open(cache_dir / "meta.json", "r") as handle:
        cached_data["_meta"] = json.load(handle)

    return cached_data


def load_analysis_cache() -> bool:
    """Restore analysis results from the disk cache into session state."""
    cache_dir = _analysis_cache_path()
    meta_path = cache_dir / "meta.json"

    if not meta_path.exists():
        return False

    try:
        modification_time = meta_path.stat().st_mtime
        cached_data = _read_analysis_cache(str(cache_dir), modification_time)

        for key in (*_LAZYFRAME_KEYS, "comparison_frame", "psm_lookup", *_PICKLE_KEYS):
            if key in cached_data:
                st.session_state[get_state_key(key)] = cached_data[key]

        metadata = cached_data.get("_meta", {})
        st.session_state[get_state_key("cache_version")] = metadata.get("cache_version", 1)
        st.session_state[get_state_key("data_ready")] = True

        logger.info("Analysis cache restored from %s", cache_dir)
        return True

    except Exception as error:
        logger.warning("Failed to restore analysis cache: %s", error)
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
        files: list = []
        for extension in extensions:
            iterator = (Path(path).rglob(f"*{extension}") if recursive
                        else Path(path).glob(f"*{extension}"))
            for found_path in iterator:
                files.append(found_path)
                if len(files) >= _max_files:
                    break
            if len(files) >= _max_files:
                break
        return sorted(str(found_path) for found_path in files)
    except (PermissionError, OSError):
        return []


# ---------------------------------------------------------------------------
# run_settings.json serialization
# ---------------------------------------------------------------------------
# The canonical field list lives in ``mspecviz.export.RUN_SETTINGS_FIELDS``.
# Streamlit stores some of those values under different session keys (e.g.
# ``idxml_files`` → ``selected_idxmls``) — the map below covers the differences.
_STREAMLIT_PARAM_KEY_OVERRIDES: dict[str, str] = {
    "idxml_files": "selected_idxmls",
    "mzml_files":  "selected_mzmls",
}


def _streamlit_param_key(run_settings_key: str) -> str:
    """Return the Streamlit session-state key for a run_settings field."""
    return _STREAMLIT_PARAM_KEY_OVERRIDES.get(run_settings_key, run_settings_key)


def build_run_settings() -> dict:
    """Capture current Streamlit settings for run_settings.json."""
    params = load_params()
    return _build_run_settings_core(
        getter=lambda key, default: get_param(params, _streamlit_param_key(key), default),
        peptides=parse_peptide_list(get_param(params, "target_peptides", "")),
    )


def apply_run_settings(run_settings: dict):
    """Load a run_settings.json dict back into session state.

    Writes both the widget key and its ``_persisted`` backup so values survive
    page navigation even after Streamlit reaps the owning widget.
    """
    updates: list[tuple[str, object]] = [
        (_streamlit_param_key(rs_key), run_settings.get(rs_key, default))
        for rs_key, default in RUN_SETTINGS_FIELDS
    ]
    updates.append(
        ("target_peptides", "\n".join(run_settings.get("peptides", [])))
    )

    for param_key, value in updates:
        st.session_state[param_key] = value
        st.session_state[_persisted_key(param_key)] = value
