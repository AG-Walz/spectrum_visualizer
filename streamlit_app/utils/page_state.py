"""Analysis-cache I/O, directory-browser listings, and run_settings.json 
(de)serialization — used by Load & Annotate (write), Explore (read), and Export (read + settings).
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
from core.helpers import parse_peptide_list
from core.export import (
    RUN_SETTINGS_FIELDS,
    build_run_settings as _build_run_settings_core,
)

logger = logging.getLogger(__name__)


ANALYSIS_CACHE_DIR = "analysis_cache"

# Keys whose values are Polars LazyFrames (saved as Parquet)
_LAZYFRAME_KEYS = ("psm_frame", "peaks_frame", "sequence_frame", "heatmap_frame")

# Keys whose values are pickle-safe Python objects
_PICKLE_KEYS = ("annotated_queries", "reference_store", "grouped_psms")

# On-disk cache layout version. Bump whenever the pickled object shape changes
# (e.g. PSMs gaining their peak arrays in the PSMList-native migration) so caches
# written by an older code version are cleared and regenerated from inputs rather
# than hydrating semantically-stale objects. Distinct from the runtime
# "cache_version" re-render counter in session state.
CACHE_SCHEMA_VERSION = 1


def _analysis_cache_path() -> Path:
    """Return the analysis-cache directory inside the active workspace.
    """
    return Path(st.session_state.workspace) / ANALYSIS_CACHE_DIR


def save_analysis_cache() -> None:
    """Write every Load-&-Annotate result to analysis_cache/ so the Explore and Export views survive a reload.
    """
    cache_dir = _analysis_cache_path()
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Polars LazyFrames -> Parquet (collect each frame to materialise it for write)
    for key in _LAZYFRAME_KEYS:
        frame = st.session_state.get(get_state_key(key))
        if frame is not None:
            frame.collect().write_parquet(cache_dir / f"{key}.parquet")

    # comparison_frame (pandas) -> Parquet; write an empty frame as a placeholder for "no metrics"
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
    if comparison_frame is not None and not comparison_frame.empty:
        comparison_frame.to_parquet(cache_dir / "comparison_frame.parquet")
    elif comparison_frame is not None:
        pd.DataFrame().to_parquet(cache_dir / "comparison_frame.parquet")

    # psm_lookup dict (spectrum-id -> PSM metadata) -> JSON
    psm_lookup = st.session_state.get(get_state_key("psm_lookup"))
    if psm_lookup is not None:
        with open(cache_dir / "psm_lookup.json", "w") as handle:
            json.dump(psm_lookup, handle)

    # Complex Python objects (annotated_queries, reference_store, grouped_psms) -> pickle
    for key in _PICKLE_KEYS:
        cached_value = st.session_state.get(get_state_key(key))
        if cached_value is not None:
            with open(cache_dir / f"{key}.pkl", "wb") as handle:
                pickle.dump(cached_value, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Metadata marker — its mtime drives @st.cache_resource invalidation in _read_analysis_cache
    metadata = {
        "saved_at": datetime.now().isoformat(),
        "cache_version": st.session_state.get(get_state_key("cache_version"), 0),
        "schema_version": CACHE_SCHEMA_VERSION,
    }
    with open(cache_dir / "meta.json", "w") as handle:
        json.dump(metadata, handle)

    logger.info("Analysis cache saved to %s", cache_dir)


@st.cache_resource(show_spinner="Restoring analysis from cache...")
def _read_analysis_cache(cache_dir_str: str, _cache_modification_time: float):
    """Read every analysis-cache file from disk; _cache_modification_time 
    busts @st.cache_resource when files change.
    """
    cache_dir = Path(cache_dir_str)
    cached_data: dict = {}

    # LazyFrames: scan_parquet each visualization frame
    for key in _LAZYFRAME_KEYS:
        parquet_path = cache_dir / f"{key}.parquet"
        if parquet_path.exists():
            cached_data[key] = pl.scan_parquet(parquet_path)

    # Comparison frame: pandas Parquet (metrics table for the Export tab)
    comparison_path = cache_dir / "comparison_frame.parquet"
    if comparison_path.exists():
        cached_data["comparison_frame"] = pd.read_parquet(comparison_path)

    # psm_lookup: JSON dict mapping spectrum ID to PSM metadata
    lookup_path = cache_dir / "psm_lookup.json"
    if lookup_path.exists():
        with open(lookup_path, "r") as handle:
            cached_data["psm_lookup"] = json.load(handle)

    # Pickled objects: annotated_queries, reference_store, grouped_psms
    for key in _PICKLE_KEYS:
        pickle_path = cache_dir / f"{key}.pkl"
        if pickle_path.exists():
            with open(pickle_path, "rb") as handle:
                cached_data[key] = pickle.load(handle)

    # Metadata marker: cache_version + saved_at timestamp
    with open(cache_dir / "meta.json", "r") as handle:
        cached_data["_meta"] = json.load(handle)

    return cached_data


def load_analysis_cache() -> bool:
    """Hydrate the on-disk analysis cache back into session_state — 
    returns False if nothing to restore.
    """
    # Discover: locate workspace cache dir, bail if no metadata marker
    cache_dir = _analysis_cache_path()
    meta_path = cache_dir / "meta.json"

    if not meta_path.exists():
        return False

    # Schema gate: a cache written by an older data layout (e.g. before PSMs
    # carried their peak arrays) is cleared and regenerated from inputs rather
    # than hydrating semantically-stale objects. Checked before any unpickling.
    try:
        with open(meta_path) as handle:
            on_disk_schema = json.load(handle).get("schema_version")
    except (OSError, ValueError):
        on_disk_schema = None
    if on_disk_schema != CACHE_SCHEMA_VERSION:
        logger.info(
            "Analysis cache schema_version=%s != current %s; clearing to regenerate.",
            on_disk_schema, CACHE_SCHEMA_VERSION,
        )
        clear_analysis_cache()
        return False

    try:
        # Read & restore: hydrate every cached frame/object back into session_state
        modification_time = meta_path.stat().st_mtime
        cached_data = _read_analysis_cache(str(cache_dir), modification_time)

        for key in (*_LAZYFRAME_KEYS, "comparison_frame", "psm_lookup", *_PICKLE_KEYS):
            if key in cached_data:
                st.session_state[get_state_key(key)] = cached_data[key]

        # Flags: mark data ready and propagate cache version
        metadata = cached_data.get("_meta", {})
        st.session_state[get_state_key("cache_version")] = metadata.get("cache_version", 1)
        st.session_state[get_state_key("data_ready")] = True

        logger.info("Analysis cache restored from %s", cache_dir)
        return True

    except (ModuleNotFoundError, AttributeError, ImportError) as error:
        # Pickle was written against an older module layout (e.g. the loader_models ->
        # loaded_spectrum rename). Cache is permanently unloadable — wipe it so the user
        # doesn't see this warning on every rerun.
        logger.warning(
            "Analysis cache from a previous code version is unreadable (%s); "
            "clearing it. Re-run Load & Annotate to regenerate.", error,
        )
        clear_analysis_cache()
        return False
    except Exception as error:
        logger.warning("Failed to restore analysis cache: %s", error)
        return False


def clear_analysis_cache() -> None:
    """Delete the workspace cache directory and invalidate the 
    in-memory @st.cache_resource entry.
    """
    cache_dir = _analysis_cache_path()
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
        _read_analysis_cache.clear()
        logger.info("Analysis cache cleared: %s", cache_dir)


def restore_analysis_if_needed():
    """Hydrate the analysis cache when a view loads on a fresh rerun without data_ready set.
    """
    if not st.session_state.get(get_state_key('data_ready'), False):
        load_analysis_cache()


@st.cache_data(ttl=10)
def list_subdirectories(path: str) -> list:
    """Directory-browser helper: list immediate subdirectories of path (cached 10s).
    """
    try:
        return sorted(str(d) for d in Path(path).iterdir() if d.is_dir())
    except (PermissionError, OSError):
        return []


@st.cache_data(ttl=10)
def list_files(path: str, extensions: tuple, recursive: bool = False, _max_files: int = 5000) -> list:
    """Directory-browser helper: list files in path matching extensions, capped at _max_files (cached 10s).
    """
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
# Canonical field list lives in core.export.RUN_SETTINGS_FIELDS. A few of those
# fields live under different session-state keys in Streamlit (e.g. idxml_files ->
# selected_idxmls); the map below names the differences.
_STREAMLIT_PARAM_KEY_OVERRIDES: dict[str, str] = {
    "idxml_files": "selected_idxmls",
    "mzml_files":  "selected_mzmls",
}


def _streamlit_param_key(run_settings_key: str) -> str:
    """Translate a run_settings.json field name to its Streamlit session-state key.
    """
    return _STREAMLIT_PARAM_KEY_OVERRIDES.get(run_settings_key, run_settings_key)


def build_run_settings() -> dict:
    """Snapshot the active Setup widgets into a run_settings.json-shaped dict for the Export tab.
    """
    params = load_params()
    return _build_run_settings_core(
        getter=lambda key, default: get_param(params, _streamlit_param_key(key), default),
        peptides=parse_peptide_list(get_param(params, "target_peptides", "")),
    )


def apply_run_settings(run_settings: dict):
    """Restore a previously exported run_settings.json into session_state — 
    writes the widget key and its _persisted backup so values survive page navigation.
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
