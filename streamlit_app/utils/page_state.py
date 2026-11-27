"""Analysis-cache I/O, directory-browser listings, and run_settings.json 
(de)serialization — used by Load & Annotate (write), Explore (read), and Export (read + settings).
"""
import json
import logging
import os
import pickle
import shutil
import uuid

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
from core.defaults import DEFAULTS
from core.helpers import parse_peptide_list
from core.export import (
    RUN_SETTINGS_FIELDS,
    build_run_settings as _build_run_settings_core,
)

logger = logging.getLogger(__name__)


ANALYSIS_CACHE_DIR = "analysis_cache"

# openms_insight component caches (Explorer + mirror pickers), per workspace. Not under
# exports/, so the export ZIP and the export history hold only real exports.
COMPONENT_CACHE_DIR = "component_cache"

# Friendly labels for the summary plot-type keys stored in run_settings —
# shared by the Export controls and the Explore plot picker.
SUMMARY_PLOT_LABELS = {
    "ridgeline": "Ridgeline density",
    "multi_boxplot": "Boxplot",
}

# Friendly labels for the peptidoform sort modes (keys are core.plotting.SUMMARY_SORT_MODES).
# Insertion order is the option order shown in the picker; the first entry (count_asc)
# is the default and reproduces the historical layout.
SUMMARY_SORT_LABELS = {
    "count_asc": "Peptide count (low → high)",
    "count_desc": "Peptide count (high → low)",
    "alpha_asc": "Alphabetical (A → Z)",
    "alpha_desc": "Alphabetical (Z → A)",
    "median_asc": "Median metric (low → high)",
    "median_desc": "Median metric (high → low)",
}

# Keys whose values are Polars LazyFrames (saved as Parquet)
_LAZYFRAME_KEYS = ("psm_frame", "peaks_frame", "sequence_frame", "heatmap_frame")

# Keys whose values are pickle-safe Python objects
_PICKLE_KEYS = ("annotated_queries", "reference_store", "peptidoform_psms")

# On-disk cache layout version. Bump whenever the pickled object shape changes
# (e.g. PSMs gaining their peak arrays in the PSMList-native migration) so caches
# written by an older code version are cleared and regenerated from inputs rather
# than hydrating semantically-stale objects. Distinct from the runtime
# "cache_version" re-render counter in session state.
# 2: PSM frame gained retention_time and the comparison frame gained rt_delta (#113).
CACHE_SCHEMA_VERSION = 2


def _analysis_cache_path() -> Path:
    """Return the analysis-cache directory inside the active workspace.
    """
    return Path(st.session_state.workspace) / ANALYSIS_CACHE_DIR


def write_analysis_cache(workspace: Path, values: dict, settings: dict | None) -> int:
    """Write one run's results as the workspace's analysis cache; return its cache_version.
    values holds the session-state keys (psm_frame, peaks_frame, ..., comparison_frame,
    psm_lookup, annotated_queries, reference_store, peptidoform_psms). No session state is
    read, so the background job thread can call it. Files go to a staging folder that then
    replaces analysis_cache/, so a reader never sees a half-written cache.
    """
    cache_dir = Path(workspace) / ANALYSIS_CACHE_DIR
    try:
        previous_version = json.loads((cache_dir / "meta.json").read_text()).get("cache_version", 0)
    except (OSError, ValueError):
        previous_version = 0
    run_id = uuid.uuid4().hex
    staging = cache_dir.with_name(f"{ANALYSIS_CACHE_DIR}.staging-{run_id[:8]}")
    # No parents=True: if the workspace was deleted during the run, fail instead of
    # bringing the deleted workspace back with only an analysis cache in it.
    staging.mkdir()
    cache_version = previous_version + 1
    # "settings" and "run_id" are additive (readers that do not know them ignore them), so
    # no schema bump. run_id is unique per run, unlike cache_version, which restarts at 1
    # in every new workspace.
    metadata = {
        "saved_at": datetime.now().isoformat(),
        "cache_version": cache_version,
        "schema_version": CACHE_SCHEMA_VERSION,
        "settings": settings,
        "run_id": run_id,
    }
    try:
        _write_cache_files(staging, values, metadata)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    # Swap: retire the previous cache, move the new one in, then delete the retired one
    retired = cache_dir.with_name(f"{ANALYSIS_CACHE_DIR}.old-{run_id[:8]}")
    if cache_dir.exists():
        os.replace(cache_dir, retired)
    os.replace(staging, cache_dir)
    shutil.rmtree(retired, ignore_errors=True)
    logger.info("Analysis cache saved to %s (version %d)", cache_dir, cache_version)
    return cache_version


def _write_cache_files(staging: Path, values: dict, metadata: dict) -> None:
    """Write every cache file of one run into staging; meta.json last."""
    # Polars LazyFrames -> Parquet (collect each frame to materialise it for write)
    for key in _LAZYFRAME_KEYS:
        frame = values.get(key)
        if frame is not None:
            frame.collect().write_parquet(staging / f"{key}.parquet")

    # comparison_frame (pandas) -> Parquet; write an empty frame as a placeholder for "no metrics"
    comparison_frame = values.get("comparison_frame")
    if comparison_frame is not None and not comparison_frame.empty:
        comparison_frame.to_parquet(staging / "comparison_frame.parquet")
    elif comparison_frame is not None:
        pd.DataFrame().to_parquet(staging / "comparison_frame.parquet")

    # psm_lookup dict (spectrum-id -> PSM metadata) -> JSON
    if values.get("psm_lookup") is not None:
        with open(staging / "psm_lookup.json", "w") as handle:
            json.dump(values["psm_lookup"], handle)

    # Complex Python objects (annotated_queries, reference_store, peptidoform_psms) -> pickle
    for key in _PICKLE_KEYS:
        if values.get(key) is not None:
            with open(staging / f"{key}.pkl", "wb") as handle:
                pickle.dump(values[key], handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Metadata marker, written last. Its mtime keys _read_analysis_cache and tells every
    # session that a new version exists.
    with open(staging / "meta.json", "w") as handle:
        json.dump(metadata, handle)


@st.cache_resource(show_spinner="Restoring analysis from cache...", max_entries=4)
def _read_analysis_cache(cache_dir_str: str, cache_modification_time: int):
    """Read every analysis-cache file from disk. cache_modification_time (meta.json's
    st_mtime_ns) is part of the cache key, so a rewritten cache is read again. It must not
    start with "_": Streamlit leaves such parameters out of the key.
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

    # Pickled objects: annotated_queries, reference_store, peptidoform_psms
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
        modification_time = meta_path.stat().st_mtime_ns
        cached_data = _read_analysis_cache(str(cache_dir), modification_time)

        for key in (*_LAZYFRAME_KEYS, "comparison_frame", "psm_lookup", *_PICKLE_KEYS):
            if key in cached_data:
                st.session_state[get_state_key(key)] = cached_data[key]

        # Flags: mark data ready, record the settings that produced the results, and set
        # cache_version to a token unique to this run. The Explorer component caches and the
        # summary-plot caches key on it, and the integer cache_version restarts at 1 in every
        # workspace. Caches written before run_id existed fall back to saved_at.
        metadata = cached_data.get("_meta", {})
        st.session_state[get_state_key("analysis_settings")] = metadata.get("settings")
        st.session_state[get_state_key("cache_version")] = (
            metadata.get("run_id") or metadata.get("saved_at") or metadata.get("cache_version", 1))
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


# Params that change what the analysis computes, with the label shown when one changed.
ANALYSIS_SETTING_LABELS: dict[str, str] = {
    "selected_idxmls": "idXML files",
    "selected_mzmls": "mzML files",
    "target_peptides": "target peptides",
    "reference_substring": "reference identifier",
    "ion_types": "ion types",
    "mz_tolerance": "m/z tolerance",
    "add_neutral_losses": "neutral losses",
    "intensity_mode": "intensity mode",
    "koina_enabled": "Koina prediction",
    "koina_model": "Koina model",
    "koina_collision_energy": "Koina collision energy",
    "koina_instrument_type": "Koina instrument",
}
# Koina details only matter while Koina is enabled
_KOINA_DETAIL_KEYS = ("koina_model", "koina_collision_energy", "koina_instrument_type")


def analysis_settings(params: dict) -> dict:
    """Snapshot the params that shape the analysis, normalised so that order, case and
    whitespace changes do not count as changes, and JSON-safe for meta.json.
    """
    defaults = {
        "selected_idxmls": [], "selected_mzmls": [], "target_peptides": "",
        "reference_substring": "",
        "ion_types": list(DEFAULTS.ion_types), "mz_tolerance": DEFAULTS.mz_tolerance,
        "add_neutral_losses": DEFAULTS.add_neutral_losses,
        "intensity_mode": DEFAULTS.intensity_mode,
        "koina_enabled": DEFAULTS.koina_enabled, "koina_model": DEFAULTS.koina_model,
        "koina_collision_energy": DEFAULTS.koina_collision_energy,
        "koina_instrument_type": DEFAULTS.koina_instrument_type,
    }
    snapshot = {key: get_param(params, key, default) for key, default in defaults.items()}
    for key in ("selected_idxmls", "selected_mzmls", "ion_types"):
        snapshot[key] = sorted(snapshot[key])
    snapshot["target_peptides"] = sorted(set(parse_peptide_list(snapshot["target_peptides"])))
    snapshot["mz_tolerance"] = float(snapshot["mz_tolerance"])
    if not snapshot["koina_enabled"]:
        for key in _KOINA_DETAIL_KEYS:
            snapshot.pop(key)
    return snapshot


def changed_analysis_settings(params: dict) -> list[str]:
    """Labels of the analysis settings that differ from those of the loaded results.
    Empty when nothing is loaded or the results predate the settings snapshot.
    """
    stored = st.session_state.get(get_state_key("analysis_settings"))
    if not stored or not st.session_state.get(get_state_key("data_ready")):
        return []
    current = analysis_settings(params)
    # A Koina detail counts only when Koina is on in both snapshots; the on/off switch
    # itself is reported as "Koina prediction".
    comparable = [
        key for key in ANALYSIS_SETTING_LABELS
        if key not in _KOINA_DETAIL_KEYS or (key in stored and key in current)
    ]
    return [ANALYSIS_SETTING_LABELS[key] for key in comparable if stored.get(key) != current.get(key)]


def show_stale_results_warning(params: dict, rerun_hint: str) -> None:
    """Warn when the loaded results were made with other settings than the current ones."""
    changed = changed_analysis_settings(params)
    if changed:
        st.warning(
            f"Settings changed since the last run: {', '.join(changed)}. "
            f"The results shown still use the previous settings. {rerun_hint}",
            icon=":material/history:",
        )


# Session keys that hold loaded analysis results (cleared when they belong to another workspace)
_ANALYSIS_SESSION_KEYS = (
    *_LAZYFRAME_KEYS, "comparison_frame", "psm_lookup", *_PICKLE_KEYS,
    "analysis_settings", "cache_version", "data_ready", "analysis_cache_loaded",
)


def restore_analysis_if_needed():
    """Load the workspace's analysis cache into session state unless this session already
    holds the version on disk. A new version appears when a background run finishes (in this
    or another session), so every page picks it up on its next rerun. Costs one stat call.
    Results loaded for another workspace (after creating a workspace or editing the
    ?workspace= URL) are dropped when the current workspace has no cache of its own.
    """
    workspace = str(st.session_state.workspace)
    loaded_key = get_state_key("analysis_cache_loaded")
    loaded = st.session_state.get(loaded_key)
    try:
        on_disk = (workspace, (_analysis_cache_path() / "meta.json").stat().st_mtime_ns)
    except OSError:
        # No cache on disk: keep this workspace's results, drop another workspace's
        if loaded is not None and loaded[0] != workspace:
            for key in _ANALYSIS_SESSION_KEYS:
                st.session_state.pop(get_state_key(key), None)
        return
    if st.session_state.get(get_state_key("data_ready")) and loaded == on_disk:
        return
    if load_analysis_cache():
        st.session_state[loaded_key] = on_disk


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
