"""On-disk idXML index: parse each idXML once, keep a slim per-file parquet, reuse it.

Parsing an idXML into psm_utils PSM objects costs about 0.5 ms and 11 kB of RAM per PSM,
and the target peptides are usually a tiny share of all PSMs. The index stores two
columns per PSM — the uppercased bare sequence and the pickled PSM — so:

- peptide counts read only the sequence column (no PSM objects are built),
- a target lookup unpickles only the matching rows,
- a file is parsed again only when its path, size or modification time changes.

Index files live in one shared cache directory, so workspaces that select the same
idXML share one index; rebuilding an index after its idXML changed deletes the old one.
Parse errors are stored next to the index file and reported again on every use, exactly
as a fresh parse would report them.

Called by streamlit_app.views.widgets.peptide_input (peptide counts) and
streamlit_app.views.load_annotate.load_annotate_actions (target lookup).
"""
import hashlib
import json
import os
import pickle
import uuid
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable, List, Optional

import polars as pl
from psm_utils import PSMList

from core.helpers import run_tasks
from core.peptide_lookup import _load_idxml_worker

# Bump when the index layout or the PSM fields stamped by _load_idxml_worker change,
# so stale index files are ignored and rebuilt rather than read.
INDEX_VERSION = 1

_INDEX_SCHEMA = {"sequence": pl.String, "psm": pl.Binary}


def _path_prefix(idxml_path: str) -> str:
    """First half of an index file name: identifies the idXML path, whatever its content."""
    return hashlib.sha1(idxml_path.encode()).hexdigest()[:20]


def index_file_for(idxml_path: str, cache_dir: Path) -> Path:
    """Return the index file path for idxml_path; raises OSError when the idXML is unreadable.
    Named <path hash>-<content hash>: the content hash covers size, modification time and
    INDEX_VERSION, so any change to the file (or to the index layout) points at a new,
    not-yet-built index file, and the path hash lets the builder delete the old one.
    """
    stat = os.stat(idxml_path)
    content_key = f"{stat.st_size}|{stat.st_mtime_ns}|v{INDEX_VERSION}"
    content_hash = hashlib.sha1(content_key.encode()).hexdigest()[:20]
    return Path(cache_dir) / f"{_path_prefix(idxml_path)}-{content_hash}.parquet"


def _errors_file(index_file: Path) -> Path:
    """Sidecar JSON that holds the parse errors of one indexed idXML (absent when clean)."""
    return index_file.with_suffix(".errors.json")


def _build_index_worker(idxml_path: str, index_file: str) -> int:
    """Per-file worker: parse one idXML and write its index; returns the PSM count.
    Partial results from a parse error are indexed too, with the errors in a sidecar file.
    Both files are written under a temporary name and renamed, so a crash never leaves a
    half-written index that later looks complete.
    """
    # Parse once with the same loader (and provenance stamps) as core.peptide_lookup
    _, psms, errors = _load_idxml_worker(idxml_path, Path(idxml_path).name)
    index_path = Path(index_file)
    # Unique per call: two threads of one server process may build the same index at once
    # (inline when workers=1), and a shared name made the second rename fail.
    staging_suffix = f".{uuid.uuid4().hex}.tmp"

    # Errors first: the index file is the "done" marker, so it must be renamed last
    if errors:
        errors_path = _errors_file(index_path)
        errors_staging = Path(f"{errors_path}{staging_suffix}")
        errors_staging.write_text(json.dumps(errors))
        os.replace(errors_staging, errors_path)

    frame = pl.DataFrame(
        {
            "sequence": [(psm.peptidoform.sequence or "").upper() for psm in psms],
            "psm": [pickle.dumps(psm, protocol=pickle.HIGHEST_PROTOCOL) for psm in psms],
        },
        schema=_INDEX_SCHEMA,
    )
    staging = Path(f"{index_file}{staging_suffix}")
    frame.write_parquet(staging)
    os.replace(staging, index_path)

    # Keep one index per idXML path: drop the files of its earlier versions
    keep = {index_path.name, _errors_file(index_path).name}
    for old_file in index_path.parent.glob(f"{_path_prefix(idxml_path)}-*"):
        if old_file.name not in keep and not old_file.name.endswith(".tmp"):
            old_file.unlink(missing_ok=True)
    return len(psms)


def missing_from_index(idxml_paths: Iterable[str], cache_dir: Path) -> List[str]:
    """Return the idXML paths that have no up-to-date index file (unreadable paths included)."""
    missing = []
    for path in idxml_paths:
        try:
            if not index_file_for(path, cache_dir).exists():
                missing.append(path)
        except OSError:
            missing.append(path)
    return missing


def ensure_idxml_index(
    idxml_paths: Iterable[str],
    cache_dir: Path,
    workers: int = 1,
    progress_callback: Optional[Callable] = None,
    error_callback: Optional[Callable] = None,
) -> List[Path]:
    """Return the index file of every readable idXML in idxml_paths, building missing ones.
    Only files without an up-to-date index are parsed (in parallel when workers > 1).
    progress_callback(done, total, message) counts the files being built; error_callback
    receives every stored or new parse error, plus one message per unreadable path.
    The returned list keeps the order of idxml_paths, so lookups are deterministic.
    """
    progress = progress_callback or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Resolve: map each readable path to its index file; queue the ones not built yet
    index_files: dict[str, Path] = {}
    tasks: list[tuple[str, str]] = []
    for path in dict.fromkeys(idxml_paths):
        try:
            index_file = index_file_for(path, cache_dir)
        except OSError as error:
            on_error(f"Cannot read {Path(path).name}: {error}")
            continue
        index_files[path] = index_file
        if not index_file.exists():
            tasks.append((path, str(index_file)))

    # Build: parse each missing idXML once; a crashed worker drops only its own file
    for done, (task, ok, result) in enumerate(
        run_tasks(tasks=tasks, worker=_build_index_worker, workers=workers), start=1,
    ):
        if not ok:
            on_error(f"Worker crashed indexing {Path(task[0]).name}: {result}")
            index_files.pop(task[0], None)
        progress(done, len(tasks), "Indexing idXML files")

    # Report: stored parse errors, for fresh and cached index files alike
    for index_file in index_files.values():
        errors_path = _errors_file(index_file)
        if errors_path.exists():
            for message in json.loads(errors_path.read_text()):
                on_error(message)
    return list(index_files.values())


def count_indexed_peptides(index_files: List[Path]) -> Counter:
    """Count PSMs per uppercased bare sequence across index_files, without building PSMs.
    Same result as core.peptide_lookup.count_psms_per_peptide on the parsed files.
    """
    if not index_files:
        return Counter()
    counts = (
        pl.scan_parquet([str(path) for path in index_files])
        .filter(pl.col("sequence") != "")
        .group_by("sequence")
        .len()
        .collect()
    )
    return Counter(dict(zip(counts["sequence"].to_list(), counts["len"].to_list())))


def load_indexed_psms(index_files: List[Path], target_peptides: Iterable[str]) -> PSMList:
    """Return the PSMs whose bare sequence is a target peptide, unpickling only those rows.
    The result is ready for core.peptide_lookup.lookup_peptides (grouping + mzML resolution).
    """
    targets = sorted({peptide.strip().upper() for peptide in target_peptides if peptide.strip()})
    if not index_files or not targets:
        return PSMList(psm_list=[])
    hits = (
        pl.scan_parquet([str(path) for path in index_files])
        .filter(pl.col("sequence").is_in(targets))
        .select("psm")
        .collect()
    )
    return PSMList(psm_list=[pickle.loads(blob) for blob in hits["psm"].to_list()])
