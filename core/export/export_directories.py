"""Export-folder allocation and renaming.
Core helper that decides where each export run writes its plots and names every
folder with a scheme {project}_{NNN}_{YYYYMMDD}/ where NNN is a zero-padded version
number. Each folder carries a run_settings.json recording the settings used, a creation 
and last-run timestamps, so the Export-history view never accumulates duplicate runs.
(called by core.export.export_pipeline and tests.core.export.test_export_directories.)
"""
import hashlib
import json
import logging
import re
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def sanitize_filename(name: str) -> str:
    """Make user provided name safe to use as a filesystem path component.
    """
    # Replace each filesystem-unsafe character with an underscore
    sanitized = re.sub(r'[<>:"/\\|?*\[\]]', '_', name)
    # Collapse consecutive underscores, then strip stray edge underscores/spaces
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_ ')
    # Return cleaned string with a maximum length of 100 characters
    return sanitized[:100]


def _settings_fingerprint(run_settings: dict) -> str:
    """Return a 16-character SHA-256 hash identifying a settings dictionary.
    Identical settings always yield the same hash regardless of key order.
    Used to detect re-runs with unchanged settings.
    """
    # Canonical JSON: sorted keys make the hash order-independent
    canonical = json.dumps(run_settings, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def _entry_fingerprint(settings_file: Path) -> str | None:
    """Read the _fingerprint value stored in a run_settings.json file.
    Returns the fingerprint string, or None if the file is missing.
    """
    try:
        with open(settings_file, encoding="utf-8") as f:
            return json.load(f).get("_fingerprint")
    except (json.JSONDecodeError, OSError) as error:
        logger.warning(f"Could not read {settings_file}: {error}")
        return None


def _scan_existing_dirs(
    base: Path, pattern: re.Pattern, fingerprint: str,
) -> tuple[int, Path | None, int | None]:
    """Scan base for existing export folders matching the project naming scheme.
    Walks the immediate children of the user-defined base folder, considering only directories whose 
    name match the regex pattern {project}_{NNN}_{YYYYMMDD}. Returns max_number (highest version number seen), 
    match_entry and match_number point matching the stored fingerprint.
    """
    max_number = 0
    match_entry = None
    match_number = None
    # A missing base directory means no prior runs — nothing to scan
    if not base.exists():
        return max_number, match_entry, match_number

    # Iterate sorted children so the fingerprint match is deterministic
    for entry in sorted(base.iterdir()):
        if not entry.is_dir():
            continue
        # Skip directories that do not follow the versioned naming scheme
        name_match = pattern.match(entry.name)
        if not name_match:
            continue
        # Track the highest version number for the fallback allocation
        number = int(name_match.group(1))
        max_number = max(max_number, number)
        # Record the first folder whose settings fingerprint matches this run
        if match_entry is None:
            settings_file = entry / "run_settings.json"
            if settings_file.exists() and _entry_fingerprint(settings_file) == fingerprint:
                match_entry = entry
                match_number = number
    return max_number, match_entry, match_number


def _rename_or_reuse_match(prefix: str, match_entry: Path,
                           match_number: int, today: str) -> Path:
    """Re-date a fingerprint-matching export folder to current (today's) date while keeping its version number.
    """
    new_path = match_entry.with_name(f"{prefix}_{match_number:03d}_{today}")
    # Folder is already dated today — nothing to rename
    if match_entry == new_path:
        logger.info(f"Reusing export directory: {new_path.name}")
        return new_path
    # Re-date the folder; on failure keep using the original (rename is not critical)
    try:
        match_entry.rename(new_path)
        logger.info(f"Renamed {match_entry.name} -> {new_path.name}")
        return new_path
    except OSError as error:
        logger.warning(f"Could not rename folder: {error}; reusing original: {match_entry.name}")
        return match_entry


def find_or_create_export_directory(
    base_dir: str,
    project_name: str,
    run_settings: dict,
) -> Path:
    """Decide the export folder for a run.
    Returns a Path named {project}_{NNN}_{YYYYMMDD}. When a prior folder
    carries the same settings fingerprint, that folder is reused (re-dated to today
    via _rename_or_reuse_match); otherwise the next free version number
    max_number + 1 is allocated. This function only plans the path — it does not
    create the folder on disk; call commit_export_directory for that.
    """
    # Setup: sanitize project name, compute fingerprint, build name regex
    base = Path(base_dir)
    prefix = sanitize_filename(project_name)
    fingerprint = _settings_fingerprint(run_settings)
    today = datetime.now().strftime("%Y%m%d")
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d{{3}})_(\d{{8}})$")

    # Scan: highest existing NNN and any fingerprint-matching folder
    max_number, match_entry, match_number = _scan_existing_dirs(base, pattern, fingerprint)

    # Reuse on match, otherwise allocate next NNN
    if match_entry is not None:
        return _rename_or_reuse_match(prefix, match_entry, match_number, today)

    next_number = max_number + 1
    return base / f"{prefix}_{next_number:03d}_{today}"


def commit_export_directory(root_dir: Path, run_settings: dict) -> None:
    """Materialise the export folder on disk and persist its run_settings.json.
    Creates root_dir (and parents) if needed. For a freshly allocated folder it
    writes a new run_settings.json stamped with the settings fingerprint and a
    creation timestamp; for a reused folder it leaves the existing settings intact
    and only refreshes the _last_run timestamp. Returns nothing — I/O failures
    are logged as warnings rather than raised, so a failed write never aborts the export.
    """
    # Ensure the export folder exists
    root_dir.mkdir(parents=True, exist_ok=True)

    # Reused folder: bump _last_run on the existing settings file
    settings_file = root_dir / "run_settings.json"
    if settings_file.exists():
        try:
            with open(settings_file, encoding="utf-8") as f:
                saved = json.load(f)
            saved["_last_run"] = datetime.now().isoformat()
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(saved, f, indent=2, default=str)
        except (json.JSONDecodeError, OSError) as error:
            logger.warning(f"Could not update run_settings.json: {error}")
    # Fresh folder: stamp fingerprint + creation time, write new file
    else:
        settings_to_save = {
            "_fingerprint": _settings_fingerprint(run_settings),
            "_created_at": datetime.now().isoformat(),
            **run_settings,
        }
        try:
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(settings_to_save, f, indent=2, default=str)
        except OSError as error:
            logger.warning(f"Could not write run_settings.json: {error}")
