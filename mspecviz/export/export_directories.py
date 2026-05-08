"""
Export utilities: filename sanitization and directory structure creation.
"""
import hashlib
import json
import logging
import re
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def sanitize_filename(name: str) -> str:
    """Strip filesystem-unsafe characters and truncate to 100 chars."""
    sanitized = re.sub(r'[<>:"/\\|?*\[\]]', '_', name)
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_ ')
    return sanitized[:100]


def _settings_fingerprint(run_settings: dict) -> str:
    """Compute a short hash of the settings dict for fast comparison."""
    canonical = json.dumps(run_settings, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def _entry_fingerprint(settings_file: Path) -> str | None:
    """Return the ``_fingerprint`` from a ``run_settings.json``, or None on parse error."""
    try:
        with open(settings_file, encoding="utf-8") as f:
            return json.load(f).get("_fingerprint")
    except (json.JSONDecodeError, OSError) as error:
        logger.warning(f"Could not read {settings_file}: {error}")
        return None


def _scan_existing_dirs(
    base: Path, pattern: re.Pattern, fingerprint: str,
) -> tuple[int, Path | None, int | None]:
    """Walk *base* looking for fingerprint-matching folders + the highest used NNN.

    Returns ``(max_number, match_entry, match_number)``.
    """
    max_number = 0
    match_entry = None
    match_number = None
    if not base.exists():
        return max_number, match_entry, match_number

    for entry in sorted(base.iterdir()):
        if not entry.is_dir():
            continue
        name_match = pattern.match(entry.name)
        if not name_match:
            continue
        number = int(name_match.group(1))
        max_number = max(max_number, number)
        if match_entry is None:
            settings_file = entry / "run_settings.json"
            if settings_file.exists() and _entry_fingerprint(settings_file) == fingerprint:
                match_entry = entry
                match_number = number
    return max_number, match_entry, match_number


def _rename_or_reuse_match(prefix: str, match_entry: Path,
                           match_number: int, today: str) -> Path:
    """Rename a fingerprint-matching folder to today's date, falling back on failure."""
    new_path = match_entry.with_name(f"{prefix}_{match_number:03d}_{today}")
    if match_entry == new_path:
        logger.info(f"Reusing export directory: {new_path.name}")
        return new_path
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
    """Plan an export dir ``{project}_{NNN}_{YYYYMMDD}/`` — reuse on fingerprint match, else bump NNN.

    Returns the root ``Path``; nothing is created until :func:`commit_export_directory`
    runs (so we don't leave empty folders on abort). On re-run with identical
    settings the matching folder is renamed to today's date, enabling
    skip-existing-plots behaviour.
    """
    base = Path(base_dir)
    prefix = sanitize_filename(project_name)
    fingerprint = _settings_fingerprint(run_settings)
    today = datetime.now().strftime("%Y%m%d")
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d{{3}})_(\d{{8}})$")

    max_number, match_entry, match_number = _scan_existing_dirs(base, pattern, fingerprint)

    if match_entry is not None:
        return _rename_or_reuse_match(prefix, match_entry, match_number, today)

    next_number = max_number + 1
    return base / f"{prefix}_{next_number:03d}_{today}"


def commit_export_directory(root_dir: Path, run_settings: dict) -> None:
    """Create export root directory and write run_settings.json."""
    root_dir.mkdir(parents=True, exist_ok=True)

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
