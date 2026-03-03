"""
Export utilities: filename sanitization and directory structure creation.
"""
import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Tuple
from datetime import datetime

logger = logging.getLogger(__name__)


def sanitize_filename(name: str) -> str:
    """
    Sanitize a string for use as a filename.

    Removes or replaces invalid characters.
    """
    # Replace problematic characters
    sanitized = re.sub(r'[<>:"/\\|?*\[\]]', '_', name)
    # Replace multiple underscores with single
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores and spaces
    sanitized = sanitized.strip('_ ')
    # Truncate if too long
    if len(sanitized) > 100:
        sanitized = sanitized[:100]
    return sanitized


def shorten_stem(stem: str, max_len: int = 30) -> str:
    """Shorten a file stem to *max_len* chars, appending a hash if truncated.

    If the stem already fits, it is returned unchanged.  Otherwise the
    result is ``stem[:max_len-9] + '_' + hash8`` where *hash8* is the
    first 8 hex chars of the SHA-256 of the full stem.  This keeps
    filenames short while guaranteeing uniqueness.
    """
    if len(stem) <= max_len:
        return stem
    digest = hashlib.sha256(stem.encode("utf-8")).hexdigest()[:8]
    return f"{stem[:max_len - 9]}_{digest}"


def create_export_directory(
    base_dir: str,
    experiment_name: str = None
) -> Tuple[Path, Path, Path]:
    """
    Create export directory structure.

    Args:
        base_dir: Base output directory
        experiment_name: Optional name for the experiment

    Returns:
        Tuple of (root_dir, png_dir, svg_dir)
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if experiment_name:
        folder_name = f"{sanitize_filename(experiment_name)}_{timestamp}"
    else:
        folder_name = f"export_{timestamp}"

    root_dir = Path(base_dir) / folder_name
    png_dir = root_dir / "png"
    svg_dir = root_dir / "svg"
    pdf_dir = root_dir / "pdf"

    root_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(exist_ok=True)
    svg_dir.mkdir(exist_ok=True)
    pdf_dir.mkdir(exist_ok=True)

    return root_dir, png_dir, svg_dir, pdf_dir


def _settings_fingerprint(run_settings: dict) -> str:
    """Compute a short hash of the settings dict for fast comparison."""
    canonical = json.dumps(run_settings, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def find_or_create_export_directory(
    base_dir: str,
    project_name: str,
    run_settings: dict,
) -> Tuple[Path, Path, Path, Path]:
    """
    Find an existing export directory matching run_settings, or create a new one.

    Naming: {project_name}_{NNN}_{YYYYMMDD}/
      e.g.  ProCa_2_run_001_20260226/

    Each folder contains a run_settings.json with the full settings + fingerprint.

    On re-run with identical settings the existing folder is renamed to today's
    date and reused (enables skipping plots already on disk).  If any setting
    differs a new numbered folder is created.

    Args:
        base_dir: Parent directory (e.g. ./exports)
        project_name: Project name (e.g. "ProCa_2_run")
        run_settings: Dict of all settings that affect output

    Returns:
        Tuple of (root_dir, png_dir, svg_dir)
    """
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)
    prefix = sanitize_filename(project_name)
    fingerprint = _settings_fingerprint(run_settings)
    today = datetime.now().strftime("%Y%m%d")

    # Pattern: {prefix}_{NNN}_{YYYYMMDD}
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d{{3}})_(\d{{8}})$")

    # Scan existing folders
    max_num = 0
    match_entry = None
    match_num = None

    for entry in sorted(base.iterdir()):
        if not entry.is_dir():
            continue
        m = pattern.match(entry.name)
        if not m:
            continue
        num = int(m.group(1))
        max_num = max(max_num, num)

        # Check settings match
        if match_entry is None:
            settings_file = entry / "run_settings.json"
            if settings_file.exists():
                try:
                    with open(settings_file) as f:
                        saved = json.load(f)
                    if saved.get("_fingerprint") == fingerprint:
                        match_entry = entry
                        match_num = num
                except Exception:
                    pass

    if match_entry is not None:
        # Rename to today's date and reuse
        new_name = f"{prefix}_{match_num:03d}_{today}"
        new_path = base / new_name
        if match_entry != new_path:
            try:
                match_entry.rename(new_path)
                logger.info(f"Renamed {match_entry.name} → {new_name}")
            except OSError as e:
                logger.warning(f"Could not rename folder: {e}")
                new_path = match_entry  # fall back to old path
        else:
            logger.info(f"Reusing export directory: {new_path.name}")

        png_dir = new_path / "png"
        svg_dir = new_path / "svg"
        pdf_dir = new_path / "pdf"
        png_dir.mkdir(exist_ok=True)
        svg_dir.mkdir(exist_ok=True)
        pdf_dir.mkdir(exist_ok=True)

        # Update timestamp in run_settings.json
        settings_file = new_path / "run_settings.json"
        try:
            with open(settings_file) as f:
                saved = json.load(f)
            saved["_last_run"] = datetime.now().isoformat()
            with open(settings_file, "w") as f:
                json.dump(saved, f, indent=2, default=str)
        except Exception:
            pass

        return new_path, png_dir, svg_dir, pdf_dir

    # No match — create next numbered folder
    next_num = max_num + 1
    folder_name = f"{prefix}_{next_num:03d}_{today}"
    root_dir = base / folder_name
    png_dir = root_dir / "png"
    svg_dir = root_dir / "svg"
    pdf_dir = root_dir / "pdf"

    root_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(exist_ok=True)
    svg_dir.mkdir(exist_ok=True)
    pdf_dir.mkdir(exist_ok=True)

    # Write run_settings.json
    settings_to_save = {
        "_fingerprint": fingerprint,
        "_created_at": datetime.now().isoformat(),
        **run_settings,
    }
    try:
        with open(root_dir / "run_settings.json", "w") as f:
            json.dump(settings_to_save, f, indent=2, default=str)
    except Exception as e:
        logger.warning(f"Could not write run_settings.json: {e}")

    return root_dir, png_dir, svg_dir, pdf_dir


