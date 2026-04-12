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

from src.common.text_utils import sanitize_filename, shorten_stem  # noqa: F401

logger = logging.getLogger(__name__)


def create_export_directory(
    base_dir: str,
    experiment_name: str = None,
    formats: list = None,
) -> Tuple[Path, Path, Path, Path]:
    """Create timestamped export directory with format subdirs. Returns (root, png, svg, pdf)."""
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

    create_all = formats is None
    if create_all or "png" in formats:
        png_dir.mkdir(exist_ok=True)
    if create_all or "svg" in formats:
        svg_dir.mkdir(exist_ok=True)
    if create_all or "pdf" in formats:
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
    Find an existing export directory matching run_settings, or plan a new one.

    Naming: {project_name}_{NNN}_{YYYYMMDD}/
      e.g.  my_project_001_20260226/

    Each folder contains a run_settings.json with the full settings + fingerprint.

    On re-run with identical settings the existing folder is renamed to today's
    date and reused (enables skipping plots already on disk).  If any setting
    differs a new numbered folder is planned.

    **No directories are created by this function.**  Call
    ``commit_export_directory(root_dir, run_settings)`` when you are ready to
    write files (i.e. after confirming there are plots to generate).

    Args:
        base_dir: Parent directory (e.g. ./exports)
        project_name: Project name (e.g. "my_project")
        run_settings: Dict of all settings that affect output

    Returns:
        Tuple of (root_dir, png_dir, svg_dir, pdf_dir) -- Path objects only,
        no directories are created.
    """
    base = Path(base_dir)
    prefix = sanitize_filename(project_name)
    fingerprint = _settings_fingerprint(run_settings)
    today = datetime.now().strftime("%Y%m%d")

    # Pattern: {prefix}_{NNN}_{YYYYMMDD}
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d{{3}})_(\d{{8}})$")

    # Scan existing folders (base may not exist yet on first run)
    max_num = 0
    match_entry = None
    match_num = None

    if base.exists():
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
                logger.info(f"Renamed {match_entry.name} -> {new_name}")
            except OSError as e:
                logger.warning(f"Could not rename folder: {e}")
                new_path = match_entry  # fall back to old path
        else:
            logger.info(f"Reusing export directory: {new_path.name}")

        return new_path, new_path / "png", new_path / "svg", new_path / "pdf"

    # No match -- plan next numbered folder (no mkdir yet)
    next_num = max_num + 1
    folder_name = f"{prefix}_{next_num:03d}_{today}"
    root_dir = base / folder_name
    return root_dir, root_dir / "png", root_dir / "svg", root_dir / "pdf"


def commit_export_directory(root_dir: Path, run_settings: dict) -> None:
    """Create export root directory and write run_settings.json."""
    root_dir.mkdir(parents=True, exist_ok=True)

    settings_file = root_dir / "run_settings.json"
    if settings_file.exists():
        # Existing dir -- update timestamp only
        try:
            with open(settings_file) as f:
                saved = json.load(f)
            saved["_last_run"] = datetime.now().isoformat()
            with open(settings_file, "w") as f:
                json.dump(saved, f, indent=2, default=str)
        except Exception as e:
            logger.warning(f"Could not update run_settings.json: {e}")
    else:
        # New dir -- write full settings
        fingerprint = _settings_fingerprint(run_settings)
        settings_to_save = {
            "_fingerprint": fingerprint,
            "_created_at": datetime.now().isoformat(),
            **run_settings,
        }
        try:
            with open(settings_file, "w") as f:
                json.dump(settings_to_save, f, indent=2, default=str)
        except Exception as e:
            logger.warning(f"Could not write run_settings.json: {e}")
