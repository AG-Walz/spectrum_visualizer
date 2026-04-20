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

from mspecviz.common.text_utils import sanitize_filename

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
    """Plan an export dir ``{project}_{NNN}_{YYYYMMDD}/`` — reuse on fingerprint match, else bump NNN.

    Returns ``(root, png, svg, pdf)`` as plain Paths; nothing is created until
    :func:`commit_export_directory` runs (so we don't leave empty folders on
    abort). On re-run with identical settings the matching folder is renamed
    to today's date, enabling skip-existing-plots behavior.
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
            name_match = pattern.match(entry.name)
            if not name_match:
                continue
            num = int(name_match.group(1))
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
                    except (json.JSONDecodeError, OSError) as e:
                        logger.warning(f"Could not read {settings_file}: {e}")

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
        except (json.JSONDecodeError, OSError) as e:
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
        except OSError as e:
            logger.warning(f"Could not write run_settings.json: {e}")
