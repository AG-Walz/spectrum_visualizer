"""
Export utilities: filename sanitization and directory structure creation.
"""
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
        folder_name = f"MSpecViz_{sanitize_filename(experiment_name)}_{timestamp}"
    else:
        folder_name = f"MSpecViz_export_{timestamp}"

    root_dir = Path(base_dir) / folder_name
    png_dir = root_dir / "png"
    svg_dir = root_dir / "svg"

    root_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(exist_ok=True)
    svg_dir.mkdir(exist_ok=True)

    return root_dir, png_dir, svg_dir


