"""Shared text, filename, and input-validation utilities."""
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, List


class DuplicateBasenameError(ValueError):
    """Two or more selected files share the same basename — refuse to start."""


def sanitize_filename(name: str) -> str:
    """
    Removes or replaces invalid characters from filename.
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


def parse_peptide_list(text: str) -> List[str]:
    """Parse newline-separated peptide list, uppercased and stripped.

    Newlines only — callers feeding comma-separated input must pre-normalise
    (e.g. ``text.replace(",", "\\n")`` as ``cli/config.py`` does for the
    inline-YAML-scalar branch).
    """
    return [p.strip().upper() for p in text.split('\n') if p.strip()]


def format_time(seconds: float) -> str:
    """Format seconds as human-readable duration (e.g. '2m 30s')."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def get_nested(source: dict, *keys, default=None):
    """Safely traverse nested dicts: ``get_nested(d, 'a', 'b')`` -> ``d['a']['b']``."""
    current = source
    for key in keys:
        if isinstance(current, dict):
            current = current.get(key)
        else:
            return default
    return current if current is not None else default


def ensure_unique_basenames(paths: Iterable[str], kind: str) -> None:
    """Raise :class:`DuplicateBasenameError` if any two ``paths`` share a basename.

    The pipeline binds idXML PSMs to mzML files by bare filename; if two selected
    files share a basename the binding is undefined and downstream artifacts
    (exports, plots) named after the file would also collide. Catch it at input.
    """
    duplicates = sorted(
        name for name, count in Counter(Path(p).name for p in paths).items()
        if count > 1
    )
    if duplicates:
        raise DuplicateBasenameError(
            f"Two or more selected {kind} files share the same basename: "
            f"{duplicates}. Rename or symlink so every {kind} basename is unique."
        )
