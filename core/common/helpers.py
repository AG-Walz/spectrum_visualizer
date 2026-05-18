"""Shared text-parsing, time-formatting, and input-validation utilities."""
from collections import Counter
from pathlib import Path
from typing import Iterable, List

class DuplicateBasenameError(ValueError):
    """Two or more selected files share the same basename — refuse to start."""


def ensure_unique_basenames(paths: Iterable[str], kind: str) -> None:
    """Raise :class:`DuplicateBasenameError` if any two ``paths`` share a basename."""
    duplicates = sorted(
        name for name, count in Counter(Path(path).name for path in paths).items()
        if count > 1
    )
    if duplicates:
        raise DuplicateBasenameError(
            f"Two or more selected {kind} files share the same basename: "
            f"{duplicates}. Rename or symlink so every {kind} basename is unique."
        )

def parse_peptide_list(text: str) -> List[str]:
    """Parse newline-separated peptide list, uppercased and stripped."""
    return [line.strip().upper() for line in text.splitlines() if line.strip()]


def format_time(seconds: float) -> str:
    """Format seconds as human-readable duration (e.g. ``"2m 30s"``)."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        minutes = int(seconds // 60)
        remaining_seconds = int(seconds % 60)
        return f"{minutes}m {remaining_seconds}s"
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    return f"{hours}h {minutes}m"

