"""
Shared text and filename utilities.
"""
import hashlib
import re
from typing import List


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


def parse_peptide_list(text: str) -> List[str]:
    """Parse newline-separated peptide list, uppercased and stripped."""
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
