"""
Shared string utilities for filename handling.

Used by src.export and src.summary_plots.
"""
import hashlib
import re


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
