"""Shared utility functions used by both CLI and Streamlit UI."""
from typing import List


def parse_peptide_list(text: str) -> List[str]:
    """Parse newline-separated peptide sequences into a normalized list.

    Strips whitespace, removes blank lines, and converts to uppercase.
    Use ``sorted(parse_peptide_list(text))`` when stable ordering is needed.
    """
    return [p.strip().upper() for p in text.split('\n') if p.strip()]


def format_time(seconds: float) -> str:
    """Format seconds as human-readable time."""
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


def get_nested(d: dict, *keys, default=None):
    """Safely get nested dict value."""
    val = d
    for k in keys:
        if isinstance(val, dict):
            val = val.get(k)
        else:
            return default
    return val if val is not None else default
