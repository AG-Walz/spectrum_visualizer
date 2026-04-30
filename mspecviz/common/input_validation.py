"""Input-side guards shared by the CLI and the Streamlit app."""
from collections import Counter
from pathlib import Path
from typing import Iterable


class DuplicateBasenameError(ValueError):
    """Two or more selected files share the same basename — refuse to start."""


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
