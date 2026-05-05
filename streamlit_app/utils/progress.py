"""Reusable ``st.status`` + progress-bar + ETA pattern.

The Search/Load/Scan/Export views all opened a ``st.status`` block, made a
``st.progress`` and a ``st.empty()`` caption, computed elapsed/remaining time
the same way, and finalised with ``status.update(label=..., state=...)``. The
context manager here folds that into one entry: ``with progress_status(...) as p:``.
"""
from __future__ import annotations

import time as _time
from contextlib import contextmanager
from typing import Iterator

import streamlit as st

from mspecviz.common.text_utils import format_time


class ProgressTracker:
    """Handle on the live progress bar + caption + status block."""

    def __init__(self, status, progress_bar, status_text, start: float) -> None:
        self._status = status
        self._progress_bar = progress_bar
        self._status_text = status_text
        self._start = start

    @property
    def elapsed(self) -> float:
        return _time.time() - self._start

    def update(self, current: int, total: int, *,
               units: str = "items", message: str | None = None) -> None:
        """Drive the progress bar + caption from a (current, total) tick."""
        if total <= 0:
            return
        fraction = min(current / total, 1.0)
        self._progress_bar.progress(fraction, message)
        if 0 < fraction < 1:
            remaining = (self.elapsed / fraction) - self.elapsed
            self._status_text.caption(
                f"{current}/{total} {units} "
                f"({format_time(self.elapsed)} elapsed, "
                f"~{format_time(remaining)} remaining)"
            )
        elif current >= total:
            self._status_text.caption(f"Done in {format_time(self.elapsed)}")

    def info(self, text: str) -> None:
        """Replace the caption text without touching the progress bar."""
        self._status_text.caption(text)

    def finalize(self, label: str, state: str = "complete") -> None:
        """Mark the status block done; called automatically on context exit if not used."""
        self._progress_bar.progress(1.0)
        self._status.update(label=label, state=state)


@contextmanager
def progress_status(label: str, *, expanded: bool = True) -> Iterator[ProgressTracker]:
    """``with progress_status("Doing X...") as p:`` — yields a :class:`ProgressTracker`."""
    with st.status(label, expanded=expanded) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        yield ProgressTracker(status, progress_bar, status_text, _time.time())
