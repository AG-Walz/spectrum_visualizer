"""One-line st.status + progress-bar + ETA helper used by Load & Annotate, Explore, and Export.
"""
from __future__ import annotations

import time as _time
from contextlib import contextmanager
from typing import Iterator

import streamlit as st

from core.helpers import format_time


class ProgressTracker:
    """Handle on the live st.status block, progress bar, and caption.
    """
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
        """Advance the bar and refresh the elapsed/remaining caption from a (current, total) tick.
        """
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
        """Overwrite the caption text without moving the progress bar.
        """
        self._status_text.caption(text)

    def finalize(self, label: str, state: str = "complete") -> None:
        """Fill the bar and mark the st.status block done.
        """
        self._progress_bar.progress(1.0)
        self._status.update(label=label, state=state)


@contextmanager
def progress_status(label: str, *, expanded: bool = True) -> Iterator[ProgressTracker]:
    """Open an st.status block and yield a ProgressTracker for with blocks.
    """
    with st.status(label, expanded=expanded) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        yield ProgressTracker(status, progress_bar, status_text, _time.time())
