"""Background analysis jobs: one per workspace, held in this process, outside session state.

A Streamlit script run stops when the user clicks anything or changes page, so a long
analysis cannot run inside it. A job runs in a daemon thread instead; the thread only
updates its AnalysisJob (never st.* or session state) and writes results to the workspace
analysis cache on disk. Every page reads the job for progress and reloads the cache when
the job finishes. Jobs live in process memory: a server restart forgets them, but results
already written to disk stay.

No Streamlit import on purpose, so streamlit_app.utils.state can import this module.
"""
from __future__ import annotations

import itertools
import threading
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

RUNNING, DONE, FAILED, CANCELLED = "running", "done", "failed", "cancelled"


class JobCancelled(Exception):
    """Raised inside the job thread (from a progress callback) when cancel was requested."""


class JobFailed(Exception):
    """An expected failure (e.g. no target peptide found); its message is shown as is,
    without a traceback."""


@dataclass
class AnalysisJob:
    """Progress and outcome of one background run; read by any session, written by its thread."""
    job_id: int
    status: str = RUNNING
    phase: str = "Starting"
    current: int = 0
    total: int = 0
    units: str = ""
    started_at: float = field(default_factory=time.time)
    finished_at: Optional[float] = None
    summary: str = ""
    error: str = ""
    error_traceback: str = ""
    warnings: list[str] = field(default_factory=list)
    cancel_requested: bool = False

    @property
    def is_running(self) -> bool:
        return self.status == RUNNING

    @property
    def elapsed(self) -> float:
        return (self.finished_at or time.time()) - self.started_at

    def progress(self, phase: str, current: int = 0, total: int = 0, units: str = "") -> None:
        """Report progress from the job thread; raises JobCancelled once cancel is requested."""
        self.phase, self.current, self.total, self.units = phase, current, total, units
        if self.cancel_requested:
            raise JobCancelled()

    def warn(self, message: str) -> None:
        self.warnings.append(message)


_JOBS: dict[str, AnalysisJob] = {}
_LOCK = threading.Lock()
_IDS = itertools.count(1)


def job_key_for(workspace) -> str:
    """Registry key of a workspace: its resolved path, so every session agrees on it."""
    return str(Path(workspace).resolve())


def get_job(key: str) -> Optional[AnalysisJob]:
    """Return the latest job for key (a workspace path), running or finished, or None."""
    return _JOBS.get(key)


def start_job(key: str, target: Callable[[AnalysisJob], str]) -> Optional[AnalysisJob]:
    """Start target(job) in a daemon thread; return the new job, or None if one is running.
    target returns the one-line summary shown on success. Any exception marks the job
    failed (with its traceback); JobCancelled marks it cancelled.
    """
    with _LOCK:
        current = _JOBS.get(key)
        if current is not None and current.is_running:
            return None
        job = AnalysisJob(job_id=next(_IDS))
        _JOBS[key] = job

    def run() -> None:
        try:
            job.summary = target(job)
            job.status = DONE
        except JobCancelled:
            job.status = CANCELLED
        except JobFailed as error:
            job.error = str(error)
            job.status = FAILED
        except Exception as error:
            job.error = str(error) or type(error).__name__
            job.error_traceback = traceback.format_exc()
            job.status = FAILED
        finally:
            job.finished_at = time.time()

    threading.Thread(target=run, name=f"analysis-job-{job.job_id}", daemon=True).start()
    return job
