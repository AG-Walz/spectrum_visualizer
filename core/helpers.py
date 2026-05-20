"""Shared text-parsing, time-formatting, input-validation, and parallel-dispatch utilities."""

from collections import Counter
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, List, Sequence, Tuple
from joblib import Parallel, delayed


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


def _carry(task: Any, worker: Callable) -> Tuple[Any, bool, Any]:
    """Wrapper that catches exceptions from a worker and returns them as values 
    ``(task, ok, result_or_error)`` with ``ok`` as a boolean flag."""
    args = task if isinstance(task, tuple) else (task,)
    try:
        return task, True, worker(*args)
    except Exception as error:
        return task, False, error


def run_tasks(
    tasks: Sequence[Any],
    worker: Callable,
    workers: int = 1,
) -> Iterator[Tuple[Any, bool, Any]]:
    """Utility that dispatches a function (worker) across a list of tasks,
     either one at a time (workers = 1) or in parallel across multiple processes;
     returns ``(task, ok, result_or_error)`` tuples."""

    # If empty input return nothing
    if not tasks:
        return
    # Cap worker count: never spawn more processes than there are tasks
    workers = max(1, min(workers, len(tasks)))

    # Sequential path: keep tracebacks readable when debugging a single worker
    if workers == 1:
        for task in tasks:
            yield _carry(task, worker)
        return

    # Parallel path: joblib streams results as soon as each worker finishes
    yield from Parallel(n_jobs=workers, return_as="generator_unordered")(
        delayed(_carry)(task, worker) for task in tasks
    )