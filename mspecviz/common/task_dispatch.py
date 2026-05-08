"""Parallel task dispatcher shared by idXML, mzML, and plot workers.

Wrapper over ``joblib.Parallel(return_as="generator_unordered")``: yields
``(task, ok, result_or_error)`` tuples in completion order. Worker exceptions
are caught per-task, resulting in an early ``break``.
"""

from typing import Any, Callable, Iterator, Sequence, Tuple
from joblib import Parallel, delayed


def _carry(task: Any, worker: Callable) -> Tuple[Any, bool, Any]:
    """Run *worker* on *task*, returning the task alongside the outcome. """
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
    """Dispatch either one or more worker(s) over tasks and yield result or error """
    if not tasks:
        return
    workers = max(1, min(workers, len(tasks)))
    if workers == 1:
        for task in tasks:
            yield _carry(task, worker)
        return
    yield from Parallel(n_jobs=workers, return_as="generator_unordered")(
        delayed(_carry)(task, worker) for task in tasks
    )