"""Parallel task dispatcher shared by ``search.idxml_workers``, 
``spectra.spectrum_streaming``, and ``plotting.plot_worker``."""

from typing import Any, Callable, Iterator, Sequence, Tuple
from joblib import Parallel, delayed


def _carry(task: Any, worker: Callable) -> Tuple[Any, bool, Any]:
    """Run a worker on one task and return ``(task, ok, result_or_error)`` 
    so failures never crash the batch."""
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
    """Dispatch ``worker`` across ``tasks`` and stream ``(task, ok, result_or_error)`` tuples; 
    ``workers=1`` runs sequentially for easier debugging."""
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