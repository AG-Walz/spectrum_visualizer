"""Parallel task dispatcher shared by idXML, mzML, and plot workers."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Iterator, Sequence, Tuple


def run_tasks(
    tasks: Sequence[Any],   # tasks to process; each is a single arg or a tuple of args
    worker: Callable,       # function to run per task
    workers: int = 1,       # number of workers; <=1 runs in-process
) -> Iterator[Tuple[Any, bool, Any]]:
    """Run tasks in parallel (or not) through worker(s) and yield ``(task, ok, result_or_error)``.  """
    
    if not tasks:
        return

    # Ensure the number of workers is between 1 and the number of tasks
    workers = max(1, min(workers, len(tasks)))

    # A bare task is auto-wrapped to a 1-tuple of args; to pass a tuple as one
    # positional arg, wrap it explicitly: ``tasks=[(t,)]``.

    # Sequential execution for a single worker
    if workers == 1:
        for task in tasks:
            args = task if isinstance(task, tuple) else (task,)
            try:
                yield task, True, worker(*args)
            except Exception as error:
                yield task, False, error
        return

    # Parallel execution for multiple workers
    with ProcessPoolExecutor(max_workers=workers) as pool:
        future_map = {
            pool.submit(worker, *(task if isinstance(task, tuple) else (task,))): task
            for task in tasks
        }
        try:
            for future in as_completed(future_map):
                task = future_map[future]
                try:
                    yield task, True, future.result()
                except Exception as error:
                    yield task, False, error
        finally:
            for future in future_map:
                future.cancel()
