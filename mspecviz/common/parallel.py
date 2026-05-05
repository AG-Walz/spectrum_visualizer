"""Parallel task dispatch shared by idXML, mzML, and plot workers."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Iterator, Tuple


def run_tasks(
    tasks: list,            # tasks to process; each is a single arg or a tuple of args
    worker: Callable,       # function to run per task
    workers: int = 1,       # number of workers; <=1 runs in-process
) -> Iterator[Tuple[Any, bool, Any]]:
    """Run tasks in parallel (depending on the number of workers),
    yielding (task, success, result or exception) for each.
    """
    if not tasks:
        return

    # Ensure the number of workers is between 1 and the number of tasks
    workers = max(1, min(workers, len(tasks)))

    def _as_args_tuple(task):
        return task if isinstance(task, tuple) else (task,)

    # Sequential execution for a single worker
    if workers == 1:
        for task in tasks:
            try:
                yield task, True, worker(*_as_args_tuple(task))
            except Exception as error:
                yield task, False, error
        return

    # Parallel execution for multiple workers
    with ProcessPoolExecutor(max_workers=workers) as pool:
        future_map = {pool.submit(worker, *_as_args_tuple(task)): task for task in tasks}
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
