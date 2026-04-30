"""Parallel task dispatch shared by idXML, mzML, and plot workers."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Iterator, Tuple


def run_tasks(
    tasks: list,            # tasks to process; each is a single arg or a tuple of args
    worker: Callable,       # function to run per task
    workers: int = 1,       # number of workers; <=1 runs in-process
) -> Iterator[Tuple[Any, bool, Any]]:
    """Yield ``(task, success, result_or_exc)`` for each task across ``workers`` processes, in finish order.

    ``workers <= 1`` runs in-process. The generator is **lazy**: callers must
    iterate it for any work to happen. Breaking out of the iteration cancels
    pending futures via the ``try/finally`` below.
    """
    if not tasks:
        return

    # Ensure the number of workers is between 1 and the number of tasks
    workers = max(1, min(workers, len(tasks)))

    def args(t):
        return t if isinstance(t, tuple) else (t,)

    # Sequential execution for a single worker
    if workers == 1:
        for task in tasks:
            try:
                yield task, True, worker(*args(task))
            except Exception as exc:
                yield task, False, exc
        return

    # Parallel execution for multiple workers
    with ProcessPoolExecutor(max_workers=workers) as pool:
        future_map = {pool.submit(worker, *args(task)): task for task in tasks}
        try:
            for future in as_completed(future_map):
                task = future_map[future]
                try:
                    yield task, True, future.result()
                except Exception as exc:
                    yield task, False, exc
        finally:
            for f in future_map:
                f.cancel()
