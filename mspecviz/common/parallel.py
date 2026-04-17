"""Sequential/parallel task dispatch shared by idXML, mzML, and plot workers.

``workers<=1`` runs in-process (debug friendly); otherwise tasks are submitted
to a ``ProcessPoolExecutor`` and drained via ``as_completed``. Tasks are argument
tuples; bare values are auto-wrapped so single-arg workers skip boilerplate.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, List, Optional


def run_tasks(
    tasks: List[Any],
    worker: Callable,
    workers: int = 1,
    on_start: Optional[Callable[[Any], None]] = None,
    on_result: Optional[Callable[[Any, Any], None]] = None,
    on_error: Optional[Callable[[Any, Exception], None]] = None,
    on_progress: Optional[Callable[[int, int], None]] = None,
    stop_when: Optional[Callable[[], bool]] = None,
) -> List[Any]:
    """Run ``worker(*task)`` across ``tasks``; results returned in finish order.

    ``workers`` is clamped to ``[1, len(tasks)]``. Callbacks fire on start,
    result, error (swallowed), and progress. ``stop_when`` predicate, if set,
    cancels pending futures and breaks iteration when it first returns True.
    """
    if not tasks:
        return []

    def _args(task):
        return task if isinstance(task, tuple) else (task,)

    effective_workers = max(1, min(workers, len(tasks)))
    results: List[Any] = []

    if effective_workers <= 1:
        for i, task in enumerate(tasks):
            if on_start:
                on_start(task)
            try:
                result = worker(*_args(task))
                results.append(result)
                if on_result:
                    on_result(task, result)
            except Exception as exc:
                if on_error:
                    on_error(task, exc)
            if on_progress:
                on_progress(i + 1, len(tasks))
            if stop_when and stop_when():
                break
        return results

    with ProcessPoolExecutor(max_workers=effective_workers) as pool:
        future_map = {}
        for task in tasks:
            if on_start:
                on_start(task)
            future_map[pool.submit(worker, *_args(task))] = task

        completed = 0
        for future in as_completed(future_map):
            completed += 1
            task = future_map[future]
            try:
                result = future.result()
                results.append(result)
                if on_result:
                    on_result(task, result)
            except Exception as exc:
                if on_error:
                    on_error(task, exc)
            if on_progress:
                on_progress(completed, len(tasks))
            if stop_when and stop_when():
                for pending in future_map:
                    pending.cancel()
                break

    return results
