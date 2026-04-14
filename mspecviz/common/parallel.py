"""
Parallel / sequential task dispatch for ProcessPoolExecutor workloads.

Consolidates the repeated pattern where an idXML/mzML/plot worker must run
either in-process (``workers<=1``, easy to debug) or across a process pool
consumed via ``as_completed``. All three former call sites
(``peptide_search.search_peptides_in_idxmls``,
``peptide_search.get_unique_peptides_from_idxmls``,
``spectrum_loader._execute_load_tasks``) now route through ``run_tasks``.

Tasks are *argument tuples*: ``run_tasks`` submits each as ``worker(*task)``.
Single-arg workers can pass plain values — they are auto-wrapped into a
one-tuple so callers don't have to write ``[(x,) for x in xs]`` boilerplate.
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
    """Run *worker* over *tasks* sequentially or in a ``ProcessPoolExecutor``.

    Args:
        tasks: List of argument tuples (or bare values for single-arg workers).
        worker: Module-level function called as ``worker(*task)``.
        workers: Effective workers is ``clamp(workers, 1, len(tasks))``. 1 = in-process.
        on_start: Callback ``(task)`` fired before each worker call (sequential)
            or during submission (parallel). Useful for UI status updates.
        on_result: Callback ``(task, result)`` after each task succeeds.
        on_error: Callback ``(task, exc)`` when a task raises. Exception is swallowed.
        on_progress: Callback ``(completed, total)`` after each task finishes.
        stop_when: Optional predicate checked after each result; when it returns
            True, remaining parallel futures are cancelled and iteration stops.

    Returns:
        List of worker results, in completion order (parallel) or submission
        order (sequential). Tasks that raise are not included in the output.
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
