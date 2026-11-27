"""Tests for the background analysis job (job_registry) and the session-free cache writer."""
import json
import threading

import pandas as pd
import polars as pl
import pytest

from streamlit_app.utils.job_registry import (
    CANCELLED, DONE, FAILED, JobFailed, get_job, job_key_for, start_job,
)
from streamlit_app.utils.page_state import (
    ANALYSIS_CACHE_DIR, _read_analysis_cache, write_analysis_cache,
)


def _wait(job, timeout=10):
    """Block until the job thread has finished (finished_at is set last)."""
    for _ in range(timeout * 100):
        if job.finished_at is not None:
            return job
        threading.Event().wait(0.01)
    raise AssertionError("job did not finish")


class TestJobRegistry:
    """start_job runs the target in a thread and records the outcome on the job."""

    def test_done_job_keeps_its_summary(self, tmp_path):
        job = _wait(start_job(job_key_for(tmp_path / "a"), lambda job: "all good"))
        assert job.status == DONE and job.summary == "all good"
        assert get_job(job_key_for(tmp_path / "a")) is job

    def test_second_start_is_refused_while_running(self, tmp_path):
        release = threading.Event()
        key = job_key_for(tmp_path / "b")
        first = start_job(key, lambda job: release.wait(5) and "done")
        assert start_job(key, lambda job: "second") is None
        release.set()
        assert _wait(first).status == DONE
        assert start_job(key, lambda job: "again") is not None

    def test_cancel_stops_at_the_next_progress_call(self, tmp_path):
        started = threading.Event()

        def target(job):
            started.set()
            while True:
                job.progress("Working")
                threading.Event().wait(0.01)

        job = start_job(job_key_for(tmp_path / "c"), target)
        started.wait(5)
        job.cancel_requested = True
        assert _wait(job).status == CANCELLED

    def test_expected_and_unexpected_failures(self, tmp_path):
        def expected(job):
            raise JobFailed("nothing found")

        def unexpected(job):
            raise ValueError("boom")

        job = _wait(start_job(job_key_for(tmp_path / "d"), expected))
        assert job.status == FAILED and job.error == "nothing found" and job.error_traceback == ""
        job = _wait(start_job(job_key_for(tmp_path / "e"), unexpected))
        assert job.status == FAILED and "ValueError" in job.error_traceback


class TestWriteAnalysisCache:
    """write_analysis_cache swaps a complete cache in and counts versions up."""

    def _values(self, rows):
        return {
            "psm_frame": pl.DataFrame({"id": list(range(rows))}).lazy(),
            "comparison_frame": pd.DataFrame(),
            "psm_lookup": {"scan=1|PEPTIDE": 0},
            "annotated_queries": {0: ("psm", "frame")},
        }

    def test_reader_returns_the_new_results_after_a_rewrite(self, tmp_path):
        # Regression: the reader's cache key ignored the meta.json mtime, so a session
        # reloading after a background run got the previous run's objects back.
        cache_dir = tmp_path / ANALYSIS_CACHE_DIR
        write_analysis_cache(tmp_path, self._values(1), {"mz_tolerance": 0.02})
        first = _read_analysis_cache(str(cache_dir), (cache_dir / "meta.json").stat().st_mtime_ns)
        write_analysis_cache(tmp_path, self._values(3), {"mz_tolerance": 0.05})
        second = _read_analysis_cache(str(cache_dir), (cache_dir / "meta.json").stat().st_mtime_ns)
        assert first["_meta"]["cache_version"] == 1
        assert second["_meta"]["cache_version"] == 2
        assert second["_meta"]["settings"] == {"mz_tolerance": 0.05}

    def test_failed_write_leaves_no_staging_folder(self, tmp_path):
        values = {**self._values(1), "annotated_queries": {0: lambda: None}}  # not picklable
        with pytest.raises(Exception):
            write_analysis_cache(tmp_path, values, None)
        assert list(tmp_path.iterdir()) == []

    def test_versions_count_up_and_no_staging_folder_remains(self, tmp_path):
        assert write_analysis_cache(tmp_path, self._values(1), {"mz_tolerance": 0.02}) == 1
        assert write_analysis_cache(tmp_path, self._values(3), {"mz_tolerance": 0.05}) == 2
        cache_dir = tmp_path / ANALYSIS_CACHE_DIR
        meta = json.loads((cache_dir / "meta.json").read_text())
        assert meta["cache_version"] == 2 and meta["settings"] == {"mz_tolerance": 0.05}
        assert pl.read_parquet(cache_dir / "psm_frame.parquet").height == 3
        assert [path.name for path in tmp_path.iterdir()] == [ANALYSIS_CACHE_DIR]
