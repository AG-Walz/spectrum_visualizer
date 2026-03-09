"""
Unit tests for src/tmpfs_staging.py
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.tmpfs_staging import StagingConfig, TmpfsStaging
from src.parallel_processor import ComparisonTask


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_task(i=0, skip=False):
    """Build a minimal ComparisonTask for staging tests."""
    return ComparisonTask(
        task_id=i,
        task_type="single",
        query_df=pd.DataFrame({"mz": [100.0], "intensity": [1.0], "annotation": ["b1"]}),
        folder_name="PEPTIDEK",
        safe_filename=f"file_{i}",
        png_dir="/fake/png",
        svg_dir="/fake/svg",
        pdf_dir="",
        export_formats=["png"],
        skip_plot=skip,
    )


# ---------------------------------------------------------------------------
# TestStagingDisabled (config.enabled=False)
# ---------------------------------------------------------------------------

class TestStagingDisabled:
    def test_iter_chunks_yields_all_tasks_as_one_chunk(self, tmp_path):
        config = StagingConfig(enabled=False)
        staging = TmpfsStaging(
            final_png_dir=str(tmp_path / "png"),
            final_svg_dir=str(tmp_path / "svg"),
            final_pdf_dir=str(tmp_path / "pdf"),
            export_formats=["png"],
            config=config,
        )
        tasks = [_make_task(i) for i in range(5)]
        with staging:
            chunks = list(staging.iter_chunks(tasks))
        assert len(chunks) == 1
        assert len(chunks[0]) == 5

    def test_context_manager_enters_exits_cleanly(self, tmp_path):
        config = StagingConfig(enabled=False)
        staging = TmpfsStaging(
            final_png_dir=str(tmp_path / "png"),
            final_svg_dir=str(tmp_path / "svg"),
            final_pdf_dir=str(tmp_path / "pdf"),
            export_formats=["png"],
            config=config,
        )
        entered = False
        with staging as s:
            entered = True
            assert s is staging
        assert entered
        assert staging.active is False

    def test_flush_is_noop_when_disabled(self, tmp_path):
        config = StagingConfig(enabled=False)
        staging = TmpfsStaging(
            final_png_dir=str(tmp_path / "png"),
            final_svg_dir=str(tmp_path / "svg"),
            final_pdf_dir=str(tmp_path / "pdf"),
            export_formats=["png"],
            config=config,
        )
        with staging:
            # flush() should not raise even when disabled
            staging.flush()


# ---------------------------------------------------------------------------
# TestStagingEnabled (requires /dev/shm)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not Path("/dev/shm").exists(), reason="no /dev/shm on this system")
class TestStagingEnabled:
    def test_staging_dir_created_in_tmpfs(self, tmp_path):
        config = StagingConfig(
            enabled=True,
            tmpfs_base="/dev/shm",
            safety_margin_gb=0.0,  # No safety margin so tests pass on low-RAM systems
            min_budget_mb=0.0,     # No minimum budget for tests
            ram_fraction=0.01,     # Small fraction to avoid hogging RAM
        )
        staging = TmpfsStaging(
            final_png_dir=str(tmp_path / "png"),
            final_svg_dir="",
            final_pdf_dir="",
            export_formats=["png"],
            config=config,
        )
        with staging as s:
            if s.active:
                # Staging dir should be under /dev/shm
                assert s._staging_dir is not None
                assert str(s._staging_dir).startswith("/dev/shm")
                assert s._staging_dir.exists()

    def test_flush_copies_files_to_final_dir(self, tmp_path):
        final_png = tmp_path / "png"
        final_png.mkdir()
        config = StagingConfig(
            enabled=True,
            tmpfs_base="/dev/shm",
            safety_margin_gb=0.0,
            min_budget_mb=0.0,
            ram_fraction=0.01,
        )
        staging = TmpfsStaging(
            final_png_dir=str(final_png),
            final_svg_dir="",
            final_pdf_dir="",
            export_formats=["png"],
            config=config,
        )
        with staging as s:
            if not s.active:
                pytest.skip("tmpfs staging not active (insufficient RAM budget)")
            # Write a file directly into the staging png directory
            staging_png = Path(s._staging_png_dir)
            test_file = staging_png / "test_file.png"
            test_file.write_text("fake png content")

            # Flush should copy to final_png
            s.flush()

        # After flush, file should be in final dir
        assert (final_png / "test_file.png").exists()
