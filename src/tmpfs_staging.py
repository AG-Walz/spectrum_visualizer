"""
RAM-backed tmpfs staging for plot generation.

Redirects matplotlib savefig() writes to /dev/shm (tmpfs) instead of
network-attached storage, then flushes to disk in bulk.  This eliminates
the small-file random-write bottleneck that dominates large export runs.

Usage:
    from src.tmpfs_staging import TmpfsStaging, StagingConfig

    config = StagingConfig(enabled=True)
    with TmpfsStaging(png_dir, svg_dir, ["svg"], config, log_fn=print) as staging:
        for chunk in staging.iter_chunks(tasks):
            run_parallel_comparisons(tasks=chunk, workers=10, ...)
            staging.flush()

Graceful degradation:
    - No /dev/shm → direct disk (no-op staging)
    - Budget too small → direct disk
    - Config disabled → direct disk
"""
import atexit
import logging
import os
import shutil
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator, List, Optional

logger = logging.getLogger(__name__)


# ============================================================================
# Configuration
# ============================================================================

@dataclass
class StagingConfig:
    """Configuration for tmpfs staging."""

    enabled: bool = True
    tmpfs_base: str = "/dev/shm"
    ram_fraction: float = 0.40
    safety_margin_gb: float = 10.0
    min_budget_mb: float = 50.0

    # Per-file size estimates (bytes) — generous headroom over measured values
    svg_size_estimate: int = 400_000    # ~400 KB (measured mean ~200 KB)
    png_size_estimate: int = 200_000    # ~200 KB


# ============================================================================
# TmpfsStaging context manager
# ============================================================================

class TmpfsStaging:
    """Context manager that optionally stages plot writes through tmpfs.

    When active, rewrites ComparisonTask.png_dir / .svg_dir to point at
    a temporary directory on /dev/shm.  After each chunk, flush() copies
    files to the final disk location and clears the staging area.

    When inactive (fallback), iter_chunks yields all tasks unchanged in a
    single chunk, and flush() is a no-op.
    """

    def __init__(
        self,
        final_png_dir: str,
        final_svg_dir: str,
        final_pdf_dir: str,
        export_formats: List[str],
        config: Optional[StagingConfig] = None,
        log_fn: Optional[Callable[[str], None]] = None,
    ):
        self._final_png_dir = final_png_dir
        self._final_svg_dir = final_svg_dir
        self._final_pdf_dir = final_pdf_dir
        self._export_formats = export_formats
        self._config = config or StagingConfig()
        self._log = log_fn or (lambda msg: None)

        # Computed at __enter__
        self._active = False
        self._staging_dir: Optional[Path] = None
        self._staging_png_dir: Optional[str] = None
        self._staging_svg_dir: Optional[str] = None
        self._staging_pdf_dir: Optional[str] = None
        self._budget_bytes: int = 0
        self._bytes_per_task: int = 0
        self._atexit_registered: bool = False

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def __enter__(self):
        if not self._config.enabled:
            self._log("tmpfs staging: disabled by config")
            return self

        tmpfs_base = Path(self._config.tmpfs_base)
        if not tmpfs_base.is_dir():
            self._log(f"tmpfs staging: {tmpfs_base} not found, using direct disk")
            return self

        # Check available space
        try:
            stat = os.statvfs(str(tmpfs_base))
            available_bytes = stat.f_bavail * stat.f_frsize
        except OSError as e:
            self._log(f"tmpfs staging: cannot stat {tmpfs_base}: {e}")
            return self

        safety_bytes = int(self._config.safety_margin_gb * 1024**3)
        usable_bytes = max(0, available_bytes - safety_bytes)
        budget_bytes = int(usable_bytes * self._config.ram_fraction)

        if budget_bytes < self._config.min_budget_mb * 1024**2:
            self._log(
                f"tmpfs staging: budget too small "
                f"({budget_bytes / 1024**2:.0f} MB < {self._config.min_budget_mb:.0f} MB), "
                f"using direct disk"
            )
            return self

        # Create staging directory
        staging_name = f"mspecviz_{uuid.uuid4().hex[:8]}"
        staging_dir = tmpfs_base / staging_name
        staging_dir.mkdir(parents=True, exist_ok=True)

        self._active = True
        self._staging_dir = staging_dir
        self._budget_bytes = budget_bytes

        # Create format-specific subdirs mirroring the final layout
        if "png" in self._export_formats and self._final_png_dir:
            png_sub = staging_dir / "png"
            png_sub.mkdir(exist_ok=True)
            self._staging_png_dir = str(png_sub)

        if "svg" in self._export_formats and self._final_svg_dir:
            svg_sub = staging_dir / "svg"
            svg_sub.mkdir(exist_ok=True)
            self._staging_svg_dir = str(svg_sub)

        if "pdf" in self._export_formats and self._final_pdf_dir:
            pdf_sub = staging_dir / "pdf"
            pdf_sub.mkdir(exist_ok=True)
            self._staging_pdf_dir = str(pdf_sub)

        # Estimate bytes per task
        bpt = 0
        if "svg" in self._export_formats:
            bpt += self._config.svg_size_estimate
        if "png" in self._export_formats:
            bpt += self._config.png_size_estimate
        if "pdf" in self._export_formats:
            bpt += self._config.png_size_estimate  # PDF ~similar size to PNG
        self._bytes_per_task = max(bpt, 1)

        # Register cleanup for crash / Ctrl+C (only once)
        if not self._atexit_registered:
            atexit.register(self._cleanup)
            self._atexit_registered = True

        avail_gb = available_bytes / 1024**3
        budget_gb = budget_bytes / 1024**3
        self._log(
            f"tmpfs staging: active on {tmpfs_base} "
            f"({avail_gb:.1f} GB available, {budget_gb:.1f} GB budget, "
            f"staging dir: {staging_dir.name})"
        )

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._cleanup()
        return False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def active(self) -> bool:
        """Whether tmpfs staging is active (vs direct-disk fallback)."""
        return self._active

    def iter_chunks(self, tasks: list) -> Iterator[list]:
        """Yield chunks of tasks with paths rewritten to staging dir.

        Tasks with skip_plot=True are NOT rewritten (they produce no files).
        """
        if not self._active:
            # Fallback: single chunk, no path rewriting
            yield tasks
            return

        # Count tasks that will actually write files
        writing_tasks = [t for t in tasks if not t.skip_plot]
        skip_tasks = [t for t in tasks if t.skip_plot]

        if not writing_tasks:
            yield tasks
            return

        # Compute chunk size based on budget
        chunk_size = max(1, self._budget_bytes // self._bytes_per_task)
        total_chunks = (len(writing_tasks) + chunk_size - 1) // chunk_size

        if total_chunks == 1:
            self._log(f"tmpfs staging: all {len(writing_tasks)} tasks fit in one chunk")
        else:
            self._log(
                f"tmpfs staging: splitting {len(writing_tasks)} tasks into "
                f"{total_chunks} chunks of ~{chunk_size} tasks"
            )

        # Yield skip_tasks in the first chunk alongside writing tasks
        for i in range(0, len(writing_tasks), chunk_size):
            chunk_writing = writing_tasks[i:i + chunk_size]

            # Rewrite paths on writing tasks to point at staging dir
            for task in chunk_writing:
                if self._staging_png_dir and task.png_dir:
                    task.png_dir = self._staging_png_dir
                if self._staging_svg_dir and task.svg_dir:
                    task.svg_dir = self._staging_svg_dir
                if self._staging_pdf_dir and task.pdf_dir:
                    task.pdf_dir = self._staging_pdf_dir

            if i == 0:
                # Include skip_tasks (metrics-only) in the first chunk
                yield skip_tasks + chunk_writing
            else:
                yield chunk_writing

    def flush(self):
        """Copy staged files to final disk location, then clear staging."""
        if not self._active or self._staging_dir is None:
            return

        # Copy PNG staging → final
        if self._staging_png_dir and self._final_png_dir:
            src = Path(self._staging_png_dir)
            if src.exists() and any(src.iterdir()):
                shutil.copytree(
                    str(src), self._final_png_dir, dirs_exist_ok=True
                )
                shutil.rmtree(str(src))
                src.mkdir(exist_ok=True)

        # Copy SVG staging → final
        if self._staging_svg_dir and self._final_svg_dir:
            src = Path(self._staging_svg_dir)
            if src.exists() and any(src.iterdir()):
                shutil.copytree(
                    str(src), self._final_svg_dir, dirs_exist_ok=True
                )
                shutil.rmtree(str(src))
                src.mkdir(exist_ok=True)

        # Copy PDF staging → final
        if self._staging_pdf_dir and self._final_pdf_dir:
            src = Path(self._staging_pdf_dir)
            if src.exists() and any(src.iterdir()):
                shutil.copytree(
                    str(src), self._final_pdf_dir, dirs_exist_ok=True
                )
                shutil.rmtree(str(src))
                src.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _cleanup(self):
        """Remove staging directory (idempotent, safe on crash)."""
        if self._staging_dir is not None and self._staging_dir.exists():
            try:
                shutil.rmtree(str(self._staging_dir))
                logger.debug("tmpfs staging: cleaned up %s", self._staging_dir)
            except OSError:
                pass
            self._staging_dir = None
        self._active = False
