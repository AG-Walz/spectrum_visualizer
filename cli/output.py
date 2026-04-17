"""
Terminal-side output helpers for the CLI.

TerminalProgress renders a one-line progress indicator with elapsed + ETA so
callers don't depend on tqdm. The *_summary helpers are pure formatters used
by the pipeline to print the config + completion banners.
"""
import time
from pathlib import Path

from mspecviz.common.text_utils import format_time


class TerminalProgress:
    """One-line progress indicator with elapsed / ETA."""

    def __init__(self, label: str, total: int, quiet: bool = False):
        self.label = label
        self.total = max(total, 1)
        self.current = 0
        self.start = time.time()
        self.quiet = quiet

    def update(self, current: int, message: str = ""):
        self.current = current
        if self.quiet:
            return
        elapsed = time.time() - self.start
        pct = min(self.current / self.total * 100, 100)
        eta_str = ""
        if 0 < self.current < self.total:
            remaining = elapsed / self.current * (self.total - self.current)
            eta_str = f", ~{format_time(remaining)} remaining"
        extra = f"  {message}" if message else ""
        print(f"\r  {self.label} [{self.current}/{self.total}] "
              f"{pct:5.1f}%  ({format_time(elapsed)} elapsed{eta_str}){extra}",
              end="", flush=True)

    def finish(self, message: str = ""):
        self.current = self.total
        if not self.quiet:
            elapsed = time.time() - self.start
            extra = f"  {message}" if message else ""
            print(f"\r  {self.label} [{self.total}/{self.total}] "
                  f"100.0%  ({format_time(elapsed)} elapsed){extra}")


def print_config_summary(config, quiet: bool) -> None:
    """Print the pre-run configuration banner."""
    if quiet:
        return

    print("\n" + "=" * 60)
    print("  MSpecViz CLI Pipeline")
    print("=" * 60)
    print(f"  idXML files:    {len(config.idxml_files)}")
    for f in config.idxml_files:
        print(f"                  {Path(f).name}")
    print(f"  mzML files:     {len(config.mzml_files)}")
    for f in config.mzml_files:
        print(f"                  {Path(f).name}")
    print(f"  Peptides:       {len(config.peptides)} ({config.peptides_mode})")
    if len(config.peptides) <= 10:
        for p in config.peptides:
            print(f"                  {p}")
    else:
        for p in config.peptides[:5]:
            print(f"                  {p}")
        print(f"                  ... and {len(config.peptides) - 5} more")
    print(f"  Reference:      {config.reference_substring or '(none)'}")
    print(f"  Ion types:      {', '.join(config.ion_types)}")
    print(f"  m/z tolerance:  {config.mz_tolerance} Da")
    print(f"  Koina:          {'enabled (' + config.koina_model + ')' if config.koina_enabled else 'disabled'}")
    print(f"  Workers:        {config.workers}")
    print(f"  Export dir:     {config.export_directory}")
    print(f"  Formats:        {', '.join(config.export_formats)}")
    print(f"  Plots:          {'skip' if config.no_plots else 'generate'}")
    print("=" * 60)
    print()


def print_completion_summary(
    *, variant_count: int, total_psms: int, annotated_count: int,
    comparison_count: int, root_dir: Path, elapsed: float, quiet: bool,
) -> None:
    """Print the post-run summary banner."""
    if quiet:
        return

    print()
    print("=" * 60)
    print("  Pipeline Complete")
    print("=" * 60)
    print(f"  Variants:     {variant_count}")
    print(f"  PSMs:         {total_psms}")
    print(f"  Spectra:      {annotated_count}")
    print(f"  Comparisons:  {comparison_count}")
    print(f"  Export dir:   {root_dir}")
    print(f"  Elapsed:      {format_time(elapsed)}")
    print("=" * 60)
    print()
