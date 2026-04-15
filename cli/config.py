"""
CLI configuration — dataclass + loaders.

Holds only the simple "read from disk" and glob helpers; the three-layer
precedence merge (settings.json < YAML < CLI) and argparse surface live in
``cli.config_builder`` so this file stays focused on I/O.
"""
import glob as globmod
import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import yaml

from mspecviz.common.text_utils import parse_peptide_list

PROJECT_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("mspecviz.cli")


@dataclass
class CLIConfig:
    """All pipeline parameters merged from settings.json, YAML, and CLI args."""

    # Input files (resolved absolute paths after glob expansion)
    idxml_files: List[str] = field(default_factory=list)
    mzml_files: List[str] = field(default_factory=list)

    # Peptides
    peptides: List[str] = field(default_factory=list)
    peptides_mode: str = ""  # "all", "file", or "inline"

    # Reference
    reference_substring: str = ""

    # Annotation
    ion_types: List[str] = field(default_factory=lambda: ["b", "y"])
    mz_tolerance: float = 0.02
    add_losses: bool = False
    add_first_prefix_ion: bool = True
    remove_precursor_peaks: bool = False
    intensity_mode: str = "Relative"

    # Export
    export_dir: str = "./exports"
    project_name: str = "MSpecViz_analysis"
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    dpi: int = 300
    figsize: Tuple[int, int] = (12, 8)
    plot_unpaired: bool = False
    no_plots: bool = False

    # Performance
    workers: int = 4

    # Koina
    koina_enabled: bool = True
    koina_model: str = "ms2pip_timsTOF2024"
    koina_ce: int = 35
    koina_instrument: str = "TIMSTOF"


def load_settings_defaults() -> dict:
    """Load defaults from settings.json."""
    settings_path = PROJECT_ROOT / "settings.json"
    if settings_path.exists():
        with open(settings_path) as f:
            return json.load(f)
    return {}


def load_yaml_config(path: str) -> dict:
    """Load YAML config file, or exit with code 2 on error."""
    config_path = Path(path)
    if not config_path.exists():
        print(f"Error: config file not found: {path}", file=sys.stderr)
        sys.exit(2)
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        print(f"Error: config file is not a valid YAML mapping: {path}", file=sys.stderr)
        sys.exit(2)
    return cfg


def expand_globs(patterns: List[str]) -> List[str]:
    """Expand glob patterns to a sorted list of unique absolute paths."""
    result: List[str] = []
    seen = set()
    for pattern in patterns:
        matches = sorted(globmod.glob(pattern, recursive=True))
        if not matches:
            logger.warning(f"No files matched pattern: {pattern}")
        for m in matches:
            abspath = str(Path(m).resolve())
            if abspath not in seen:
                seen.add(abspath)
                result.append(abspath)
    return result


def load_peptides(
    peptide_spec,
    idxml_files: List[str],
    workers: int,
    quiet: bool,
    progress_factory: Optional[Callable] = None,
) -> Tuple[List[str], str]:
    """
    Resolve peptides from spec (inline list, path, or "all").

    ``progress_factory(label, total, quiet)`` returns a progress indicator when
    scanning idXMLs; defaults to ``cli.output.TerminalProgress`` when None.

    Returns (peptide_list, mode_string).
    """
    if isinstance(peptide_spec, list):
        peps = [str(p).strip().upper() for p in peptide_spec if str(p).strip()]
        return peps, "inline"

    peptide_spec = str(peptide_spec).strip()

    if peptide_spec.lower() == "all":
        # Lazy import to avoid CLI output module churn on bare config loads
        from mspecviz.search import get_unique_peptides_from_idxmls
        if progress_factory is None:
            from cli.output import TerminalProgress as progress_factory  # noqa: N813

        if not quiet:
            print("  Extracting all peptides from idXML files...")
        progress = progress_factory("Scanning", len(idxml_files), quiet)

        def on_progress(current, total, msg):
            progress.update(current, msg)

        peps = get_unique_peptides_from_idxmls(
            idxml_files, limit=10000,
            progress_callback=on_progress, workers=workers,
        )
        progress.finish(f"{len(peps)} unique peptides")
        return peps, "all"

    pep_path = Path(peptide_spec)
    if pep_path.exists() and pep_path.is_file():
        return parse_peptide_list(pep_path.read_text()), "file"

    # Comma/newline-separated inline list
    peps = [p.strip().upper() for p in peptide_spec.replace(",", "\n").split("\n") if p.strip()]
    return peps, "inline"
