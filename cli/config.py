"""CLI config dataclass + disk/glob loaders (argparse + merge live in ``cli.config_builder``)."""

import glob as globmod
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Tuple, Union
import yaml

from cli.output import TerminalProgress
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.helpers import parse_peptide_list
from mspecviz.search import get_unique_peptides_from_idxmls

PROJECT_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("mspecviz.cli")


@dataclass
class CLIConfig:
    """Pipeline parameters merged from settings.json, YAML, and CLI args.

    Defaults delegate to :data:`mspecviz.common.defaults.DEFAULTS` — edit
    there, not here, to keep a single source of truth.
    """

    # Input files (resolved absolute paths after glob expansion)
    idxml_files: List[str] = field(default_factory=list)
    mzml_files: List[str] = field(default_factory=list)

    # Peptides — can come from ``peptide_source`` (either a peptides.txt, 
    # or an inline list in YAML). ``load_peptides`` resolves it into the 
    # ``peptides`` list at runtime.
    peptide_source: Union[str, List[str]] = ""
    peptides: List[str] = field(default_factory=list)
    peptides_mode: str = ""                          # "all", "file", or "inline"

    # Reference
    reference_substring: str = ""

    # Annotation
    ion_types: List[str] = field(default_factory=lambda: list(DEFAULTS.ion_types))
    mz_tolerance: float = DEFAULTS.mz_tolerance
    add_neutral_losses: bool = DEFAULTS.add_neutral_losses
    intensity_mode: str = DEFAULTS.intensity_mode

    # Export
    export_directory: str = DEFAULTS.export_directory
    project_name: str = DEFAULTS.project_name
    export_formats: List[str] = field(default_factory=lambda: list(DEFAULTS.export_formats))
    export_resolution: int = DEFAULTS.export_resolution
    figure_size: Tuple[int, int] = DEFAULTS.figure_size
    plot_unpaired_queries: bool = DEFAULTS.plot_unpaired_queries
    no_plots: bool = False

    # Summary plots
    summary_metric: str = DEFAULTS.summary_metric
    summary_plot_types: List[str] = field(
        default_factory=lambda: list(DEFAULTS.summary_plot_types)
    )

    # Performance
    workers: int = DEFAULTS.workers

    # Koina
    koina_enabled: bool = DEFAULTS.koina_enabled
    koina_model: str = DEFAULTS.koina_model
    koina_collision_energy: int = DEFAULTS.koina_collision_energy
    koina_instrument_type: str = DEFAULTS.koina_instrument_type


def load_yaml_config(path: str) -> dict:
    """Load YAML config file, or exit with code 2 on error."""
    config_path = Path(path)
    if not config_path.exists():
        print(f"Error: config file not found: {path}", file=sys.stderr)
        sys.exit(2)
    with open(config_path, encoding="utf-8") as f:
        loaded = yaml.safe_load(f)
    if not isinstance(loaded, dict):
        print(f"Error: config file is not a valid YAML mapping: {path}", file=sys.stderr)
        sys.exit(2)
    return loaded


def expand_globs(patterns: List[str]) -> List[str]:
    """Expand shell-style glob patterns into a flat list of absolute file paths.

    Used to normalise the ``idxml_files`` and ``mzml_files`` entries before
    they land on :class:`CLIConfig`. Accepts literal paths, simple wildcards,
    and recursive ``**`` globs; every pattern is fed through
    :func:`glob.glob(recursive=True) <glob.glob>` so ``**/*.mzML`` walks
    subdirectories.
    """
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
        peptide_source, 
        idxml_files: List[str], 
        workers: int, 
        quiet: bool,
        progress_factory: Optional[Callable] = None,
) -> Tuple[List[str], str]:
    """Resolve a peptide source (inline list / ``"all"`` / file path) into ``(peptides, mode)``.
    
    ``mode`` is recorded in ``run_settings.json`` so reruns are reproducible.
    """
    # --- Option 1: YAML already gave us a Python list (flow-style [...]) ---
    if isinstance(peptide_source, list):
        peptides = [str(p).strip().upper() for p in peptide_source if str(p).strip()]
        return peptides, "inline"

    peptide_source = str(peptide_source).strip()

    # --- Option 2: keyword "all" — scan every idXML for unique peptides ---
    if peptide_source.lower() == "all":
        if progress_factory is None:
            progress_factory = TerminalProgress

        if not quiet:
            print("  Extracting all peptides from idXML files...")
        progress = progress_factory("Scanning", len(idxml_files), quiet)

        def on_progress(current, total, message):
            progress.update(current, message)

        peptides = get_unique_peptides_from_idxmls(
            idxml_files, limit=10000,
            progress_callback=on_progress, workers=workers,
        )
        progress.finish(f"{len(peptides)} unique peptides")
        return peptides, "all"

    # --- Option 3: a text file with comma/newline-separated sequences ---
    peptide_path = Path(peptide_source)
    if peptide_path.exists() and peptide_path.is_file():
        return parse_peptide_list(peptide_path.read_text()), "file"

    # --- Option 4: bare inline string, e.g. "PEPTIDE, SAMPLE" from YAML scalar ---
    peptides = [p.strip().upper() for p in peptide_source.replace(",", "\n").split("\n") if p.strip()]
    return peptides, "inline"
