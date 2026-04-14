"""
CLI configuration loading + argument parsing.

Three-layer merge: settings.json defaults < YAML config < CLI flags. The table
at the top of ``build_config`` documents the full precedence chain so adding a
new parameter only requires one table entry instead of four scattered branches.
"""
import argparse
import glob as globmod
import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import yaml

from mspecviz.common.text_utils import get_nested, parse_peptide_list

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
        from mspecviz.peptide_search import get_unique_peptides_from_idxmls
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


# --- YAML → CLIConfig mapping tables ---------------------------------------
# Each entry: (yaml_section, yaml_key, cfg_attr, optional_cast)
_YAML_OVERRIDES = [
    ("annotation", "ion_types", "ion_types", None),
    ("annotation", "mz_tolerance", "mz_tolerance", float),
    ("annotation", "add_losses", "add_losses", bool),
    ("annotation", "add_first_prefix_ion", "add_first_prefix_ion", bool),
    ("annotation", "remove_precursor_peaks", "remove_precursor_peaks", bool),
    ("annotation", "intensity_mode", "intensity_mode", None),
    ("export", "directory", "export_dir", None),
    ("export", "project_name", "project_name", None),
    ("export", "formats", "export_formats", None),
    ("export", "dpi", "dpi", int),
    ("export", "figsize", "figsize", tuple),
    ("export", "plot_unpaired", "plot_unpaired", bool),
    ("performance", "workers", "workers", int),
    ("koina", "enabled", "koina_enabled", bool),
    ("koina", "model", "koina_model", None),
    ("koina", "collision_energy", "koina_ce", int),
    ("koina", "instrument_type", "koina_instrument", None),
]

# Each entry: (args_attr, cfg_attr, optional_transform)
# Transforms run only when the attr has a meaningful value (not None / falsy flag).
_CLI_OVERRIDES = [
    ("idxml", "idxml_files", expand_globs),
    ("mzml", "mzml_files", expand_globs),
    ("reference_substring", "reference_substring", None),  # None sentinel → skip
    ("ion_types", "ion_types", None),
    ("mz_tolerance", "mz_tolerance", None),
    ("export_dir", "export_dir", None),
    ("project_name", "project_name", None),
    ("formats", "export_formats", None),
    ("workers", "workers", None),
    ("koina_model", "koina_model", None),
    ("koina_ce", "koina_ce", None),
    ("koina_instrument", "koina_instrument", None),
]


def _apply_settings_defaults(cfg: CLIConfig, defaults: dict) -> None:
    """Populate cfg from settings.json (layer 1)."""
    cfg.ion_types = get_nested(defaults, "ion_types", "default", default=["b", "y"])
    cfg.mz_tolerance = get_nested(defaults, "mz_tolerance", "default", default=0.02)
    cfg.add_losses = get_nested(defaults, "annotation", "add_losses", default=False)
    cfg.add_first_prefix_ion = get_nested(defaults, "annotation", "add_first_prefix_ion", default=True)
    cfg.remove_precursor_peaks = get_nested(defaults, "annotation", "remove_precursor_peaks", default=False)
    cfg.export_dir = get_nested(defaults, "export", "default_directory", default="./exports")
    cfg.project_name = get_nested(defaults, "export", "project_name", default="MSpecViz_analysis")
    cfg.export_formats = get_nested(defaults, "export", "formats", default=["png"])
    cfg.dpi = get_nested(defaults, "export", "dpi", default=300)
    figsize = get_nested(defaults, "export", "figsize", default=[12, 8])
    cfg.figsize = tuple(figsize) if isinstance(figsize, list) else (12, 8)
    cfg.workers = get_nested(defaults, "performance", "workers", default=4)
    cfg.koina_enabled = get_nested(defaults, "koina", "enabled", default=True)
    cfg.koina_model = get_nested(defaults, "koina", "model", default="ms2pip_timsTOF2024")
    cfg.koina_ce = get_nested(defaults, "koina", "collision_energy", default=35)
    cfg.koina_instrument = get_nested(defaults, "koina", "instrument_type", default="TIMSTOF")


def _apply_yaml_overrides(cfg: CLIConfig, yaml_cfg: dict) -> None:
    """Apply YAML config overrides (layer 2)."""
    cfg.idxml_files = expand_globs(get_nested(yaml_cfg, "input", "idxml_files", default=[]))
    cfg.mzml_files = expand_globs(get_nested(yaml_cfg, "input", "mzml_files", default=[]))
    cfg.reference_substring = get_nested(yaml_cfg, "reference", "substring", default="")

    for section_name, yaml_key, attr, cast in _YAML_OVERRIDES:
        section = yaml_cfg.get(section_name, {})
        if yaml_key in section:
            val = cast(section[yaml_key]) if cast else section[yaml_key]
            setattr(cfg, attr, val)


def _apply_cli_overrides(cfg: CLIConfig, args: argparse.Namespace) -> None:
    """Apply argparse CLI overrides (layer 3, highest priority)."""
    for arg_attr, cfg_attr, transform in _CLI_OVERRIDES:
        val = getattr(args, arg_attr, None)
        if val is None or val == []:
            continue
        setattr(cfg, cfg_attr, transform(val) if transform else val)

    # Boolean flags with non-uniform semantics
    if args.no_plots:
        cfg.no_plots = True
    if args.no_koina:
        cfg.koina_enabled = False


def build_config(yaml_cfg: dict, args: argparse.Namespace) -> CLIConfig:
    """Merge settings.json defaults, YAML config, and CLI overrides."""
    cfg = CLIConfig()
    _apply_settings_defaults(cfg, load_settings_defaults())
    _apply_yaml_overrides(cfg, yaml_cfg)
    _apply_cli_overrides(cfg, args)
    return cfg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mspecviz-cli",
        description="MSpecViz CLI -- Terminal-based MS/MS spectrum analysis pipeline.",
        epilog="Example: python -m cli --config cli/config/config_example.yaml --verbose",
    )

    parser.add_argument("-c", "--config", required=True,
                        help="Path to YAML config file")

    parser.add_argument("--idxml", nargs="+", metavar="PATTERN",
                        help="idXML file paths or glob patterns (overrides config)")
    parser.add_argument("--mzml", nargs="+", metavar="PATTERN",
                        help="mzML file paths or glob patterns (overrides config)")
    parser.add_argument("--peptides", metavar="PATH_OR_ALL",
                        help='Path to peptides.txt or "all" (overrides config)')
    parser.add_argument("--reference-substring", metavar="STR",
                        help="Reference file substring (overrides config)")

    parser.add_argument("--ion-types", nargs="+", metavar="ION",
                        help="Ion types, e.g. b y (overrides config)")
    parser.add_argument("--mz-tolerance", type=float, metavar="DA",
                        help="m/z tolerance in Da (overrides config)")

    parser.add_argument("--export-dir", metavar="DIR",
                        help="Export directory (overrides config)")
    parser.add_argument("--project-name", metavar="NAME",
                        help="Project name (overrides config)")
    parser.add_argument("--formats", nargs="+", metavar="FMT",
                        choices=["png", "svg", "pdf", "html"],
                        help="Export formats: png svg pdf html (overrides config)")
    parser.add_argument("--workers", type=int, metavar="N",
                        help="Parallel workers (overrides config)")

    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plot generation, export CSV/Parquet only")
    parser.add_argument("--no-koina", action="store_true",
                        help="Disable Koina AI fallback predictions")

    parser.add_argument("--koina-model", metavar="ID",
                        help="Koina model ID (overrides config)")
    parser.add_argument("--koina-ce", type=int, metavar="EV",
                        help="Koina collision energy (overrides config)")
    parser.add_argument("--koina-instrument", metavar="TYPE",
                        help="Koina instrument type (overrides config)")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output (DEBUG logging)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress non-error output")

    return parser
