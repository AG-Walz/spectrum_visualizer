"""
Three-layer CLI config merger: ``settings.json`` < YAML < argparse.

Each layer's mapping is a small table so adding a new parameter is a
one-line addition in the relevant table, not a fourth scattered branch.
``build_config`` runs the three layers in order; ``build_parser`` is the
companion argparse surface.
"""
import argparse

from mspecviz.common.text_utils import get_nested

from cli.config import CLIConfig, expand_globs, load_settings_defaults


def _as_tuple(val):
    """Cast list → tuple; pass tuples through unchanged."""
    return tuple(val) if isinstance(val, list) else val


# --- settings.json → CLIConfig mapping table (layer 1 defaults) ------------
# Each entry: (cfg_attr, (settings_key_path), default_value, optional_cast)
_SETTINGS_DEFAULTS = [
    ("ion_types", ("ion_types", "default"), ["b", "y"], None),
    ("mz_tolerance", ("mz_tolerance", "default"), 0.02, None),
    ("add_losses", ("annotation", "add_losses"), False, None),
    ("add_first_prefix_ion", ("annotation", "add_first_prefix_ion"), True, None),
    ("remove_precursor_peaks", ("annotation", "remove_precursor_peaks"), False, None),
    ("export_dir", ("export", "default_directory"), "./exports", None),
    ("project_name", ("export", "project_name"), "MSpecViz_analysis", None),
    ("export_formats", ("export", "formats"), ["png"], None),
    ("dpi", ("export", "dpi"), 300, None),
    ("figsize", ("export", "figsize"), (12, 8), _as_tuple),
    ("workers", ("performance", "workers"), 4, None),
    ("koina_enabled", ("koina", "enabled"), True, None),
    ("koina_model", ("koina", "model"), "ms2pip_timsTOF2024", None),
    ("koina_ce", ("koina", "collision_energy"), 35, None),
    ("koina_instrument", ("koina", "instrument_type"), "TIMSTOF", None),
]


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
    """Populate cfg from settings.json (layer 1).

    Driven by :data:`_SETTINGS_DEFAULTS` — add a new parameter there rather
    than extending this function.
    """
    for attr, path, default, cast in _SETTINGS_DEFAULTS:
        val = get_nested(defaults, *path, default=default)
        setattr(cfg, attr, cast(val) if cast else val)


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
