"""
Two-layer CLI config merger: YAML < argparse.

Defaults come from :data:`mspecviz.common.defaults.DEFAULTS` (baked into
:class:`CLIConfig`'s field defaults). YAML config and CLI flags override,
in that order. ``build_config`` runs the layers; ``build_parser`` is the
companion argparse surface.
"""
import argparse

from mspecviz.common.text_utils import get_nested

from cli.config import CLIConfig, expand_globs


def _as_tuple(val):
    """Cast list → tuple; pass tuples through unchanged."""
    return tuple(val) if isinstance(val, list) else val


# --- YAML → CLIConfig mapping table ---------------------------------------
# Each entry: (yaml_section, yaml_key, cfg_attr, optional_cast).
# YAML uses snake_case throughout.
_YAML_OVERRIDES = [
    ("annotation", "ion_types", "ion_types", None),
    ("annotation", "mz_tolerance", "mz_tolerance", float),
    ("annotation", "add_neutral_losses", "add_neutral_losses", bool),
    ("annotation", "add_first_prefix_ion", "add_first_prefix_ion", bool),
    ("annotation", "remove_precursor_peaks", "remove_precursor_peaks", bool),
    ("annotation", "intensity_mode", "intensity_mode", None),
    ("export", "directory", "export_directory", None),
    ("export", "project_name", "project_name", None),
    ("export", "formats", "export_formats", None),
    ("export", "export_dpi", "export_dpi", int),
    ("export", "figure_size", "figure_size", _as_tuple),
    ("export", "plot_unpaired_queries", "plot_unpaired_queries", bool),
    ("summary", "metric", "summary_metric", None),
    ("summary", "plot_types", "summary_plot_types", None),
    ("performance", "workers", "workers", int),
    ("koina", "enabled", "koina_enabled", bool),
    ("koina", "model", "koina_model", None),
    ("koina", "collision_energy", "koina_collision_energy", int),
    ("koina", "instrument_type", "koina_instrument_type", None),
]

# --- argparse → CLIConfig mapping table -----------------------------------
# Each entry: (argparse_dest, cfg_attr, optional_transform).
# Transforms run only when the attr has a meaningful value (not None / falsy).
_CLI_OVERRIDES = [
    ("idxml", "idxml_files", expand_globs),
    ("mzml", "mzml_files", expand_globs),
    ("reference_substring", "reference_substring", None),
    ("ion_types", "ion_types", None),
    ("mz_tolerance", "mz_tolerance", None),
    ("export_directory", "export_directory", None),
    ("project_name", "project_name", None),
    ("formats", "export_formats", None),
    ("workers", "workers", None),
    ("koina_model", "koina_model", None),
    ("koina_collision_energy", "koina_collision_energy", None),
    ("koina_instrument_type", "koina_instrument_type", None),
    ("summary_metric", "summary_metric", None),
    ("summary_plot_types", "summary_plot_types", None),
]


def _apply_yaml_overrides(cfg: CLIConfig, yaml_cfg: dict) -> None:
    """Apply YAML config overrides."""
    cfg.idxml_files = expand_globs(get_nested(yaml_cfg, "input", "idxml_files", default=[]))
    cfg.mzml_files = expand_globs(get_nested(yaml_cfg, "input", "mzml_files", default=[]))
    cfg.reference_substring = get_nested(yaml_cfg, "reference", "substring", default="")

    for section_name, yaml_key, attr, cast in _YAML_OVERRIDES:
        # A YAML key with no value parses to None, not {} — coerce to empty dict.
        section = yaml_cfg.get(section_name) or {}
        if yaml_key in section:
            val = cast(section[yaml_key]) if cast else section[yaml_key]
            setattr(cfg, attr, val)


def _apply_cli_overrides(cfg: CLIConfig, args: argparse.Namespace) -> None:
    """Apply argparse CLI overrides (highest priority)."""
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
    """Merge defaults, YAML config, and CLI overrides into a CLIConfig."""
    cfg = CLIConfig()
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

    parser.add_argument("--export-directory", metavar="DIR",
                        help="Export directory (overrides config)")
    parser.add_argument("--project-name", metavar="NAME",
                        help="Project name (overrides config)")
    parser.add_argument("--formats", nargs="+", metavar="FMT",
                        choices=["png", "svg", "pdf", "html"],
                        help="Export formats: png svg pdf html (overrides config)")
    parser.add_argument("--workers", type=int, metavar="N",
                        help="Parallel workers (overrides config)")

    parser.add_argument("--summary-metric", metavar="NAME",
                        help="Default summary plot metric (overrides config)")
    parser.add_argument("--summary-plot-types", nargs="+", metavar="TYPE",
                        choices=["dual_axis", "ridgeline", "heatmap"],
                        help="Summary plot types to render (overrides config)")

    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plot generation, export CSV/Parquet only")
    parser.add_argument("--no-koina", action="store_true",
                        help="Disable Koina AI fallback predictions")

    parser.add_argument("--koina-model", metavar="ID",
                        help="Koina model ID (overrides config)")
    parser.add_argument("--koina-collision-energy", type=int, metavar="EV",
                        help="Koina collision energy in eV (overrides config)")
    parser.add_argument("--koina-instrument-type", metavar="TYPE",
                        help="Koina instrument type (overrides config)")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output (DEBUG logging)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress non-error output")

    return parser
