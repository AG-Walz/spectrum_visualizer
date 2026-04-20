"""Config merger: :data:`mspecviz.common.defaults.DEFAULTS` < YAML < argparse."""

import argparse
from dataclasses import fields
from cli.config import CLIConfig, expand_globs


def _copy_into_config(config: CLIConfig, source: dict) -> None:
    """Copy every ``source`` entry whose key matches a CLIConfig field.

    ``None`` and empty-list values mean "layer did not supply this key" and are
    skipped (left as default) so they can't overwrite a lower-precedence value.
    """
    config_fields = {f.name for f in fields(config)}
    for field_name, field_value in source.items():
        if field_name not in config_fields:
            continue
        if field_value is None or field_value == []:
            continue
        # ``figure_size`` is typed as Tuple but YAML/argparse deliver a list.
        if field_name == "figure_size" and isinstance(field_value, list):
            field_value = tuple(field_value)
        if field_name in ("idxml_files", "mzml_files"):
            field_value = expand_globs(field_value)
        setattr(config, field_name, field_value)


def build_config(yaml_config: dict, args: argparse.Namespace) -> CLIConfig:
    """Merge defaults, YAML config, and CLI overrides into a CLIConfig."""
    config = CLIConfig()
    _copy_into_config(config, yaml_config)
    _copy_into_config(config, vars(args))

    # --no-koina is a negated toggle with no matching CLIConfig field; it flips
    # ``koina_enabled`` only when the user passed the flag (None = absent).
    if args.no_koina:
        config.koina_enabled = False

    return config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mspecviz-cli",
        description="MSpecViz CLI -- Terminal-based MS/MS spectrum analysis pipeline.",
        epilog="Example: python -m cli --config cli/config/config_example.yaml --verbose",
    )

    parser.add_argument("-c", "--config", required=True,
                        help="Path to YAML config file")

    # ``dest=`` is set whenever argparse's default name doesn't match the
    # target CLIConfig field destination.
    parser.add_argument("--idxml", dest="idxml_files", nargs="+", metavar="PATTERN",
                        help="idXML file paths or glob patterns (overrides config)")
    parser.add_argument("--mzml", dest="mzml_files", nargs="+", metavar="PATTERN",
                        help="mzML file paths or glob patterns (overrides config)")
    parser.add_argument("--peptides", dest="peptide_source", metavar="PATH_OR_ALL",
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
    parser.add_argument("--formats", dest="export_formats", nargs="+", metavar="FMT",
                        choices=["png", "svg", "pdf", "html"],
                        help="Export formats: png svg pdf html (overrides config)")
    parser.add_argument("--workers", type=int, metavar="N",
                        help="Parallel workers (overrides config)")

    parser.add_argument("--summary-metric", metavar="NAME",
                        help="Default summary plot metric (overrides config)")
    parser.add_argument("--summary-plot-types", nargs="+", metavar="TYPE",
                        choices=["dual_axis", "ridgeline", "heatmap"],
                        help="Summary plot types to render (overrides config)")

    # ``default=None`` turns argparse's usual False into "not supplied" so the
    # copy loop's None-filter lets a YAML-set ``no_plots: true`` survive when
    # the flag is absent. Present → True → overrides; absent → None → skipped.
    parser.add_argument("--no-plots", action="store_true", default=None,
                        help="Skip plot generation, export CSV/Parquet only")
    parser.add_argument("--no-koina", action="store_true", default=None,
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
