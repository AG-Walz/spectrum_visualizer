#!/usr/bin/env python3
"""MSpecViz CLI — ``python -m cli --config config.yaml [overrides]``.

Reuses the same backend as the Streamlit app; argparse/config/output/plots/
pipeline live in sibling ``cli.*`` modules. Exit codes: 0 ok, 1 pipeline
failure, 2 config/input error.
"""
import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cli.config import load_peptides, load_yaml_config
from cli.config_builder import build_config, build_parser
from cli.pipeline import run_pipeline


def main():
    args = build_parser().parse_args()

    log_level = logging.ERROR if args.quiet else (
        logging.DEBUG if args.verbose else logging.WARNING)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    yaml_config = load_yaml_config(args.config)
    config = build_config(yaml_config, args)

    if not config.peptide_source:
        print("Error: No peptides specified. Use 'peptide_source' in config or --peptides.", file=sys.stderr)
        sys.exit(2)
    config.peptides, config.peptides_mode = load_peptides(
        config.peptide_source, config.idxml_files, config.workers, args.quiet,
    )
    if not config.peptides:
        print("Error: Peptide list is empty.", file=sys.stderr)
        sys.exit(2)

    sys.exit(run_pipeline(config, quiet=args.quiet))


if __name__ == "__main__":
    main()
