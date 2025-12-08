#!/usr/bin/env python3
"""
mzml_filter_cli.py - Command-line interface for batch mzML filtering operations.
Provides argparse-based CLI for filtering mzML files using peptide lists from files or direct input.
Supports batch processing of multiple idXML/mzML pairs with flexible peptide specification options.
"""
import sys; from pathlib import Path as P
sys.path.insert(0, str(P(__file__).parent.parent / "utils"))
sys.path.insert(0, str(P(__file__).parent.parent / "core"))
from shared_imports import logging, Path, argparse
from mzml_processor import process_idxml_batch


def load_peptides_from_file(file_path: str) -> list:
    """Load peptides from a text file (one per line)."""
    peptides = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if line and not line.startswith('#'):
                peptides.append(line)
    return peptides


def main():
    parser = argparse.ArgumentParser(
        description="Filter mzML files based on peptide matches in idXML files",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument(
        '-i', '--idxml',
        required=True,
        help='Directory containing idXML files'
    )
    
    parser.add_argument(
        '-m', '--mzml',
        required=True,
        help='Directory containing original mzML files'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for filtered mzML files'
    )
    
    parser.add_argument(
        '-p', '--peptides',
        required=True,
        help='Peptides to filter (file path or comma-separated list)'
    )
    
    parser.add_argument(
        '--no-append',
        action='store_true',
        help='Overwrite existing filtered files instead of appending'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging (very verbose)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG
    elif args.verbose:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    logger = logging.getLogger(__name__)
    
    # Parse peptides
    peptides_arg = args.peptides
    
    if Path(peptides_arg).is_file():
        # Load from file
        logger.info(f"Loading peptides from file: {peptides_arg}")
        target_peptides = load_peptides_from_file(peptides_arg)
        logger.info(f"Loaded {len(target_peptides)} peptides")
    else:
        # Parse as comma-separated
        logger.info(f"Parsing peptides from command line")
        target_peptides = [p.strip() for p in peptides_arg.split(',') if p.strip()]
        logger.info(f"Found {len(target_peptides)} peptides")
    
    if not target_peptides:
        logger.error("No peptides specified!")
        sys.exit(1)
    
    # Get idXML files
    idxml_dir = Path(args.idxml)
    if not idxml_dir.is_dir():
        logger.error(f"idXML directory not found: {idxml_dir}")
        sys.exit(1)
    
    idxml_files = list(idxml_dir.glob("*.idXML"))
    if not idxml_files:
        logger.error(f"No idXML files found in {idxml_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(idxml_files)} idXML files")
    
    # Check mzML directory
    mzml_dir = Path(args.mzml)
    if not mzml_dir.is_dir():
        logger.error(f"mzML directory not found: {mzml_dir}")
        sys.exit(1)
    
    mzml_files = list(mzml_dir.glob("*.mzML"))
    logger.info(f"Found {len(mzml_files)} mzML files")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Process
    append_mode = not args.no_append
    logger.info(f"Append mode: {append_mode}")
    
    print("\n" + "="*70)
    print("mzML Filter - Batch Processing")
    print("="*70)
    print(f"idXML files:      {len(idxml_files)}")
    print(f"mzML files:       {len(mzml_files)}")
    print(f"Target peptides:  {len(target_peptides)}")
    print(f"Output directory: {output_dir}")
    print(f"Append mode:      {append_mode}")
    print("="*70 + "\n")
    
    try:
        filtered_files = process_idxml_batch(
            idxml_paths=[str(f) for f in idxml_files],
            target_peptides=target_peptides,
            mzml_dir=str(mzml_dir),
            output_dir=str(output_dir),
            append_mode=append_mode
        )
        
        print("\n" + "="*70)
        print("Processing Complete!")
        print("="*70)
        print(f"Successfully filtered {len(filtered_files)} mzML files")
        print(f"\nFiltered files:")
        for original, filtered in filtered_files.items():
            filtered_size = Path(filtered).stat().st_size / (1024 * 1024)
            print(f"  {original:40s} -> {Path(filtered).name:40s} ({filtered_size:.2f} MB)")
        print("="*70)
        
        return 0
        
    except Exception as e:
        logger.error(f"Error during processing: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())