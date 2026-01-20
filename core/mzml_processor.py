"""mzML Filter - Extract spectra for target peptides from large files.

Filters mzML/idXML files to contain only spectra matching target peptides,
enabling memory-efficient processing of large datasets (reduces size by 90%+).
Key functions: extract_peptide_relevant_spectra(), filter_mzml(), filter_idxml()
"""
import sys; from pathlib import Path as P; sys.path.insert(0, str(P(__file__).parent.parent / "utils"))
from shared_imports import logging, oms, Path, Dict, List, Tuple, Set, defaultdict

logger = logging.getLogger(__name__)


def extract_peptide_relevant_spectra(
    idxml_path: str,
    target_peptides: List[str],
    mzml_dir: str = None
) -> Dict[str, List[str]]:
    """
    Extract spectrum IDs for target peptides (including modified versions).
    
    Args:
        idxml_path: Path to idXML file
        target_peptides: List of unmodified peptide sequences to search for
        mzml_dir: Optional directory containing mzML files (for validation)
    
    Returns:
        Dictionary mapping mzML filename -> list of spectrum IDs
        Format: {
            "sample1.mzML": ["index=123", "index=456"],
            "sample2.mzML": ["index=789"]
        }
    """
    logger.info(f"Processing idXML: {Path(idxml_path).name}")
    logger.info(f"Searching for {len(target_peptides)} target peptides")
    
    # Convert target peptides to set for faster lookup (with quotes for exact matching)
    target_set = set(f'"{seq}"' for seq in target_peptides)
    logger.debug(f"Target peptides (quoted): {list(target_set)[:5]}...")
    
    # Load idXML file
    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)
    
    logger.info(f"Loaded {len(peptide_ids)} peptide identifications from idXML")
    
    # Create a map of identifier -> mzML filename from protein_ids for fast lookup
    run_map = {}
    for protein_id in protein_ids:
        mzml_file = None
        if protein_id.metaValueExists("spectra_data"):
            # Handle bytes/list logic safely
            val = protein_id.getMetaValue("spectra_data")
            if isinstance(val, list) and val: val = val[0]
            if isinstance(val, bytes): val = val.decode()
            mzml_file = Path(val).name
        
        if mzml_file:
            run_map[protein_id.getIdentifier()] = mzml_file

    # Dictionary to store results: mzML -> [spectrum_ids]
    mzml_to_spectra = defaultdict(list)
    
    matches_found = 0
    unique_mzmls = set()
    
    # Process each peptide identification
    for peptide_id in peptide_ids:
        # Get the spectrum reference (native ID)
        spectrum_id = peptide_id.getMetaValue("spectrum_reference") if peptide_id.metaValueExists("spectrum_reference") else ""
        
        # Extract mzML filename from peptide_id metadata
        # This is usually stored in the run or file description
        mzml_file = None
        if peptide_id.metaValueExists("file_origin"):
            mzml_file = Path(peptide_id.getMetaValue("file_origin")).name
        
        # If no file_origin, we need to extract from the PeptideIdentification itself
        # OpenMS stores this in the identifier field or we can infer from protein_ids
        if not mzml_file and protein_ids:
            identifier = peptide_id.getIdentifier()
            
            # Try exact match
            if identifier in run_map:
                mzml_file = run_map[identifier]
            else:
                # Try finding a run that matches the prefix of the peptide identifier
                for run_id, run_file in run_map.items():
                    if identifier.startswith(run_id):
                        mzml_file = run_file
                        break
        
        # Fallback: extract from spectrum_reference if it contains path
        if not mzml_file and spectrum_id:
            # Sometimes spectrum_reference contains file path
            if "/" in spectrum_id or "\\" in spectrum_id:
                mzml_file = Path(spectrum_id.split()[0]).name
        
        if not mzml_file:
            logger.warning(f"Could not determine mzML file for spectrum {spectrum_id}")
            continue
        
        # Ensure .mzML extension
        if not mzml_file.endswith('.mzML'):
            mzml_file = f"{mzml_file}.mzML"
        
        # CLEANUP: Clean spectrum_id (strip filename prefix) to ensure matching with mzML Native IDs
        if spectrum_id:
            # Ensure string
            if isinstance(spectrum_id, bytes):
                spectrum_id = spectrum_id.decode()
            
            clean_id = spectrum_id.strip()
            
            # Heuristic 1: If it contains "scan=", "index=", or "spectrum=", assume that starts the ID
            for tag in ["scan=", "index=", "spectrum=", "controllerType="]:
                if tag in clean_id:
                    idx = clean_id.find(tag)
                    # If the tag is not at the start, check if it's preceded by space or colon
                    if idx > 0 and clean_id[idx-1] in [' ', ':']:
                        clean_id = clean_id[idx:]
                        break
            
            # Heuristic 2: If we didn't find a tag but it looks like "file.mzML 123", strip the file
            if clean_id == spectrum_id.strip() and (".mzML" in clean_id or ".mzXML" in clean_id):
                 parts = clean_id.split()
                 if len(parts) > 1:
                     # Check if first part looks like a file
                     if parts[0].endswith(".mzML") or parts[0].endswith(".mzXML") or parts[0].endswith(".mzML\"") or parts[0].endswith(".mzXML\""):
                         clean_id = " ".join(parts[1:])
            
            spectrum_id = clean_id
        
        # Process each peptide hit
        hits = peptide_id.getHits()
        for hit in hits:
            # Get the peptide sequence
            sequence = hit.getSequence()
            
            # Convert to unmodified string using AASequence
            unmodified = sequence.toUnmodifiedString()
            
            # Create quoted version for exact matching
            quoted_unmodified = f'"{unmodified}"'
            
            # Check if this matches any target peptide
            if quoted_unmodified in target_set:
                mzml_to_spectra[mzml_file].append(spectrum_id)
                matches_found += 1
                unique_mzmls.add(mzml_file)
                logger.debug(f"Match found: {unmodified} in {mzml_file} (spectrum: {spectrum_id})")
                break  # Only count once per peptide_id
    
    logger.info(f"✅ Found {matches_found} matching spectra across {len(unique_mzmls)} mzML file(s)")
    logger.info(f"   mzML files: {list(unique_mzmls)[:5]}...")
    
    # Convert defaultdict to regular dict
    return dict(mzml_to_spectra)


def filter_mzml(
    original_mzml_path: str,
    spectrum_ids: List[str],
    output_dir: str,
    append_if_exists: bool = True
) -> str:
    """
    Create or append to a filtered mzML file containing only specified spectra.
    
    Args:
        original_mzml_path: Path to original mzML file
        spectrum_ids: List of spectrum IDs (native IDs) to include
        output_dir: Directory to save filtered mzML files
        append_if_exists: If True, append to existing file; if False, overwrite
    
    Returns:
        Path to the filtered mzML file
    """
    original_name = Path(original_mzml_path).name
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    output_file = output_path / original_name
    
    logger.info(f"Filtering {original_name}: {len(spectrum_ids)} spectra")
    
    # Load original experiment
    original_exp = oms.MSExperiment()
    oms.MzMLFile().load(original_mzml_path, original_exp)
    
    logger.info(f"  Original file has {original_exp.getNrSpectra()} spectra")
    
    # Create spectrum lookup for efficient searching
    spectrum_lookup = oms.SpectrumLookup()
    
    # Detect spectrum ID format
    if original_exp.getNrSpectra() > 0:
        first_spectrum = original_exp.getSpectrum(0)
        native_id = first_spectrum.getNativeID()
        
        if 'scan=' in native_id:
            pattern = "scan=(?<SCAN>\\d+)"
        elif 'index=' in native_id:
            pattern = "index=(?<SCAN>\\d+)"
        elif 'spectrum=' in native_id:
            pattern = "spectrum=(?<SCAN>\\d+)"
        else:
            pattern = "(?<SCAN>\\d+)"
        
        logger.info(f"  Using pattern: {pattern}")
        spectrum_lookup.readSpectra(original_exp, pattern)
    
    # Check if output file already exists and we're appending
    if output_file.exists() and append_if_exists:
        logger.info(f"  Appending to existing file: {output_file.name}")
        
        # Load existing filtered experiment
        filtered_exp = oms.MSExperiment()
        oms.MzMLFile().load(str(output_file), filtered_exp)
        
        # Get existing spectrum IDs to avoid duplicates
        existing_ids = set()
        for i in range(filtered_exp.getNrSpectra()):
            existing_ids.add(filtered_exp.getSpectrum(i).getNativeID())
        
        logger.info(f"  Existing file has {len(existing_ids)} spectra")
        
        # Add new spectra (excluding duplicates)
        added_count = 0
        for spec_id in spectrum_ids:
            if spec_id not in existing_ids:
                try:
                    idx = spectrum_lookup.findByNativeID(spec_id)
                    spectrum = original_exp.getSpectrum(idx)
                    
                    # ENFORCE MS2
                    if spectrum.getMSLevel() != 2:
                        # logger.debug(f"  Skipping non-MS2 spectrum {spec_id} (Level={spectrum.getMSLevel()})")
                        continue
                        
                    filtered_exp.addSpectrum(spectrum)
                    added_count += 1
                except Exception as e:
                    logger.warning(f"  Could not find spectrum {spec_id}: {e}")
        
        logger.info(f"  Added {added_count} new spectra (skipped {len(spectrum_ids) - added_count} duplicates)")
        
    else:
        # Create new filtered experiment
        logger.info(f"  Creating new filtered file: {output_file.name}")
        
        filtered_exp = oms.MSExperiment()
        
        # Copy experimental settings (metadata) from original
        # Note: setExperimentalSettings doesn't exist in newer PyOpenMS
        # Instead, we copy relevant metadata manually
        try:
            # Copy contact information
            filtered_exp.setContacts(original_exp.getContacts())
        except:
            pass
        
        try:
            # Copy instrument information
            filtered_exp.setInstrument(original_exp.getInstrument())
        except:
            pass
        
        try:
            # Copy source files
            filtered_exp.setSourceFiles(original_exp.getSourceFiles())
        except:
            pass
        
        try:
            # Copy sample information
            filtered_exp.setSample(original_exp.getSample())
        except:
            pass
        
        # Add only the specified spectra
        added_count = 0
        for spec_id in spectrum_ids:
            try:
                idx = spectrum_lookup.findByNativeID(spec_id)
                spectrum = original_exp.getSpectrum(idx)
                
                # ENFORCE MS2
                if spectrum.getMSLevel() != 2:
                    # logger.debug(f"  Skipping non-MS2 spectrum {spec_id} (Level={spectrum.getMSLevel()})")
                    continue
                    
                filtered_exp.addSpectrum(spectrum)
                added_count += 1
            except Exception as e:
                logger.warning(f"  Could not find spectrum {spec_id}: {e}")
        
        logger.info(f"  Added {added_count} spectra")
    
    # Sort spectra by retention time (important for downstream processing)
    filtered_exp.sortSpectra(True)
    
    # Save filtered experiment
    oms.MzMLFile().store(str(output_file), filtered_exp)
    
    logger.info(f"✅ Saved filtered mzML: {output_file.name} ({filtered_exp.getNrSpectra()} total spectra)")
    
    return str(output_file)


def process_idxml_batch(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_dir: str,
    output_dir: str,
    append_mode: bool = True
) -> Dict[str, str]:
    """
    Process multiple idXML files and create filtered mzML files.
    
    Args:
        idxml_paths: List of paths to idXML files
        target_peptides: List of unmodified peptide sequences to filter
        mzml_dir: Directory containing original mzML files
        output_dir: Directory to save filtered mzML files
        append_mode: If True, append to existing filtered files
    
    Returns:
        Dictionary mapping original mzML name -> filtered mzML path
    """
    logger.info("="*60)
    logger.info("Starting batch processing")
    logger.info("="*60)
    logger.info(f"idXML files: {len(idxml_paths)}")
    logger.info(f"Target peptides: {len(target_peptides)}")
    logger.info(f"mzML directory: {mzml_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Append mode: {append_mode}")
    
    # Step 1: Extract all relevant spectra from all idXML files
    all_mzml_to_spectra = defaultdict(list)
    
    for i, idxml_path in enumerate(idxml_paths, 1):
        logger.info(f"\n[{i}/{len(idxml_paths)}] Processing {Path(idxml_path).name}")
        
        try:
            mzml_spectra_map = extract_peptide_relevant_spectra(
                idxml_path=idxml_path,
                target_peptides=target_peptides,
                mzml_dir=mzml_dir
            )
            
            # Merge results (accumulate spectrum IDs per mzML)
            for mzml_file, spec_ids in mzml_spectra_map.items():
                all_mzml_to_spectra[mzml_file].extend(spec_ids)
            
        except Exception as e:
            logger.error(f"Error processing {Path(idxml_path).name}: {e}", exc_info=True)
            continue
    
    # Remove duplicates from spectrum ID lists
    for mzml_file in all_mzml_to_spectra:
        all_mzml_to_spectra[mzml_file] = list(set(all_mzml_to_spectra[mzml_file]))
    
    logger.info("\n" + "="*60)
    logger.info("Extraction complete")
    logger.info("="*60)
    logger.info(f"Total mzML files to filter: {len(all_mzml_to_spectra)}")
    for mzml_file, spec_ids in all_mzml_to_spectra.items():
        logger.info(f"  {mzml_file}: {len(spec_ids)} spectra")
    
    # Step 2: Create filtered mzML files
    logger.info("\n" + "="*60)
    logger.info("Creating filtered mzML files")
    logger.info("="*60)
    
    filtered_files = {}
    
    for i, (mzml_file, spec_ids) in enumerate(all_mzml_to_spectra.items(), 1):
        logger.info(f"\n[{i}/{len(all_mzml_to_spectra)}] Filtering {mzml_file}")
        
        # Find original mzML file
        original_path = Path(mzml_dir) / mzml_file
        
        if not original_path.exists():
            logger.warning(f"  Original mzML not found: {original_path}")
            continue
        
        try:
            filtered_path = filter_mzml(
                original_mzml_path=str(original_path),
                spectrum_ids=spec_ids,
                output_dir=output_dir,
                append_if_exists=append_mode
            )
            
            filtered_files[mzml_file] = filtered_path
            
        except Exception as e:
            logger.error(f"Error filtering {mzml_file}: {e}", exc_info=True)
            continue
    
    logger.info("\n" + "="*60)
    logger.info("Batch processing complete")
    logger.info("="*60)
    logger.info(f"Successfully created {len(filtered_files)} filtered mzML files")
    logger.info(f"Output directory: {output_dir}")
    
    return filtered_files


if __name__ == "__main__":
    # Example usage
    import sys
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Example: Process idXML files to create filtered mzMLs
    example_peptides = [
        "PEPTIDEK",
        "SEQUENCEM",
        "PROTEINR"
    ]
    
    example_idxml_dir = "/path/to/idxml/files"
    example_mzml_dir = "/path/to/mzml/files"
    example_output_dir = "/path/to/filtered_mzml"
    
    # Get all idXML files
    idxml_files = list(Path(example_idxml_dir).glob("*.idXML"))
    
    if idxml_files:
        filtered_files = process_idxml_batch(
            idxml_paths=[str(f) for f in idxml_files],
            target_peptides=example_peptides,
            mzml_dir=example_mzml_dir,
            output_dir=example_output_dir,
            append_mode=True
        )
        
        print("\nFiltered mzML files created:")
        for original, filtered in filtered_files.items():
            print(f"  {original} -> {filtered}")
    else:
        print(f"No idXML files found in {example_idxml_dir}")
