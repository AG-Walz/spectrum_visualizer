"""idXML scanning workers, which open and read the peptide identifications out of them, 
designed to be dispatched across multiple processes by the run_tasks utility in helpers.py
(used by core.search.peptide_search).
"""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Set, Tuple
from psm_utils.io.idxml import IdXMLReader

logger = logging.getLogger(__name__)

def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Simple reader for the Setup page (before peptides are chosen) - Count PSMs per base sequence in one idXML. 
    Returns a dict mapping unmodified peptide sequence to count, e.g. {"PEPTIDE": 42}.
    """
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            # Only count if there's a sequence (skip empty PSMs) and ignore modifications
            if psm.peptidoform.sequence:
                counts[psm.peptidoform.sequence.upper()] += 1
    # On parse errors, return whatever it counted before the error and log a warning 
    except Exception as error:
        logger.warning(f"count_psms_per_peptide failed for {Path(idxml_path).name}: {error}")
    return counts


def search_single_idxml(
    idxml_path: str, idxml_name: str,
    target_set: Set[str], mzml_by_stem: Dict[str, str],
) -> Tuple[str, list, list]:
    """Scan one idXML file for target peptides; return the matching PSMs
     with extracted metadata attached (name, matches, errors).
    """
    # Open psm_utils' idxml-reader, which yields PSMs one at a time
    matches: list = []
    try:
        reader = IdXMLReader(str(idxml_path))
    # Catch file-level errors, returning empty results with an error message instead of crashing the whole search
    except Exception as load_error:
        return idxml_name, matches, [f"No peptide identifications in {idxml_name}: {load_error}"]

    # Fill a list with all PSMs, to surface parse errors as partial results
    psms: list = []
    parse_error: str | None = None
    try:
        for psm in reader:
            psms.append(psm)
    except Exception as error:
        parse_error = str(error)

    # Filter & annotate: keep target hits, stamp resolution metadata on each PSM
    for psm in psms:
        # Skip if no sequence or not in target set 
        if psm.peptidoform.sequence not in target_set:
            continue
        # Annotate with idXML file 
        psm.metadata["idxml_file"] = idxml_name
        psm.metadata["idxml_resolved_path"] = str(idxml_path)
        # Annotate with mzML path if possible (match by run name, case-insensitive)
        if psm.run:
            resolved = mzml_by_stem.get(psm.run.lower(), "")
            if resolved:
                psm.metadata["mzml_resolved_path"] = resolved
        # Append modified sequence and metadata-enriched PSM to results (later becomes the variant_key)
        matches.append((psm.peptidoform.modified_sequence, psm))

    # Surface partial-read errors so the orchestrator can warn the user
    errors = []
    if parse_error is not None:
        errors.append(
            f"{idxml_name}: stopped after {len(psms)} PSMs — parse error: {parse_error}. "
            f"Try _pout_filtered.idXML for complete results."
        )
    return idxml_name, matches, errors