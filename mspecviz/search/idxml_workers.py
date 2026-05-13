"""ProcessPoolExecutor-picklable idXML scanning workers (dispatched by ``search.peptide_search``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Set, Tuple

from psm_utils.io.idxml import IdXMLReader

logger = logging.getLogger(__name__)


def search_single_idxml(
    idxml_path: str, idxml_name: str,
    target_set: Set[str], mzml_by_stem: Dict[str, str],
) -> Tuple[str, list, list]:
    """Worker: scan one idXML for target peptides; return ``(name, matches, errors)``."""
    # Open reader: bail with a single error if the file can't even be parsed
    matches: list = []
    try:
        reader = IdXMLReader(str(idxml_path))
    except Exception as load_error:
        return idxml_name, matches, [f"No peptide identifications in {idxml_name}: {load_error}"]

    # Drain reader: isolate erroneous PSMs from the rest
    psms: list = []
    parse_error: str | None = None
    try:
        for psm in reader:
            psms.append(psm)
    except Exception as error:
        logger.debug(f"PSM parse error in {idxml_name} after {len(psms)} PSMs: {error}")
        parse_error = str(error)

    # Filter & annotate: keep target hits, stamp resolution metadata on each PSM
    for psm in psms:
        if psm.peptidoform.sequence.upper() not in target_set:
            continue
        psm.metadata["idxml_file"] = idxml_name
        psm.metadata["idxml_resolved_path"] = str(idxml_path)
        if psm.run:
            resolved = mzml_by_stem.get(psm.run.lower(), "")
            if resolved:
                psm.metadata["mzml_resolved_path"] = resolved
        matches.append((psm.peptidoform.modified_sequence, psm))

    # Surface partial-read errors so the orchestrator can warn the user
    errors = []
    if parse_error is not None:
        errors.append(
            f"{idxml_name}: stopped after {len(psms)} PSMs — parse error: {parse_error}. "
            f"Try _pout_filtered.idXML for complete results."
        )
    return idxml_name, matches, errors


def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Count PSMs per base sequence in one idXML; partial counts on parse error."""
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            if psm.peptidoform.sequence:
                counts[psm.peptidoform.sequence.upper()] += 1
    except Exception as error:
        logger.warning(f"count_psms_per_peptide failed for {Path(idxml_path).name}: {error}")
    return counts
