"""ProcessPoolExecutor-picklable idXML scanning workers (dispatched by ``search.peptide_search``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

from psm_utils.io.idxml import IdXMLReader

if TYPE_CHECKING:
    from psm_utils import PSM

logger = logging.getLogger(__name__)


def _iterate_psms_safely(
    reader, idxml_name: str,
) -> Tuple[List["PSM"], int, Optional[str]]:
    """Drain a psm_utils IdXMLReader, isolating one bad PSM from the rest.

    psm_utils' ``__iter__`` is a generator that dies on a single malformed PSM
    (e.g. null rescoring feature in unfiltered ``_pout.idXML``). We catch the
    generator failure so earlier PSMs in the same file are still usable.
    Returns ``(psms, psms_read, parse_error)``.
    """
    psms, psms_read = [], 0
    psm_iter = iter(reader)
    while True:
        try:
            psms.append(next(psm_iter))
            psms_read += 1
        except StopIteration:
            return psms, psms_read, None
        except Exception as error:
            logger.debug(f"PSM parse error in {idxml_name} after {psms_read} PSMs: {error}")
            return psms, psms_read, str(error)


def _search_single_idxml(
    idxml_path: str, idxml_name: str,
    target_set: Set[str], mzml_by_stem: Dict[str, str],
) -> Tuple[str, list, list]:
    """Worker: scan one idXML for target peptides; return ``(name, matches, errors)``."""
    matches: list = []
    try:
        reader = IdXMLReader(str(idxml_path))
    except Exception as load_error:
        return idxml_name, matches, [f"No peptide identifications in {idxml_name}: {load_error}"]

    psms, psms_read, parse_error = _iterate_psms_safely(reader, idxml_name)
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

    errors = []
    if parse_error is not None:
        errors.append(
            f"{idxml_name}: stopped after {psms_read} PSMs — parse error: {parse_error}. "
            f"Try _pout_filtered.idXML for complete results."
        )
    return idxml_name, matches, errors


def _iter_base_sequences(idxml_path: str, caller: str):
    """Yield uppercase base sequences from every PSM in *idxml_path*; warn on parse error."""
    try:
        for psm in IdXMLReader(str(idxml_path)):
            if psm.peptidoform.sequence:
                yield psm.peptidoform.sequence.upper()
    except Exception as error:
        logger.warning(f"{caller} failed for {Path(idxml_path).name}: {error}")


def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Count PSMs per base sequence in one idXML; ``{}`` on parse error."""
    return Counter(_iter_base_sequences(idxml_path, "count_psms_per_peptide"))


def _scan_single_idxml_peptides(idxml_path: str) -> Tuple[str, Set[str]]:
    """Worker: return ``(idxml_path, {base_sequence, ...})`` for one idXML."""
    return idxml_path, set(_iter_base_sequences(idxml_path, "_scan_single_idxml_peptides"))
