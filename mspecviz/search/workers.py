"""ProcessPoolExecutor-picklable idXML scanning workers (dispatched by ``search.orchestration``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from psm_utils.io.idxml import IdXMLReader

logger = logging.getLogger(__name__)


def _enrich_with_metadata(
    psm, idxml_path: str, idxml_name: str, mzml_by_stem: Dict[str, str],
) -> None:
    """Stuff mspecviz-specific resolution metadata into ``psm.metadata`` in place."""
    psm.metadata["idxml_file"] = idxml_name
    psm.metadata["idxml_resolved_path"] = str(idxml_path)
    if psm.run:
        resolved = mzml_by_stem.get(psm.run.lower(), "")
        if resolved:
            psm.metadata["mzml_resolved_path"] = resolved


def _process_psm_against_targets(
    psm, idxml_path: str, idxml_name: str,
    target_set: Set[str], mzml_by_stem: Dict[str, str], matches: list,
) -> None:
    """Append ``(variant_key, psm)`` to *matches* if the PSM targets a peptide."""
    try:
        if psm.peptidoform.sequence.upper() not in target_set:
            return
        _enrich_with_metadata(psm, idxml_path, idxml_name, mzml_by_stem)
        variant_key = psm.peptidoform.modified_sequence
        matches.append((variant_key, psm))
    except Exception as error:
        logger.debug(f"Error processing PSM: {error}")


def _iterate_psms_safely(
    reader, idxml_name: str,
) -> Tuple[list, int, Optional[str]]:
    """Manually drain the psm_utils generator, isolating one bad PSM from the rest.

    Returns ``(psms, psms_read, parse_error)``. We iterate manually because
    psm_utils' ``__iter__`` is a generator: one bad PSM (e.g. null rescoring
    feature in unfiltered ``_pout.idXML``) kills the stream entirely.
    """
    psms: list = []
    psm_iter = iter(reader)
    psms_read = 0
    while True:
        try:
            psm = next(psm_iter)
            psms_read += 1
            psms.append(psm)
        except StopIteration:
            return psms, psms_read, None
        except Exception as error:
            logger.debug(
                f"PSM parse error in {idxml_name} after {psms_read} PSMs: {error}"
            )
            return psms, psms_read, str(error)


def _search_single_idxml(
    idxml_path: str,
    idxml_name: str,
    target_set: Set[str],
    mzml_by_stem: Dict[str, str],
) -> Tuple[str, list, list]:
    """Worker: scan one idXML for target peptides via psm_utils PSM iteration."""
    matches: list = []
    errors: list = []

    try:
        try:
            reader = IdXMLReader(str(idxml_path))
        except Exception as load_error:
            errors.append(f"No peptide identifications in {idxml_name}: {load_error}")
            return idxml_name, matches, errors

        psms, psms_read, parse_error = _iterate_psms_safely(reader, idxml_name)
        for psm in psms:
            _process_psm_against_targets(
                psm, idxml_path, idxml_name, target_set, mzml_by_stem, matches,
            )

        if parse_error is not None:
            errors.append(
                f"{idxml_name}: stopped after {psms_read} PSMs — a PSM could not be "
                f"parsed (null rescoring feature value). This happens with unfiltered "
                f"Percolator output (_pout.idXML) that contains PSMs missing metadata. "
                f"Switch to _pout_filtered.idXML for complete, reliable results."
            )
    except Exception as error:
        errors.append(f"Error loading {Path(idxml_path).name}: {error}")

    return idxml_name, matches, errors


def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Count PSMs per base sequence in one idXML; ``{}`` on parse error."""
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            seq = psm.peptidoform.sequence.upper()
            if seq:
                counts[seq] += 1
    except Exception as error:
        logger.warning(f"count_psms_per_peptide failed for {Path(idxml_path).name}: {error}")
    return dict(counts)


def _scan_single_idxml_peptides(idxml_path: str) -> Tuple[str, Set[str]]:
    """Worker: return ``(idxml_path, {base_sequence, ...})`` for one idXML."""
    peptides: Set[str] = set()
    try:
        reader = IdXMLReader(str(idxml_path))
        for psm in reader:
            base_sequence = psm.peptidoform.sequence.upper()
            if base_sequence:
                peptides.add(base_sequence)
    except Exception as error:
        logger.warning(f"_scan_single_idxml_peptides failed for {Path(idxml_path).name}: {error}")

    return idxml_path, peptides
