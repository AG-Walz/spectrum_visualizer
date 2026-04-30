"""ProcessPoolExecutor-picklable idXML scanning workers (dispatched by ``search.orchestration``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from psm_utils.io.idxml import IdXMLReader

from mspecviz.search.idxml_parsing import mod_key_from_peptidoform
from mspecviz.search.models import PeptideMatch

logger = logging.getLogger(__name__)


def _build_match_from_psm(
    psm,
    idxml_path: str,
    idxml_name: str,
    mzml_by_name: Dict[str, str],
) -> Optional[Tuple[str, PeptideMatch]]:
    """Build a (variant_key, PeptideMatch) from a psm_utils PSM, or None.

    ``mzml_by_name`` maps lowercase basename → full path. The unique-basename
    invariant (enforced at input) makes this lookup unambiguous.
    """
    base_seq = psm.peptidoform.sequence.upper()

    mzml_filename = (psm.run + ".mzML") if psm.run else ""
    resolved_path = mzml_by_name.get(mzml_filename.lower(), "") if mzml_filename else ""

    mod_key = mod_key_from_peptidoform(psm.peptidoform)
    variant_key = f"{base_seq}|{mod_key}" if mod_key else base_seq

    match = PeptideMatch(
        peptidoform=psm.peptidoform.proforma.split("/")[0],
        base_sequence=base_seq,
        spectrum_id=psm.spectrum_id or "",
        charge=psm.get_precursor_charge() or 0,
        score=psm.score or 0.0,
        idxml_file=idxml_name,
        mzml_filename=mzml_filename,
        id_merge_index=-1,
        mzml_resolved_path=resolved_path,
        idxml_resolved_path=idxml_path,
    )
    return variant_key, match


def _search_single_idxml(
    idxml_path: str,
    idxml_name: str,
    target_set: Set[str],
    mzml_by_name: Dict[str, str],
) -> Tuple[str, list, list]:
    """Worker: scan one idXML for target peptides via psm_utils PSM iteration.

    Returns ``(idxml_name, [(variant_key, PeptideMatch)], errors)``. Module-level
    so ProcessPoolExecutor can pickle it; ``psm.run`` already carries the
    resolved mzML basename and ``mzml_by_name`` (lowercase basename → full path)
    binds it to a concrete file.
    """
    matches: list = []
    errors: list = []

    try:
        try:
            reader = IdXMLReader(str(idxml_path))
        except Exception as load_error:
            errors.append(f"No peptide identifications in {idxml_name}: {load_error}")
            return idxml_name, matches, errors

        # Manual iteration: psm_utils' __iter__ is a generator, so a single bad PSM
        # (e.g. null rescoring feature in unfiltered _pout.idXML) kills the stream.
        psm_iter = iter(reader)
        psms_read = 0
        parse_error: Optional[str] = None
        while True:
            try:
                psm = next(psm_iter)
                psms_read += 1
            except StopIteration:
                break
            except Exception as e:
                parse_error = str(e)
                logger.debug(f"PSM parse error in {idxml_name} after {psms_read} PSMs: {e}")
                break  # generator is dead — no point continuing

            try:
                if psm.peptidoform.sequence.upper() not in target_set:
                    continue
                built = _build_match_from_psm(psm, idxml_path, idxml_name, mzml_by_name)
                if built is not None:
                    matches.append(built)
            except Exception as e:
                logger.debug(f"Error processing PSM: {e}")
                continue

        if parse_error is not None:
            errors.append(
                f"{idxml_name}: stopped after {psms_read} PSMs — a PSM could not be "
                f"parsed (null rescoring feature value). This happens with unfiltered "
                f"Percolator output (_pout.idXML) that contains PSMs missing metadata. "
                f"Switch to _pout_filtered.idXML for complete, reliable results."
            )

    except Exception as e:
        errors.append(f"Error loading {Path(idxml_path).name}: {e}")

    return idxml_name, matches, errors


def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Count PSMs per base sequence in one idXML; ``{}`` on parse error."""
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            seq = psm.peptidoform.sequence.upper()
            if seq:
                counts[seq] += 1
    except Exception as e:
        logger.warning(f"count_psms_per_peptide failed for {Path(idxml_path).name}: {e}")
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
    except Exception as e:
        logger.warning(f"_scan_single_idxml_peptides failed for {Path(idxml_path).name}: {e}")

    return idxml_path, peptides
