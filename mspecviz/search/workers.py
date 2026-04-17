"""Module-level workers for ProcessPoolExecutor-based idXML scanning.

Each worker takes a single idXML path (plus context) and returns a result that
is serialisable and process-safe. They live as top-level functions so they can
be pickled when shipped to a worker process. The orchestration layer
(``mspecviz.search.orchestration``) submits these via ``run_tasks``.
"""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from psm_utils.io.idxml import IdXMLReader

from mspecviz.search.idxml_parsing import (
    extract_spectra_data_map,
    mod_key_from_peptidoform,
)
from mspecviz.search.models import PeptideMatch
from mspecviz.path_utils import resolve_mzml_path

logger = logging.getLogger(__name__)


def _build_run_to_metadata_map(protein_ids) -> Dict[str, str]:
    """Return {run_stem: full_path} for mzML disambiguation."""
    _, run_map_full = extract_spectra_data_map(protein_ids)
    run_to_metadata: Dict[str, str] = {}
    for full_paths in run_map_full.values():
        for full_path in full_paths:
            stem = Path(full_path).stem
            if stem not in run_to_metadata:
                run_to_metadata[stem] = full_path
    return run_to_metadata


def _build_match_from_psm(
    psm,
    idxml_path: str,
    idxml_name: str,
    run_to_metadata: Dict[str, str],
    mzml_path_strs: List[str],
) -> Optional[Tuple[str, PeptideMatch]]:
    """Build a (variant_key, PeptideMatch) from a psm_utils PSM, or None."""
    base_seq = psm.peptidoform.sequence.upper()

    mzml_filename = (psm.run + ".mzML") if psm.run else ""
    metadata_path = run_to_metadata.get(psm.run, "") if psm.run else ""

    resolved_path = ""
    if mzml_path_strs and mzml_filename:
        resolved = resolve_mzml_path(
            mzml_filename, metadata_path, mzml_path_strs, idxml_path,
        )
        if resolved:
            resolved_path = resolved

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
    mzml_path_strs: List[str],
) -> Tuple[str, List[Tuple[str, PeptideMatch]], List[str]]:
    """
    Worker: search one idXML file for target peptides using psm_utils PSM iteration.

    psm_utils resolves id_merge_index → run name internally, so each psm.run
    already contains the correct mzML basename (without extension).

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_name, list of (variant_key, PeptideMatch), list of errors)
    """
    matches: List[Tuple[str, PeptideMatch]] = []
    errors: List[str] = []

    try:
        try:
            _reader = IdXMLReader(str(idxml_path))
        except Exception as _load_err:
            errors.append(f"No peptide identifications in {idxml_name}: {_load_err}")
            return idxml_name, matches, errors

        run_to_metadata = _build_run_to_metadata_map(_reader.protein_ids)

        # Manual iteration so a parse error inside psm_utils __iter__ (e.g. float(None) on a
        # missing rescoring feature) is reported clearly. psm_utils uses a generator — once
        # it raises, iteration stops. Affected files are typically unfiltered Percolator
        # output (_pout.idXML).
        psm_iter = iter(_reader)
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
                built = _build_match_from_psm(
                    psm, idxml_path, idxml_name, run_to_metadata, mzml_path_strs,
                )
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
    """Count PSMs per base peptide sequence in one idXML file.

    Module-level so it can be submitted to a ProcessPoolExecutor.
    Returns an empty dict on any parsing error.
    """
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
    """
    Worker: extract all unique base sequences from one idXML file.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_path, set of base sequences)
    """
    peptides: Set[str] = set()
    try:
        _reader = IdXMLReader(str(idxml_path))
        for _psm in _reader:
            _base_seq = _psm.peptidoform.sequence.upper()
            if _base_seq:
                peptides.add(_base_seq)
    except Exception as e:
        logger.warning(f"_scan_single_idxml_peptides failed for {Path(idxml_path).name}: {e}")

    return idxml_path, peptides
