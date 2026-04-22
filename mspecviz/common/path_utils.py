"""Path dedup + shortest-unique-display-name helpers (``exp1/data.mzML`` vs ``exp2/data.mzML``)."""
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def disambiguate_paths(paths: List[str]) -> Dict[str, str]:
    """Map each full path to the shortest unique display name (prepend parents until unique)."""
    unique_paths = list(dict.fromkeys(paths))

    if not unique_paths:
        return {}

    # Split each path into its components
    path_parts = {p: Path(p).parts for p in unique_paths}

    result = {}
    unresolved = set(unique_paths)
    depth = 1

    while unresolved:
        # Build display name using the last `depth` components of each path
        current_names = {}
        for p in unresolved:
            parts = path_parts[p]
            n = min(depth, len(parts))
            current_names[p] = str(Path(*parts[-n:]))

        # Find which names are now unique at this depth
        name_counts = Counter(current_names.values())
        newly_resolved = set()

        for p in unresolved:
            if name_counts[current_names[p]] == 1:
                result[p] = current_names[p]
                newly_resolved.add(p)

        unresolved -= newly_resolved
        depth += 1

        # Safety: if we've exhausted all path components, use the full path
        if unresolved:
            max_depth = max(len(path_parts[p]) for p in unresolved)
            if depth > max_depth:
                for p in unresolved:
                    result[p] = str(Path(*path_parts[p]))
                break

    return result


def resolve_mzml_path(
    mzml_bare_name: str,
    metadata_full_path: str,
    mzml_paths: List[str],
    idxml_path: Optional[str] = None,
) -> Optional[str]:
    """Resolve a bare mzML filename to the best full path in ``mzml_paths``.

    Strategy: (1) unique bare-name match, (2) longest trailing-component match
    against ``metadata_full_path``, (3) directory proximity to ``idxml_path``.
    Returns ``None`` when nothing matches.
    """
    if not mzml_bare_name or not mzml_paths:
        return None

    # Build candidates by bare name match (case-insensitive)
    bare_lower = mzml_bare_name.lower()
    candidates = [p for p in mzml_paths if Path(p).name.lower() == bare_lower]

    if not candidates:
        return None

    if len(candidates) == 1:
        return candidates[0]

    # Multiple candidates — use metadata path for suffix matching
    if metadata_full_path:
        meta_parts = Path(metadata_full_path).parts
        best_match = None
        best_depth = 0

        for cand in candidates:
            cand_parts = Path(cand).parts
            # Count how many trailing path components match
            match_depth = 0
            for m, c in zip(reversed(meta_parts), reversed(cand_parts)):
                if m.lower() == c.lower():
                    match_depth += 1
                else:
                    break

            if match_depth > best_depth:
                best_depth = match_depth
                best_match = cand

        # Accept if we matched more than just the filename
        if best_match and best_depth > 1:
            return best_match

    # Fallback: directory proximity to idXML
    if idxml_path:
        idxml_dir = str(Path(idxml_path).parent)
        for cand in candidates:
            if str(Path(cand).parent) == idxml_dir:
                return cand

    # Ultimate fallback: first candidate
    logger.warning(
        f"Ambiguous mzML resolution for '{mzml_bare_name}': "
        f"{len(candidates)} candidates, using first match"
    )
    return candidates[0]
