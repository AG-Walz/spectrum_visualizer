"""
Path disambiguation and deduplication utilities.

Handles files with the same name from different directories by computing
shortest unique display names (e.g., "exp1/data.mzML" vs "exp2/data.mzML").
"""
import logging
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def deduplicate_paths(paths: List[str]) -> List[str]:
    """
    Remove exact duplicate paths, preserving order.

    Args:
        paths: List of file path strings

    Returns:
        List with duplicates removed (first occurrence kept)
    """
    return list(dict.fromkeys(paths))


def disambiguate_paths(paths: List[str]) -> Dict[str, str]:
    """
    Map full paths to shortest unique display names.

    For unique filenames, returns just the filename.
    For duplicates, prepends parent directories until all names are unique.
    Exact duplicate paths are deduplicated first.

    Examples:
        ["/a/exp1/data.mzML", "/a/exp2/data.mzML", "/a/exp1/other.mzML"]
        -> {
            "/a/exp1/data.mzML": "exp1/data.mzML",
            "/a/exp2/data.mzML": "exp2/data.mzML",
            "/a/exp1/other.mzML": "other.mzML",
        }

    Args:
        paths: List of full file path strings

    Returns:
        Dict mapping full_path -> shortest unique display name
    """
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
    """
    Resolve a bare mzML filename to the best matching full path.

    Uses multiple strategies:
    1. If only one selected path matches the bare name -> use it
    2. If multiple match -> use metadata full path for suffix matching
    3. Fallback -> directory proximity to the idXML file

    Args:
        mzml_bare_name: Bare filename from idXML metadata (e.g., "data.mzML")
        metadata_full_path: Original full path from idXML metadata (may be
            from a different machine, but directory structure often matches)
        mzml_paths: List of user-selected mzML full paths
        idxml_path: Full path of the idXML file (for directory proximity)

    Returns:
        Best matching full path from mzml_paths, or None if no match
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
