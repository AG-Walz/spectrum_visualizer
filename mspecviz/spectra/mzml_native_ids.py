"""mzML native-ID parsing helpers.

mzML converters emit spectrum IDs in a handful of formats (``scan=N``,
``index=N``, ``spectrum=N``, plus Thermo's ``controllerType=…``). Centralising
the regexes + cleaning logic here keeps ``mzml_lookup`` focused on indexing.
"""
import re
from typing import Optional

# --- Native-ID tag constants ------------------------------------------------
# mzML converters emit spectrum IDs in a handful of formats. Listing the tags
# here (and reusing them in regex + substring checks) means adding a new
# dialect touches a single location.
SCAN_TAG = "scan="
INDEX_TAG = "index="
SPECTRUM_TAG = "spectrum="
NATIVE_ID_TAGS = (SCAN_TAG, INDEX_TAG, SPECTRUM_TAG)

# Regex matching a ``(scan|index|spectrum)=N`` fragment inside a native ID.
NATIVE_ID_NUM_RE = re.compile(r'(?:scan|index|spectrum)=(\d+)')
# Fallback regex: any digit run.
ANY_DIGITS_RE = re.compile(r'(\d+)')


def _extract_id_number(native_id: str) -> Optional[int]:
    """Extract the numeric value from a native ID string."""
    keyed_match = NATIVE_ID_NUM_RE.search(native_id)
    if keyed_match:
        return int(keyed_match.group(1))
    digit_match = ANY_DIGITS_RE.search(native_id)
    if digit_match:
        return int(digit_match.group(1))
    return None


def clean_spectrum_id(spectrum_id: str) -> str:
    """Clean a spectrum ID for lookup (strips whitespace, extracts canonical form)."""
    if not spectrum_id:
        return ""

    if isinstance(spectrum_id, bytes):
        spectrum_id = spectrum_id.decode()

    clean_id = spectrum_id.strip()

    # If it contains "controllerType" (Thermo format), keep it
    if "controllerType" in clean_id:
        return clean_id

    # For simple scan=N format, ensure clean extraction
    for tag in NATIVE_ID_TAGS:
        if tag in clean_id:
            idx = clean_id.find(tag)
            rest = clean_id[idx:]
            parts = rest.split()
            return parts[0]

    return clean_id
