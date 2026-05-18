"""mzML native-ID parsing helpers — regexes + canonicalisation for ``scan=N`` / ``index=N`` / ``spectrum=N`` / Thermo ``controllerType=…`` (used by ``mzml_lookup``)."""
import re
from typing import Optional

SCAN_TAG = "scan="
INDEX_TAG = "index="
SPECTRUM_TAG = "spectrum="
NATIVE_ID_TAGS = (SCAN_TAG, INDEX_TAG, SPECTRUM_TAG)

NATIVE_ID_NUMBER_PATTERN = re.compile(r'(?:scan|index|spectrum)=(\d+)')
BARE_DIGITS_PATTERN = re.compile(r'(\d+)')
_TAGGED_TOKEN_PATTERN = re.compile(r'(?:scan|index|spectrum)=\S+')


def extract_id_number(native_id: str) -> Optional[int]:
    """Pull the integer out of a native ID (``scan=42`` -> ``42``); falls back to first digit run."""
    keyed_match = NATIVE_ID_NUMBER_PATTERN.search(native_id)
    if keyed_match:
        return int(keyed_match.group(1))
    digit_match = BARE_DIGITS_PATTERN.search(native_id)
    if digit_match:
        return int(digit_match.group(1))
    return None


def clean_spectrum_id(spectrum_id: str) -> str:
    """Canonicalise a spectrum ID to its first ``(scan|index|spectrum)=N`` token.

    Thermo ``controllerType=…`` IDs pass through unchanged (their ``scan=N`` is a sub-fragment)."""
    if not spectrum_id:
        return ""
    clean_id = spectrum_id.strip()
    if "controllerType" in clean_id:
        return clean_id
    token = _TAGGED_TOKEN_PATTERN.search(clean_id)
    return token.group(0) if token else clean_id