"""mzML native-ID parsing helpers.

mzML converters emit spectrum IDs in a handful of formats (``scan=N``,
``index=N``, ``spectrum=N``, plus Thermo's ``controllerType=…``). Centralising
the regexes + cleaning logic here keeps ``mzml_lookup`` focused on indexing.
"""
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
    """Extract the numeric value from a native ID string."""
    keyed_match = NATIVE_ID_NUMBER_PATTERN.search(native_id)
    if keyed_match:
        return int(keyed_match.group(1))
    digit_match = BARE_DIGITS_PATTERN.search(native_id)
    if digit_match:
        return int(digit_match.group(1))
    return None


def clean_spectrum_id(spectrum_id: str) -> str:
    """Canonicalise a spectrum ID for lookup.

    Thermo IDs (``controllerType=…``) are returned as-is — their embedded
    ``scan=N`` is a sub-fragment, not the whole ID. Otherwise we strip down to
    the first ``(scan|index|spectrum)=N`` token.
    """
    if not spectrum_id:
        return ""
    clean_id = spectrum_id.strip()
    if "controllerType" in clean_id:
        return clean_id
    token = _TAGGED_TOKEN_PATTERN.search(clean_id)
    return token.group(0) if token else clean_id