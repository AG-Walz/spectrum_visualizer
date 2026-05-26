"""mzML native-ID parsing helpers for mass-spectrometry spectrum identifiers.
Provides regexes plus two functions that turn raw mzML nativeID strings into a
canonical form and a plain integer 
(imported by core.spectra.spectrum_streaming for its scan-vs-index detection, ID matching
and binary search, and streamlit_app/views/explore/mirror_plots.py to derive a scan value for display.)
"""
import re
from typing import Optional

# Native-ID tag prefixes — the three conventions an mzML nativeID may use to key a
# spectrum number. Consumed by spectrum_streaming to decide which numbering scheme a run uses.
SCAN_TAG = "scan="
INDEX_TAG = "index="
SPECTRUM_TAG = "spectrum="
NATIVE_ID_TAGS = (SCAN_TAG, INDEX_TAG, SPECTRUM_TAG)

# Captures the integer that immediately follows a scan= / index= / spectrum= key.
NATIVE_ID_NUMBER_PATTERN = re.compile(r'(?:scan|index|spectrum)=(\d+)')
# Fallback — the first run of digits anywhere in the string, used when no keyed tag is present.
BARE_DIGITS_PATTERN = re.compile(r'(\d+)')
# Matches a whole keyed token (key plus its non-whitespace value) for canonicalisation.
_TAGGED_TOKEN_PATTERN = re.compile(r'(?:scan|index|spectrum)=\S+')


def extract_id_number(native_id: str) -> Optional[int]:
    """Pull the spectrum number out of an mzML native ID.
    Returns the integer keyed by scan= / index= / spectrum. When no keyed tag is present, 
    falls back to the first run of digits anywhere in the string; returns None when the ID contains no digits at all.
    The returned value carries no base assumption — an index=N number may be
    0-based or 1-based depending on the instrument vendor (see module docstring).
    """
    # Prefer a keyed match — the digits attached to a recognised tag are unambiguous.
    keyed_match = NATIVE_ID_NUMBER_PATTERN.search(native_id)
    if keyed_match:
        return int(keyed_match.group(1))
    # No keyed tag — accept any bare digit run as a best-effort spectrum number.
    digit_match = BARE_DIGITS_PATTERN.search(native_id)
    if digit_match:
        return int(digit_match.group(1))
    return None


def clean_spectrum_id(spectrum_id: str) -> str:
    """Canonicalise a spectrum ID down to a single comparable token.
    Returns the first whole (scan|index|spectrum)=N token found in the ID, so that
    IDs carrying extra surrounding fields still compare equal on the part that matters.
    An empty or falsy input returns ""; an ID with no recognised tag is returned
    stripped but otherwise unchanged.
    Thermo controllerType=… IDs are passed through whole — there the scan=N is a
    sub-fragment of a larger compound key and must not be isolated, or distinct spectra
    would collapse onto the same canonical ID.
    """
    if not spectrum_id:
        return ""
    clean_id = spectrum_id.strip()
    # Thermo compound IDs — keep intact; the scan=N part alone does not identify the spectrum.
    if "controllerType" in clean_id:
        return clean_id
    # Reduce to the first keyed token; if none exists, hand back the stripped ID as-is.
    token = _TAGGED_TOKEN_PATTERN.search(clean_id)
    return token.group(0) if token else clean_id