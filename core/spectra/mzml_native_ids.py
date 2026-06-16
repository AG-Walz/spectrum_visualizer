"""mzML native-ID parsing helper for mass-spectrometry spectrum identifiers.
Provides regexes plus extract_id_number, which turns a raw mzML nativeID string into a plain
integer (imported by core.spectra.spectrum_streaming for its scan-number fallback lookup, and
streamlit_app/views/explore/mirror_plots.py to derive a scan value for display).
"""
import re
from typing import Optional

# Captures the integer that immediately follows a scan= / index= / spectrum= key.
NATIVE_ID_NUMBER_PATTERN = re.compile(r'(?:scan|index|spectrum)=(\d+)')
# Fallback — the first run of digits anywhere in the string, used when no keyed tag is present.
BARE_DIGITS_PATTERN = re.compile(r'(\d+)')


def extract_id_number(native_id: str) -> Optional[int]:
    """Pull the spectrum number out of an mzML native ID.
    Returns the integer keyed by scan= / index= / spectrum=, falling back to the first
    run of digits, or None when the ID has no digits. The value carries no base
    assumption — an index=N number may be 0- or 1-based per vendor.
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
