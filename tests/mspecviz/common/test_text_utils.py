"""Unit tests for mspecviz/common/text_utils.py"""
import pytest

from mspecviz.common.helpers import (
    DuplicateBasenameError,
    ensure_unique_basenames,
    format_time,
    parse_peptide_list,
)


# ---------------------------------------------------------------------------
# TestParsePeptideList
# ---------------------------------------------------------------------------

class TestParsePeptideList:
    @pytest.mark.parametrize("text, expected", [
        ("", []),
        ("PEPTIDE", ["PEPTIDE"]),
        ("PEPA\nPEPB\nPEPC", ["PEPA", "PEPB", "PEPC"]),
        ("  pepa  \n\n PEPB \n", ["PEPA", "PEPB"]),
        ("pepa\r\npepb", ["PEPA", "PEPB"]),
        ("PEPA, PEPB", ["PEPA, PEPB"]),
    ], ids=["empty", "single", "multiple", "strips_and_drops_blank", "crlf", "comma_not_separator"])
    def test_parse(self, text, expected):
        assert parse_peptide_list(text) == expected


# ---------------------------------------------------------------------------
# TestFormatTime
# ---------------------------------------------------------------------------

class TestFormatTime:
    @pytest.mark.parametrize("seconds, expected", [
        (0, "0s"),
        (30, "30s"),
        (59, "59s"),
        (60, "1m 0s"),
        (90, "1m 30s"),
        (3599, "59m 59s"),
        (3600, "1h 0m"),
        (3725, "1h 2m"),
        (7200, "2h 0m"),
    ], ids=["zero", "sub_minute", "boundary_59s", "exact_minute",
            "minute_and_seconds", "boundary_59m_59s", "exact_hour",
            "hour_and_minutes", "multi_hour"])
    def test_format(self, seconds, expected):
        assert format_time(seconds) == expected


# ---------------------------------------------------------------------------
# TestEnsureUniqueBasenames
# ---------------------------------------------------------------------------

class TestEnsureUniqueBasenames:
    def test_all_unique_does_not_raise(self):
        ensure_unique_basenames(["/a/x.idXML", "/b/y.idXML", "/c/z.idXML"], "idXML")

    def test_empty_does_not_raise(self):
        ensure_unique_basenames([], "idXML")

    def test_duplicate_basename_in_different_dirs_raises(self):
        with pytest.raises(DuplicateBasenameError, match="dup.idXML"):
            ensure_unique_basenames(["/a/dup.idXML", "/b/dup.idXML"], "idXML")

    def test_message_lists_all_duplicates_sorted(self):
        with pytest.raises(DuplicateBasenameError) as exc:
            ensure_unique_basenames(
                ["/a/banana.mzML", "/x/apple.mzML", "/b/banana.mzML", "/y/apple.mzML"],
                "mzML",
            )
        assert "['apple.mzML', 'banana.mzML']" in str(exc.value)

    def test_kind_interpolated_in_message(self):
        with pytest.raises(DuplicateBasenameError, match="selected mzML files"):
            ensure_unique_basenames(["/a/dup.mzML", "/b/dup.mzML"], "mzML")

    def test_is_value_error_subclass(self):
        assert issubclass(DuplicateBasenameError, ValueError)
