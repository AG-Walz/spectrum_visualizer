"""Unit tests for core.helpers (peptide parsing, time formatting, basename uniqueness)."""
import pytest

from core.helpers import (
    DuplicateBasenameError,
    basenames_shared_across_dirs,
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

    def test_duplicate_basename_in_different_dirs_allowed(self):
        # Same basename in different folders is now allowed — source folders
        # disambiguate them downstream (illumina/ vs pac/ comparison use case).
        ensure_unique_basenames(["/illumina/dup.idXML", "/pac/dup.idXML"], "idXML")

    def test_exact_same_file_twice_raises(self):
        with pytest.raises(DuplicateBasenameError, match="dup.idXML"):
            ensure_unique_basenames(["/a/dup.idXML", "/a/dup.idXML"], "idXML")

    def test_kind_interpolated_in_message(self):
        with pytest.raises(DuplicateBasenameError, match="same mzML file"):
            ensure_unique_basenames(["/a/dup.mzML", "/a/dup.mzML"], "mzML")


# ---------------------------------------------------------------------------
# TestBasenamesSharedAcrossDirs
# ---------------------------------------------------------------------------

class TestBasenamesSharedAcrossDirs:
    def test_returns_cross_dir_duplicates(self):
        shared = basenames_shared_across_dirs(
            ["/illumina/s.idXML", "/pac/s.idXML", "/illumina/only.idXML"]
        )
        assert set(shared) == {"s.idXML"}
        assert sorted(shared["s.idXML"]) == ["/illumina/s.idXML", "/pac/s.idXML"]

    def test_empty_when_all_unique(self):
        assert basenames_shared_across_dirs(["/a/x.idXML", "/b/y.idXML"]) == {}
