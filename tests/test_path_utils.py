"""
Unit tests for src/path_utils.py
"""
import pytest

from src.path_utils import deduplicate_paths, disambiguate_paths, resolve_mzml_path


# ---------------------------------------------------------------------------
# TestDeduplicatePaths
# ---------------------------------------------------------------------------

class TestDeduplicatePaths:
    @pytest.mark.parametrize("input_paths, expected", [
        (["/a", "/a", "/b"], ["/a", "/b"]),
        (["/b", "/a", "/b", "/c", "/a"], ["/b", "/a", "/c"]),
        ([], []),
        (["/x/a.mzML", "/y/b.mzML", "/z/c.mzML"], ["/x/a.mzML", "/y/b.mzML", "/z/c.mzML"]),
    ])
    def test_deduplication(self, input_paths, expected):
        assert deduplicate_paths(input_paths) == expected


# ---------------------------------------------------------------------------
# TestDisambiguatePaths
# ---------------------------------------------------------------------------

class TestDisambiguatePaths:
    def test_unique_filenames_returns_stems(self):
        paths = ["/exp1/data_a.mzML", "/exp2/data_b.mzML"]
        result = disambiguate_paths(paths)
        assert result["/exp1/data_a.mzML"] == "data_a.mzML"
        assert result["/exp2/data_b.mzML"] == "data_b.mzML"

    def test_same_filename_prepends_parent(self):
        paths = ["/exp1/data.mzML", "/exp2/data.mzML"]
        result = disambiguate_paths(paths)
        values = set(result.values())
        # Both must end with data.mzML but differ by parent
        assert len(values) == 2
        for v in values:
            assert v.endswith("data.mzML")

    def test_empty_list_returns_empty_dict(self):
        assert disambiguate_paths([]) == {}

    def test_single_path(self):
        result = disambiguate_paths(["/some/dir/file.mzML"])
        assert result == {"/some/dir/file.mzML": "file.mzML"}

    def test_three_levels_deep(self):
        # All three have the same filename and same parent -- need grandparent to disambiguate
        paths = ["/A/shared/data.mzML", "/B/shared/data.mzML", "/C/shared/data.mzML"]
        result = disambiguate_paths(paths)
        # All display names must be unique
        assert len(set(result.values())) == 3
        # Each value must start with A, B, or C
        for path, display in result.items():
            assert "/" in display  # at least two components


# ---------------------------------------------------------------------------
# TestResolveMzmlPath
# ---------------------------------------------------------------------------

class TestResolveMzmlPath:
    def test_single_candidate_exact_match(self):
        candidates = ["/data/run1.mzML"]
        result = resolve_mzml_path("run1.mzML", "/data/run1.mzML", candidates)
        assert result == "/data/run1.mzML"

    def test_no_candidate_returns_none(self):
        candidates = ["/data/other.mzML"]
        result = resolve_mzml_path("run1.mzML", "/data/run1.mzML", candidates)
        assert result is None

    def test_case_insensitive_match(self):
        candidates = ["/data/data.mzml"]
        result = resolve_mzml_path("Data.mzML", "/path/Data.mzML", candidates)
        assert result == "/data/data.mzml"

    def test_multiple_candidates_suffix_match(self):
        candidates = ["/exp1/project/run1.mzML", "/exp2/project/run1.mzML"]
        # metadata path shares /exp1/project/ prefix with first candidate
        result = resolve_mzml_path(
            "run1.mzML",
            "/exp1/project/run1.mzML",
            candidates,
        )
        assert result == "/exp1/project/run1.mzML"

    def test_multiple_candidates_idxml_proximity(self, tmp_path):
        # Two candidates with the same name, no metadata suffix match
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()
        (dir1 / "run.mzML").write_text("")
        (dir2 / "run.mzML").write_text("")
        idxml_path = str(dir1 / "search.idXML")
        candidates = [str(dir1 / "run.mzML"), str(dir2 / "run.mzML")]
        result = resolve_mzml_path(
            "run.mzML",
            "/completely/different/machine/run.mzML",  # no suffix match
            candidates,
            idxml_path=idxml_path,
        )
        assert result == str(dir1 / "run.mzML")

    def test_multiple_candidates_ultimate_fallback(self):
        candidates = ["/x/run1.mzML", "/y/run1.mzML"]
        # No metadata suffix match (no overlapping directory components)
        # No idxml_path proximity match
        result = resolve_mzml_path(
            "run1.mzML",
            "/completely/different/path/run1.mzML",
            candidates,
            idxml_path=None,
        )
        # Must return one of the candidates
        assert result in candidates
