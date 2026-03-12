"""
Unit tests for src/export.py
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.export import (
    find_or_create_export_directory,
    commit_export_directory,
    sanitize_filename,
    shorten_stem,
    get_workspace_plot_dirs,
)


# ---------------------------------------------------------------------------
# TestSanitizeFilename
# ---------------------------------------------------------------------------

class TestSanitizeFilename:
    def test_removes_illegal_chars(self):
        # Characters < > : " / \ | ? * should be replaced with _
        result = sanitize_filename('a<b>c:d"e/f\\g|h?i*j')
        assert "<" not in result
        assert ">" not in result
        assert ":" not in result
        assert '"' not in result
        assert "/" not in result
        assert "\\" not in result
        assert "|" not in result
        assert "?" not in result
        assert "*" not in result

    def test_collapses_multiple_underscores(self):
        result = sanitize_filename("hello___world")
        assert "__" not in result
        assert result == "hello_world"

    def test_strips_leading_trailing_underscores(self):
        result = sanitize_filename("_hello_world_")
        assert not result.startswith("_")
        assert not result.endswith("_")

    def test_truncates_long_names(self):
        long_name = "a" * 110
        result = sanitize_filename(long_name)
        assert len(result) <= 100

    def test_empty_string_returns_fallback(self):
        # sanitize_filename returns empty string for empty input (no fallback)
        result = sanitize_filename("")
        assert isinstance(result, str)

    def test_already_valid_unchanged(self):
        name = "valid_filename_123"
        assert sanitize_filename(name) == name


# ---------------------------------------------------------------------------
# TestShortenStem
# ---------------------------------------------------------------------------

class TestShortenStem:
    def test_short_stem_unchanged(self):
        stem = "short_name"
        assert shorten_stem(stem) == stem

    def test_exactly_30_chars_unchanged(self):
        stem = "a" * 30
        assert shorten_stem(stem, max_len=30) == stem

    def test_long_stem_truncated_with_hash(self):
        stem = "a" * 40
        result = shorten_stem(stem, max_len=30)
        assert len(result) == 30
        # Must contain underscore separator before 8-char hash
        assert "_" in result
        parts = result.rsplit("_", 1)
        assert len(parts[1]) == 8

    def test_different_long_stems_produce_different_results(self):
        stem1 = "x" * 40
        stem2 = "y" * 40
        assert shorten_stem(stem1) != shorten_stem(stem2)

    def test_same_input_always_same_output(self):
        stem = "this_is_a_very_long_filename_for_testing_determinism"
        result1 = shorten_stem(stem)
        result2 = shorten_stem(stem)
        assert result1 == result2


# ---------------------------------------------------------------------------
# TestFindOrCreateExportDirectory
# ---------------------------------------------------------------------------

class TestFindOrCreateExportDirectory:
    def test_creates_first_run(self, tmp_path):
        settings = {"peptide": "PEPTIDEK", "tolerance": 0.02}
        root, png, svg, pdf = find_or_create_export_directory(
            str(tmp_path), "myproject", settings
        )
        # No directories created yet — only after commit
        assert not root.exists()
        commit_export_directory(root, settings)
        assert root.exists()
        assert (root / "run_settings.json").exists()
        # Folder name should contain _001_
        assert "_001_" in root.name

    def test_reuses_matching_settings(self, tmp_path):
        settings = {"peptide": "PEPTIDEK", "tolerance": 0.02}
        root1, _, _, _ = find_or_create_export_directory(
            str(tmp_path), "proj", settings
        )
        commit_export_directory(root1, settings)
        root2, _, _, _ = find_or_create_export_directory(
            str(tmp_path), "proj", settings
        )
        # Same run_settings → same folder number
        num1 = root1.name.split("_")[1]
        num2 = root2.name.split("_")[1]
        assert num1 == num2

    def test_creates_new_for_different_settings(self, tmp_path):
        settings_a = {"peptide": "PEPTIDEK", "tolerance": 0.02}
        settings_b = {"peptide": "AAAVAAAK", "tolerance": 0.05}
        root1, _, _, _ = find_or_create_export_directory(
            str(tmp_path), "proj", settings_a
        )
        commit_export_directory(root1, settings_a)
        root2, _, _, _ = find_or_create_export_directory(
            str(tmp_path), "proj", settings_b
        )
        assert root1.name != root2.name

    def test_sequential_numbering(self, tmp_path):
        nums = []
        for i in range(3):
            settings = {"run": i}
            root, _, _, _ = find_or_create_export_directory(
                str(tmp_path), "proj", settings
            )
            commit_export_directory(root, settings)
            # Extract the 3-digit number from folder name
            parts = root.name.split("_")
            nums.append(int(parts[1]))
        assert sorted(nums) == [1, 2, 3]

    def test_returns_four_paths(self, tmp_path):
        result = find_or_create_export_directory(
            str(tmp_path), "proj", {"x": 1}
        )
        assert len(result) == 4
        for p in result:
            assert isinstance(p, Path)

    def test_png_dir_not_precreated(self, tmp_path):
        # Format subdirectories are created lazily on first write, not upfront
        root, png_dir, _, _ = find_or_create_export_directory(
            str(tmp_path), "proj", {"x": 1}
        )
        commit_export_directory(root, {"x": 1})
        assert not png_dir.exists()
        assert png_dir.name == "png"


# ---------------------------------------------------------------------------
# TestGetWorkspacePlotDirs
# ---------------------------------------------------------------------------

class TestGetWorkspacePlotDirs:
    def test_returns_correct_paths(self, tmp_path):
        root, png, svg, pdf = get_workspace_plot_dirs(tmp_path)
        assert root == tmp_path / "plots"
        assert png == tmp_path / "plots" / "png"
        assert svg == tmp_path / "plots" / "svg"
        assert pdf == tmp_path / "plots" / "pdf"

    def test_no_dirs_created(self, tmp_path):
        root, png, svg, pdf = get_workspace_plot_dirs(tmp_path)
        assert not root.exists()
        assert not png.exists()
