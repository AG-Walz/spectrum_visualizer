"""
Unit tests for mspecviz/export/directories.py
"""
from pathlib import Path

import pytest

from mspecviz.export.directories import (
    find_or_create_export_directory,
    commit_export_directory,
    sanitize_filename,
)


# ---------------------------------------------------------------------------
# TestSanitizeFilename
# ---------------------------------------------------------------------------

class TestSanitizeFilename:
    def test_removes_illegal_chars(self):
        result = sanitize_filename('a<b>c:d"e/f\\g|h?i*j')
        for char in '<>:"/\\|?*':
            assert char not in result

    @pytest.mark.parametrize("input_name, check", [
        ("hello___world", lambda r: "__" not in r and r == "hello_world"),
        ("_hello_world_", lambda r: not r.startswith("_") and not r.endswith("_")),
        ("a" * 110, lambda r: len(r) <= 100),
        ("", lambda r: isinstance(r, str)),
        ("valid_filename_123", lambda r: r == "valid_filename_123"),
    ], ids=["collapses_underscores", "strips_edges", "truncates_long", "empty_fallback", "valid_unchanged"])
    def test_sanitize_rules(self, input_name, check):
        assert check(sanitize_filename(input_name))


# ---------------------------------------------------------------------------
# TestFindOrCreateExportDirectory
# ---------------------------------------------------------------------------

class TestFindOrCreateExportDirectory:
    def test_creates_first_run(self, tmp_path):
        settings = {"peptide": "PEPTIDEK", "tolerance": 0.02}
        root, png, svg, pdf = find_or_create_export_directory(
            str(tmp_path), "myproject", settings
        )
        # No directories created yet -- only after commit
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
        # Same run_settings -> same folder number
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
