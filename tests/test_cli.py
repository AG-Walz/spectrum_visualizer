"""Unit tests for CLI pure-Python helpers (cli.py).

All tests use synthetic data — no real mzML/idXML files needed.
"""
import argparse
import tempfile
from pathlib import Path

import pytest
import yaml

# Ensure project root is importable
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cli import (
    CLIConfig,
    _expand_file_globs,
    _filter_files_by_substring,
    _load_peptides,
    _strip_modifications,
    create_parser,
    load_config,
)


# ---------------------------------------------------------------------------
# _strip_modifications
# ---------------------------------------------------------------------------

class TestStripModifications:
    def test_single_mod(self):
        assert _strip_modifications("ATYSVM(Oxidation)QK") == "ATYSVMQK"

    def test_multiple_mods(self):
        assert _strip_modifications("M(Ox)PEPTM(Ox)IDE") == "MPEPTMIDE"

    def test_no_mods(self):
        assert _strip_modifications("LGGNEQVTR") == "LGGNEQVTR"

    def test_empty(self):
        assert _strip_modifications("") == ""


# ---------------------------------------------------------------------------
# _load_peptides
# ---------------------------------------------------------------------------

class TestLoadPeptides:
    def test_from_file(self, tmp_path):
        pep_file = tmp_path / "peptides.txt"
        pep_file.write_text("LGGNEQVTR\nGAGSSEPVTGLDAK\natysvgmqk\n")
        result = _load_peptides(str(pep_file))
        assert result == ["LGGNEQVTR", "GAGSSEPVTGLDAK", "ATYSVGMQK"]

    def test_from_file_strips_mods(self, tmp_path):
        pep_file = tmp_path / "peptides.txt"
        pep_file.write_text("ATYSVM(Oxidation)QK\nLGGNEQVTR\n")
        result = _load_peptides(str(pep_file))
        assert "ATYSVMQK" in result
        assert "LGGNEQVTR" in result

    def test_from_file_deduplicates(self, tmp_path):
        pep_file = tmp_path / "peptides.txt"
        pep_file.write_text("AAA\naaa\nBBB\nAAA\n")
        result = _load_peptides(str(pep_file))
        assert result == ["AAA", "BBB"]

    def test_from_inline_list(self):
        result = _load_peptides(["PEP1", "pep2", "PEP1"])
        assert result == ["PEP1", "PEP2"]

    def test_all_keyword(self):
        assert _load_peptides("all") == []
        assert _load_peptides("ALL") == []
        assert _load_peptides("  All  ") == []

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            _load_peptides("/nonexistent/peptides.txt")

    def test_empty_lines_skipped(self, tmp_path):
        pep_file = tmp_path / "peptides.txt"
        pep_file.write_text("\n\nLGGNEQVTR\n\n\n")
        result = _load_peptides(str(pep_file))
        assert result == ["LGGNEQVTR"]


# ---------------------------------------------------------------------------
# _filter_files_by_substring
# ---------------------------------------------------------------------------

class TestFilterFilesBySubstring:
    def test_empty_filter_returns_all(self, tmp_path):
        files = ["/data/run1.mzML", "/data/run2.mzML"]
        assert _filter_files_by_substring(files, "") == files

    def test_case_insensitive(self):
        files = ["/data/Standard_001.mzML", "/data/query_001.mzML"]
        assert _filter_files_by_substring(files, "standard") == ["/data/Standard_001.mzML"]
        assert _filter_files_by_substring(files, "STANDARD") == ["/data/Standard_001.mzML"]

    def test_matches_full_path_including_directory(self):
        files = ["/run1/data.mzML", "/data/run2.mzML"]
        # "run1" appears in the directory component of the first file — should match
        result = _filter_files_by_substring(files, "run1")
        assert result == ["/run1/data.mzML"]

    def test_no_match_returns_empty(self):
        files = ["/data/run1.mzML", "/data/run2.mzML"]
        assert _filter_files_by_substring(files, "xyz") == []

    def test_partial_match(self):
        files = ["/data/ProCa_run1.idXML", "/data/ProCa_run2.idXML", "/data/other.idXML"]
        result = _filter_files_by_substring(files, "ProCa")
        assert len(result) == 2
        assert "/data/other.idXML" not in result


# ---------------------------------------------------------------------------
# _expand_file_globs
# ---------------------------------------------------------------------------

class TestExpandFileGlobs:
    def test_glob_expansion(self, tmp_path):
        for name in ["a.idXML", "b.idXML", "c.mzML"]:
            (tmp_path / name).touch()
        result = _expand_file_globs([str(tmp_path / "*.idXML")])
        assert len(result) == 2
        assert all(p.endswith(".idXML") for p in result)

    def test_literal_path(self, tmp_path):
        f = tmp_path / "test.idXML"
        f.touch()
        result = _expand_file_globs([str(f)])
        assert len(result) == 1
        assert result[0] == str(f.resolve())

    def test_deduplication(self, tmp_path):
        f = tmp_path / "test.idXML"
        f.touch()
        result = _expand_file_globs([str(f), str(f)])
        assert len(result) == 1

    def test_no_match_keeps_literal(self):
        result = _expand_file_globs(["/nonexistent/fake.idXML"])
        assert len(result) == 1


# ---------------------------------------------------------------------------
# CLIConfig defaults
# ---------------------------------------------------------------------------

class TestCLIConfigDefaults:
    def test_default_values(self):
        cfg = CLIConfig()
        assert cfg.ion_types == ["b", "y"]
        assert cfg.mz_tolerance == 0.02
        assert cfg.export_formats == ["png"]
        assert cfg.dpi == 300
        assert cfg.workers == 4
        assert cfg.cache_enabled is True
        assert cfg.normalize is True

    def test_raw_intensity(self):
        cfg = CLIConfig(intensity_mode="Raw")
        assert cfg.normalize is False


# ---------------------------------------------------------------------------
# create_parser
# ---------------------------------------------------------------------------

class TestCreateParser:
    def test_config_required(self):
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])  # --config is required

    def test_basic_args(self):
        parser = create_parser()
        args = parser.parse_args(["--config", "test.yaml"])
        assert args.config == "test.yaml"

    def test_annotation_flags(self):
        parser = create_parser()
        args = parser.parse_args([
            "--config", "test.yaml",
            "--mz-tolerance", "0.05",
            "--add-losses",
            "--ion-types", "a", "b", "y",
        ])
        assert args.mz_tolerance == 0.05
        assert args.add_losses is True
        assert args.ion_types == ["a", "b", "y"]

    def test_verbosity_flags(self):
        parser = create_parser()
        args = parser.parse_args(["--config", "t.yaml", "--verbose"])
        assert args.verbose is True

        args = parser.parse_args(["--config", "t.yaml", "--quiet"])
        assert args.quiet is True


# ---------------------------------------------------------------------------
# load_config — 3-layer merge
# ---------------------------------------------------------------------------

class TestLoadConfig:
    def _make_yaml(self, tmp_path, data):
        yaml_path = tmp_path / "config.yaml"
        yaml_path.write_text(yaml.dump(data))
        return str(yaml_path)

    def _make_cli_args(self, **overrides):
        """Build a Namespace mimicking argparse output with all defaults None/False."""
        defaults = dict(
            config="dummy.yaml",
            idxml=None, idxml_filter=None, mzml=None, mzml_filter=None,
            parquet=None, peptides=None,
            reference_substring=None,
            ion_types=None, mz_tolerance=None,
            add_losses=False, no_first_prefix_ion=False,
            remove_precursor_peaks=False, raw_intensity=False,
            export_dir=None, project_name=None, formats=None,
            dpi=None, plot_unpaired=False,
            no_cache=False, workers=None, no_tmpfs=False,
            verbose=False, quiet=False,
        )
        defaults.update(overrides)
        return argparse.Namespace(**defaults)

    def test_yaml_overrides_defaults(self, tmp_path):
        yaml_path = self._make_yaml(tmp_path, {
            "annotation": {"mz_tolerance": 0.1, "ion_types": ["a", "b"]},
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert cfg.mz_tolerance == 0.1
        assert cfg.ion_types == ["a", "b"]

    def test_cli_overrides_yaml(self, tmp_path):
        yaml_path = self._make_yaml(tmp_path, {
            "annotation": {"mz_tolerance": 0.1},
        })
        cli_args = self._make_cli_args(mz_tolerance=0.05)
        cfg = load_config(yaml_path, cli_args)
        assert cfg.mz_tolerance == 0.05

    def test_export_settings(self, tmp_path):
        yaml_path = self._make_yaml(tmp_path, {
            "export": {"project_name": "test_proj", "dpi": 150},
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert cfg.project_name == "test_proj"
        assert cfg.dpi == 150

    def test_peptides_all(self, tmp_path):
        yaml_path = self._make_yaml(tmp_path, {"peptides": "all"})
        cfg = load_config(yaml_path, self._make_cli_args())
        assert cfg.peptides == []

    def test_reference_substring(self, tmp_path):
        yaml_path = self._make_yaml(tmp_path, {
            "reference": {"substring": "Reference"},
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert cfg.reference_substring == "Reference"

    def test_yaml_idxml_filter(self, tmp_path):
        # Create two idXML files; filter should keep only the one matching "run1"
        (tmp_path / "run1_sample.idXML").touch()
        (tmp_path / "run2_sample.idXML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {
                "idxml_files": [str(tmp_path / "*.idXML")],
                "idxml_filter": "run1",
            },
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert len(cfg.idxml_files) == 1
        assert "run1" in cfg.idxml_files[0]

    def test_yaml_mzml_filter(self, tmp_path):
        (tmp_path / "Standard_001.mzML").touch()
        (tmp_path / "query_001.mzML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {
                "mzml_files": [str(tmp_path / "*.mzML")],
                "mzml_filter": "Standard",
            },
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert len(cfg.mzml_files) == 1
        assert "Standard" in cfg.mzml_files[0]

    def test_cli_filter_overrides_yaml(self, tmp_path):
        # YAML has no filter; CLI narrows files already loaded from YAML
        (tmp_path / "run1.idXML").touch()
        (tmp_path / "run2.idXML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {"idxml_files": [str(tmp_path / "*.idXML")]},
        })
        cfg = load_config(yaml_path, self._make_cli_args(idxml_filter="run1"))
        assert len(cfg.idxml_files) == 1
        assert "run1" in cfg.idxml_files[0]

    def test_cli_filter_applies_when_idxml_also_given(self, tmp_path):
        (tmp_path / "run1.mzML").touch()
        (tmp_path / "run2.mzML").touch()
        yaml_path = self._make_yaml(tmp_path, {})
        cfg = load_config(
            yaml_path,
            self._make_cli_args(
                mzml=[str(tmp_path / "*.mzML")],
                mzml_filter="run2",
            ),
        )
        assert len(cfg.mzml_files) == 1
        assert "run2" in cfg.mzml_files[0]

    def test_yaml_filter_with_cli_files(self, tmp_path):
        # YAML sets the filter; CLI provides --idxml (not YAML idxml_files).
        # Filter must apply to the CLI-expanded files.
        (tmp_path / "run1.idXML").touch()
        (tmp_path / "run2.idXML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {"idxml_filter": "run1"},  # filter only, no files
        })
        cfg = load_config(
            yaml_path,
            self._make_cli_args(idxml=[str(tmp_path / "*.idXML")]),
        )
        assert len(cfg.idxml_files) == 1
        assert "run1" in cfg.idxml_files[0]

    def test_cli_filter_overrides_yaml_filter(self, tmp_path):
        # YAML has filter="run1" and files; CLI changes filter to "run2".
        # The final file list must reflect the CLI filter applied to the full
        # (unfiltered) glob, not "run2 applied to already-run1-filtered files".
        (tmp_path / "run1.idXML").touch()
        (tmp_path / "run2.idXML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {
                "idxml_files": [str(tmp_path / "*.idXML")],
                "idxml_filter": "run1",
            },
        })
        cfg = load_config(yaml_path, self._make_cli_args(idxml_filter="run2"))
        assert len(cfg.idxml_files) == 1
        assert "run2" in cfg.idxml_files[0]

    def test_filter_zero_results_returns_empty_list(self, tmp_path):
        (tmp_path / "run1.idXML").touch()
        yaml_path = self._make_yaml(tmp_path, {
            "input": {
                "idxml_files": [str(tmp_path / "*.idXML")],
                "idxml_filter": "nomatch",
            },
        })
        cfg = load_config(yaml_path, self._make_cli_args())
        assert cfg.idxml_files == []
