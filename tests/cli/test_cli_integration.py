"""Integration tests for cli.main — subprocess-driven against the small_dataset fixture."""
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLI_ARGS = [sys.executable, "-m", "cli"]


@pytest.mark.integration
@pytest.mark.slow
class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run(CLI_ARGS, cwd=PROJECT_ROOT,
                                capture_output=True, text=True, timeout=30,
                                encoding="utf-8", errors="replace")
        assert result.returncode != 0
        assert "Usage:" in (result.stderr + result.stdout)

    def test_missing_config_exits_two(self):
        result = subprocess.run(CLI_ARGS + ["/nonexistent.yaml"], cwd=PROJECT_ROOT,
                                capture_output=True, text=True, timeout=30,
                                encoding="utf-8", errors="replace")
        assert result.returncode != 0

    def test_full_pipeline_with_test_data(self, tmp_path, idxml_path, mzml_paths, target_peptides):
        # Build the config from the in-repo examples/ fixtures rather than cloning
        # config_example.yaml (whose paths point at developer-specific data). All
        # other keys fall back to DEFAULTS. export_directory -> tmp_path so the
        # test never writes into the developer's real exports folder.
        config_data = {
            "idxml_files": [idxml_path],
            "mzml_files": list(mzml_paths),
            "peptide_source": list(target_peptides),
            "reference_substring": "",
            "export_formats": ["png"],
            "export_directory": str(tmp_path),
        }
        config_path = tmp_path / "test_config.yaml"
        config_path.write_text(yaml.safe_dump(config_data), encoding="utf-8")

        env = {**os.environ, "MSPECVIZ_LOG_LEVEL": "ERROR"}
        result = subprocess.run(
            CLI_ARGS + [str(config_path)], cwd=PROJECT_ROOT,
            capture_output=True, text=True, timeout=120,
            encoding="utf-8", errors="replace", env=env,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        export_dirs = [p for p in tmp_path.iterdir() if p.is_dir()]
        assert len(export_dirs) >= 1
