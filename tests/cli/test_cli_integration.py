"""Integration tests for the CLI pipeline."""
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

    def test_full_pipeline_with_test_data(self, tmp_path):
        source_config = PROJECT_ROOT / "cli" / "config" / "config_example.yaml"
        if not source_config.exists():
            pytest.skip("config_example.yaml not found")

        # Copy the example config but redirect the export directory to tmp_path.
        config_data = yaml.safe_load(source_config.read_text(encoding="utf-8")) or {}
        config_data["export_directory"] = str(tmp_path)
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
