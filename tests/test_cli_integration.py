"""Integration tests for the CLI pipeline."""
import pytest
import subprocess
import sys
from pathlib import Path

CLI_PATH = Path(__file__).parent.parent / "cli.py"
PYTHON = "/mnt/volume/bin/miniconda/envs/msinsight/bin/python"


@pytest.mark.integration
class TestCLI:
    def test_help_exits_zero(self):
        result = subprocess.run([PYTHON, str(CLI_PATH), "--help"],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0
        assert "MSInsight CLI" in result.stdout

    def test_missing_config_exits_two(self):
        result = subprocess.run([PYTHON, str(CLI_PATH), "-c", "/nonexistent.yaml"],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 2

    def test_full_pipeline_with_test_data(self, tmp_path):
        config = Path(__file__).parent.parent / "config" / "config_example.yaml"
        if not config.exists():
            pytest.skip("config_example.yaml not found")

        result = subprocess.run(
            [PYTHON, str(CLI_PATH), "-c", str(config),
             "--export-dir", str(tmp_path), "-q"],
            capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Check outputs were created
        export_dirs = list(tmp_path.iterdir())
        assert len(export_dirs) >= 1
