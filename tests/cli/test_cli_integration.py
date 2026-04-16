"""Integration tests for the CLI pipeline."""
import pytest
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLI_ARGS = [sys.executable, "-m", "cli"]


@pytest.mark.integration
class TestCLI:
    def test_help_exits_zero(self):
        result = subprocess.run(CLI_ARGS + ["--help"], cwd=PROJECT_ROOT,
                                capture_output=True, text=True, timeout=30,
                                encoding="utf-8", errors="replace")
        assert result.returncode == 0
        assert "MSpecViz CLI" in result.stdout

    def test_missing_config_exits_two(self):
        result = subprocess.run(CLI_ARGS + ["-c", "/nonexistent.yaml"], cwd=PROJECT_ROOT,
                                capture_output=True, text=True, timeout=30,
                                encoding="utf-8", errors="replace")
        assert result.returncode == 2

    def test_full_pipeline_with_test_data(self, tmp_path):
        config = PROJECT_ROOT / "cli" / "config" / "config_example.yaml"
        if not config.exists():
            pytest.skip("config_example.yaml not found")

        result = subprocess.run(
            CLI_ARGS + ["-c", str(config), "--export-directory", str(tmp_path), "-q"],
            cwd=PROJECT_ROOT,
            capture_output=True, text=True, timeout=120,
            encoding="utf-8", errors="replace",
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Check outputs were created
        export_dirs = list(tmp_path.iterdir())
        assert len(export_dirs) >= 1
