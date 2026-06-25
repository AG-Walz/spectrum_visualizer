"""Bug-report helpers: build a compact, privacy-scrubbed session snapshot and a
pre-filled GitHub issue URL for the AG-Walz bug_report.yml form.

Pure Python, no Streamlit imports, so it stays unit-testable in isolation.
"""
import platform
import sys
import urllib.parse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping

# Packages whose versions are worth recording in a bug report (the Python
# equivalent of R's sessionInfo() that the AG-Walz template asks for).
REPORTED_PACKAGES = (
    "streamlit",
    "openms_insight",
    "pyopenms",
    "psm_utils",
    "polars",
    "pandas",
    "numpy",
    "plotly",
)

# Settings worth capturing in the snapshot — the inputs that shape a run.
# Everything else in session_state (DataFrames, peak arrays, PSM objects) is too
# large or diagnostically irrelevant for a bug report.
SNAPSHOT_SETTING_KEYS = (
    "selected_idxmls",
    "selected_mzmls",
    "target_peptides",
    "reference_substring",
    "ion_types",
    "mz_tolerance",
    "add_neutral_losses",
    "intensity_mode",
    "workers",
    "koina_enabled",
    "koina_model",
    "koina_collision_energy",
    "koina_instrument_type",
    "image_format",
)

# Path-list settings collapse to a count + basenames (drop absolute dirs: smaller
# and avoids leaking server paths); long peptide lists collapse to a sample.
PATH_LIST_KEYS = ("selected_idxmls", "selected_mzmls")
SNAPSHOT_SAMPLE_SIZE = 10
SNAPSHOT_SCHEMA = "mspecviz-bugreport/1"

# Conservative ceiling on the pre-filled issue URL. GitHub's edge rejects URLs
# past ~8 KB with HTTP 414, so we stay well under to leave room for encoding.
MAX_PREFILL_URL_BYTES = 7000


def collect_environment(app_name: str, app_version: str) -> dict[str, Any]:
    """Collect app, Python, OS, and key package versions for the report."""
    packages: dict[str, str] = {}
    for package in REPORTED_PACKAGES:
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = "not installed"
    return {
        "app_name": app_name,
        "app_version": app_version,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": packages,
    }


def format_system_info(environment: Mapping[str, Any]) -> str:
    """Render the environment as the plain-text 'System information' field."""
    lines = [
        f"- {environment.get('app_name', 'MSpecViz')} version: {environment.get('app_version', '?')}",
        f"- Python version: {environment.get('python_version', '?')}",
        f"- Platform: {environment.get('platform', '?')}",
        "- Package versions:",
    ]
    for package, package_version in environment.get("packages", {}).items():
        lines.append(f"    - {package}: {package_version}")
    return "\n".join(lines)


def _summarize_paths(paths: Any) -> Any:
    """Collapse a list of file paths to a count plus basenames (capped)."""
    if not isinstance(paths, (list, tuple)):
        return paths
    names = [Path(path).name for path in paths]
    if len(names) > SNAPSHOT_SAMPLE_SIZE:
        names = names[:SNAPSHOT_SAMPLE_SIZE] + [f"...+{len(paths) - SNAPSHOT_SAMPLE_SIZE} more"]
    return {"count": len(paths), "names": names}


def _summarize_peptides(peptides: Any) -> Any:
    """Collapse target peptides (list or newline-separated string) to count + sample."""
    if isinstance(peptides, str):
        peptides = [line.strip() for line in peptides.splitlines() if line.strip()]
    if not isinstance(peptides, (list, tuple)):
        return peptides
    return {"count": len(peptides), "sample": list(peptides[:SNAPSHOT_SAMPLE_SIZE])}


def build_compact_snapshot(
    environment: Mapping[str, Any],
    settings: Mapping[str, Any],
    status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble a small, inline-friendly snapshot (environment + curated settings).

    Path lists collapse to basenames and long peptide lists to a sample, so the
    JSON stays a few KB. Exotic values are passed through untouched — callers
    serialize with ``json.dumps(..., default=str)``, which stringifies anything
    that is not JSON-native.
    """
    captured: dict[str, Any] = {}
    for key in SNAPSHOT_SETTING_KEYS:
        if key not in settings:
            continue
        if key in PATH_LIST_KEYS:
            captured[key] = _summarize_paths(settings[key])
        elif key == "target_peptides":
            captured[key] = _summarize_peptides(settings[key])
        else:
            captured[key] = settings[key]
    snapshot: dict[str, Any] = {
        "schema": SNAPSHOT_SCHEMA,
        "environment": dict(environment),
        "settings": captured,
    }
    if status:
        snapshot["status"] = dict(status)
    return snapshot


def build_bug_report_url(
    repo: str,
    system_info: str,
    snapshot_json: str | None = None,
    max_url_bytes: int = MAX_PREFILL_URL_BYTES,
) -> tuple[str, bool]:
    """Build a pre-filled URL for the AG-Walz bug_report.yml issue form.

    Query keys match the template's field ids (`command_used`, `files`,
    `system`). The required `description` field is left empty so its template
    placeholder prompts the user (and blocks accidental stub submissions). The
    snapshot is inlined into `command_used` (fenced as a ```json block) when it
    fits; if the URL would risk GitHub's ~8 KB limit it is dropped and the
    returned flag is False. NOTE: the form must live on the repo's DEFAULT
    branch — GitHub resolves `?template=` only there, otherwise it shows a blank
    form and drops every field. `labels` is omitted (the template declares it).
    """
    fields = {
        "template": "bug_report.yml",
        "system": system_info,
        "files": "Drag any relevant files here (e.g. the downloaded mspecviz_bug_report.json).",
    }
    base = f"https://github.com/{repo}/issues/new"
    if snapshot_json:
        fenced = f"```json\n{snapshot_json}\n```"
        url = f"{base}?{urllib.parse.urlencode(dict(fields, command_used=fenced))}"
        if len(url.encode("utf-8")) <= max_url_bytes:
            return url, True
    return f"{base}?{urllib.parse.urlencode(fields)}", False
