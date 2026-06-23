"""Bug-report helpers: sanitize the Streamlit session state into JSON and
collect system information for the AG-Walz bug_report.yml issue form.

Pure Python, no Streamlit imports, so it stays unit-testable in isolation.
"""
import json
import platform
import sys
import urllib.parse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping

import requests

# GitHub's documented limits / endpoints, kept here so callers don't hard-code them.
GIST_API_URL = "https://api.github.com/gists"
GIST_TIMEOUT_SECONDS = 15

# Packages whose versions are worth recording in a bug report (the Python
# equivalent of R's sessionInfo() that the AG-Walz template asks for).
REPORTED_PACKAGES = (
    "streamlit",
    "pyopenms",
    "psm_utils",
    "polars",
    "pandas",
    "numpy",
    "plotly",
)

# Lists/tuples longer than this are summarized rather than expanded element by
# element — a bug report does not need thousands of "<PSM>" placeholders.
MAX_CONTAINER_ITEMS = 100


def sanitize_for_json(value: Any) -> Any:
    """Return a JSON-serializable view of value, summarizing anything exotic.

    Containers are walked recursively; unserializable objects (DataFrames,
    pyopenms/psm_utils objects, ...) become a short ``<type ...>`` summary
    string rather than raising.
    """
    # Fast path: json-native scalars.
    if value is None or isinstance(value, (bool, int, float, str)):
        return value

    # Paths serialize as their string form.
    if isinstance(value, Path):
        return str(value)

    # Mappings: sanitize each value, coercing keys to strings.
    if isinstance(value, Mapping):
        return {str(key): sanitize_for_json(item) for key, item in value.items()}

    # Ordered containers (summarize very large ones).
    if isinstance(value, (list, tuple)):
        if len(value) > MAX_CONTAINER_ITEMS:
            return _summarize(value)
        return [sanitize_for_json(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted(sanitize_for_json(item) for item in value)

    # numpy scalars expose .item() without importing numpy here.
    item_method = getattr(value, "item", None)
    if callable(item_method) and not hasattr(value, "__len__"):
        try:
            return sanitize_for_json(item_method())
        except Exception:
            pass

    # Anything else: try a direct JSON round-trip, else summarize it.
    try:
        json.dumps(value)
        return value
    except (TypeError, ValueError):
        return _summarize(value)


def _summarize(value: Any) -> str:
    """Build a ``<ClassName ...>`` summary for an unserializable object."""
    type_name = type(value).__name__
    shape = getattr(value, "shape", None)
    if shape is not None:
        return f"<{type_name} shape={tuple(shape)}>"
    length = None
    try:
        length = len(value)
    except (TypeError, AttributeError):
        pass
    if length is not None:
        return f"<{type_name} len={length}>"
    return f"<{type_name}>"


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


def build_session_snapshot(
    session_state: Mapping[str, Any],
    environment: Mapping[str, Any],
) -> dict[str, Any]:
    """Assemble the downloadable bug-report payload (environment + state)."""
    return {
        "environment": dict(environment),
        "session_state": sanitize_for_json(dict(session_state)),
    }


def create_snapshot_gist(
    report_json: str,
    token: str,
    description: str = "MSpecViz session snapshot",
) -> str:
    """POST the snapshot as a SECRET gist and return its html_url. Raises on failure.

    GitHub has no API to attach a file to an issue, so the report is hosted as a
    gist and linked from the issue instead. Needs a token with only the gist scope.
    """
    response = requests.post(
        GIST_API_URL,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.github+json",
        },
        json={
            "description": description,
            "public": False,
            "files": {"mspecviz_bug_report.json": {"content": report_json}},
        },
        timeout=GIST_TIMEOUT_SECONDS,
    )
    response.raise_for_status()
    return response.json()["html_url"]


def build_bug_report_url(
    repo: str,
    system_info: str,
    gist_url: str | None = None,
) -> str:
    """Build a pre-filled GitHub URL for the AG-Walz bug_report.yml issue form.

    Query keys match the template's field ids (`description`, `system`, `files`).
    When gist_url is given, the `files` field links the auto-uploaded snapshot;
    otherwise it instructs the reporter to drag in the downloaded file.
    """
    if gist_url:
        files = f"Full session snapshot (auto-uploaded): {gist_url}"
    else:
        files = "Attach the downloaded mspecviz_bug_report.json here (drag & drop)."
    query = urllib.parse.urlencode(
        {
            "template": "bug_report.yml",
            "labels": "bug",
            "description": "Describe the bug here.\n",
            "system": system_info,
            "files": files,
        }
    )
    return f"https://github.com/{repo}/issues/new?{query}"
