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
    "openms_insight",
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

# Settings worth capturing in the compact snapshot — the inputs that shape a run.
# Everything else in session_state (DataFrames, peak arrays, PSM objects, widget
# state) is either too large or diagnostically irrelevant.
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
    """Assemble the full downloadable bug-report payload (environment + state)."""
    return {
        "environment": dict(environment),
        "session_state": sanitize_for_json(dict(session_state)),
    }


def _summarize_paths(paths: Any) -> Any:
    """Collapse a list of file paths to a count plus basenames (capped)."""
    if not isinstance(paths, (list, tuple)):
        return sanitize_for_json(paths)
    names = [Path(path).name for path in paths]
    if len(names) > SNAPSHOT_SAMPLE_SIZE:
        names = names[:SNAPSHOT_SAMPLE_SIZE] + [f"...+{len(paths) - SNAPSHOT_SAMPLE_SIZE} more"]
    return {"count": len(paths), "names": names}


def _summarize_peptides(peptides: Any) -> Any:
    """Collapse target peptides (list or newline-separated string) to count + sample."""
    if isinstance(peptides, str):
        peptides = [line.strip() for line in peptides.splitlines() if line.strip()]
    if not isinstance(peptides, (list, tuple)):
        return sanitize_for_json(peptides)
    return {"count": len(peptides), "sample": list(peptides[:SNAPSHOT_SAMPLE_SIZE])}


def build_compact_snapshot(
    environment: Mapping[str, Any],
    settings: Mapping[str, Any],
    status: Mapping[str, Any] | None = None,
    error: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble a small, inline-friendly snapshot (environment + curated settings).

    Path lists collapse to basenames and long peptide lists to a sample, so the
    JSON stays a few KB — small enough to embed directly in the issue URL.
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
            captured[key] = sanitize_for_json(settings[key])
    snapshot: dict[str, Any] = {
        "schema": SNAPSHOT_SCHEMA,
        "environment": dict(environment),
        "settings": captured,
    }
    if status:
        snapshot["status"] = sanitize_for_json(dict(status))
    if error:
        snapshot["error"] = sanitize_for_json(dict(error))
    return snapshot


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
    snapshot_json: str | None = None,
    gist_url: str | None = None,
    max_url_bytes: int = MAX_PREFILL_URL_BYTES,
) -> tuple[str, bool]:
    """Build a pre-filled URL for the AG-Walz bug_report.yml issue form.

    Query keys match the template's field ids (`command_used`, `files`,
    `system`). The required `description` field is left empty on purpose — its
    template `placeholder` prompts the user, and an unfilled required field
    blocks accidental stub submissions. The compact snapshot is inlined into
    `command_used` (fenced as a ```json block so it renders monospace) when it
    fits; if the URL would risk GitHub's ~8 KB limit it is dropped and the
    returned flag is False. NOTE: the form must live on the repo's DEFAULT
    branch — GitHub resolves `?template=` only there, otherwise it shows a blank
    form and silently drops every field. `labels` is omitted (template declares it).
    """
    fields = {
        "template": "bug_report.yml",
        "system": system_info,
    }
    if gist_url:
        fields["files"] = f"Full session snapshot (auto-uploaded): {gist_url}"
    else:
        fields["files"] = "Drag any relevant files here (e.g. the downloaded snapshot)."

    base = f"https://github.com/{repo}/issues/new"
    if snapshot_json:
        fenced = f"```json\n{snapshot_json}\n```"
        with_snapshot = dict(fields, command_used=fenced)
        url = f"{base}?{urllib.parse.urlencode(with_snapshot)}"
        if len(url.encode("utf-8")) <= max_url_bytes:
            return url, True
    return f"{base}?{urllib.parse.urlencode(fields)}", False
