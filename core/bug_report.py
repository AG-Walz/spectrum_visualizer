"""Build pre-filled GitHub issue URLs (bug report / feature request) carrying a
tiny session snapshot: workspace, page, the last captured error, and any errors
logged to the terminal this session. Settings live in the workspace's params.json
on the server, so the report only points at the workspace. No Streamlit imports."""
import json
import urllib.parse
from typing import Any

# Conservative ceiling on the pre-filled issue URL. GitHub's edge rejects URLs
# past ~8 KB with HTTP 414, so we stay well under to leave room for encoding.
MAX_PREFILL_URL_BYTES = 7000

# Snapshot fields dropped (largest first) when the inlined URL overflows, so the
# workspace pointer still survives instead of the whole snapshot being discarded.
_DROPPABLE_FIELDS = ("terminal_log", "error")


def build_snapshot(
    workspace: str,
    page: str,
    error: str | None = None,
    terminal_log: list[str] | None = None,
) -> dict[str, Any]:
    """Tiny report payload: workspace (a pointer to params.json on the server),
    current page, the last error traceback if one was captured, and any
    ERROR-level log records seen in the terminal this session."""
    snapshot: dict[str, Any] = {"workspace": workspace, "page": page}
    if error:
        snapshot["error"] = error
    if terminal_log:
        snapshot["terminal_log"] = list(terminal_log)
    return snapshot


def build_issue_url(repo: str, template: str, snapshot: dict[str, Any]) -> tuple[str, bool]:
    """Inline the snapshot into the form's `context` textarea, fenced as a
    ```json block. If the URL would pass GitHub's ~8 KB limit, drop the bulkiest
    optional fields (terminal_log, then error) one at a time so the workspace
    pointer still survives. Returns (url, complete); complete is False when any
    field was dropped. The form must live on the repo's DEFAULT branch — GitHub
    resolves `?template=` only there, and `context` must match the field id
    declared in both issue-template YAMLs.
    """
    base = f"https://github.com/{repo}/issues/new"
    snapshot = dict(snapshot)
    # Try the full snapshot first, then shed the heaviest fields until it fits.
    for dropped, field in enumerate((None, *_DROPPABLE_FIELDS)):
        if field is not None:
            snapshot.pop(field, None)
        if snapshot:
            fenced = f"```json\n{json.dumps(snapshot, indent=2, default=str)}\n```"
            query = {"template": template, "context": fenced}
        else:
            query = {"template": template}
        url = f"{base}?{urllib.parse.urlencode(query)}"
        if len(url.encode("utf-8")) <= MAX_PREFILL_URL_BYTES:
            return url, dropped == 0
    # Even workspace + page overflowed (shouldn't happen): bare template link.
    return f"{base}?{urllib.parse.urlencode({'template': template})}", False
