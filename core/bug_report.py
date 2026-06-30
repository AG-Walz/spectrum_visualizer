"""Build pre-filled GitHub issue URLs (bug report / feature request) carrying a
tiny session snapshot: workspace, page, and the last captured error. Settings
live in the workspace's params.json on the server, so the report only points at
the workspace. No Streamlit imports."""
import urllib.parse
from typing import Any

# Conservative ceiling on the pre-filled issue URL. GitHub's edge rejects URLs
# past ~8 KB with HTTP 414, so we stay well under to leave room for encoding.
MAX_PREFILL_URL_BYTES = 7000


def build_snapshot(workspace: str, page: str, error: str | None = None) -> dict[str, Any]:
    """Tiny report payload: workspace (a pointer to params.json on the server),
    current page, and the last error traceback if one was captured."""
    snapshot: dict[str, Any] = {"workspace": workspace, "page": page}
    if error:
        snapshot["error"] = error
    return snapshot


def build_issue_url(repo: str, template: str, snapshot_json: str) -> tuple[str, bool]:
    """Inline the snapshot into the form's `context` textarea, fenced as a
    ```json block. Returns (url, inlined); inlined is False when the snapshot is
    dropped to stay under GitHub's ~8 KB URL limit. The form must live on the
    repo's DEFAULT branch — GitHub resolves `?template=` only there, and `context`
    must match the field id declared in both issue-template YAMLs.
    """
    base = f"https://github.com/{repo}/issues/new"
    fenced = f"```json\n{snapshot_json}\n```"
    url = f"{base}?{urllib.parse.urlencode({'template': template, 'context': fenced})}"
    if len(url.encode("utf-8")) <= MAX_PREFILL_URL_BYTES:
        return url, True
    return f"{base}?{urllib.parse.urlencode({'template': template})}", False
