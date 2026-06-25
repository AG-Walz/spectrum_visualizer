"""Build pre-filled GitHub issue URLs (bug report / feature request) carrying a
tiny session snapshot: workspace, page, and the last captured error. Settings
live in the workspace's params.json on the server, so the report only points at
the workspace. No Streamlit imports."""
import urllib.parse
from typing import Any

SNAPSHOT_SCHEMA = "mspecviz-bugreport/1"

# Conservative ceiling on the pre-filled issue URL. GitHub's edge rejects URLs
# past ~8 KB with HTTP 414, so we stay well under to leave room for encoding.
MAX_PREFILL_URL_BYTES = 7000


def build_snapshot(workspace: str, page: str, error: str | None = None) -> dict[str, Any]:
    """Tiny report payload: workspace (a pointer to params.json on the server),
    current page, and the last error traceback if one was captured."""
    snapshot: dict[str, Any] = {"schema": SNAPSHOT_SCHEMA, "workspace": workspace, "page": page}
    if error:
        snapshot["error"] = error
    return snapshot


def build_issue_url(
    repo: str,
    template: str,
    snapshot_json: str,
    field_id: str = "context",
    max_url_bytes: int = MAX_PREFILL_URL_BYTES,
) -> tuple[str, bool]:
    """Inline the snapshot into the form's `field_id` textarea, fenced as a
    ```json block. Returns (url, inlined); inlined is False when the snapshot is
    dropped to stay under GitHub's ~8 KB URL limit. The form must live on the
    repo's DEFAULT branch — GitHub resolves `?template=` only there.
    """
    base = f"https://github.com/{repo}/issues/new"
    fenced = f"```json\n{snapshot_json}\n```"
    url = f"{base}?{urllib.parse.urlencode({'template': template, field_id: fenced})}"
    if len(url.encode("utf-8")) <= max_url_bytes:
        return url, True
    return f"{base}?{urllib.parse.urlencode({'template': template})}", False
