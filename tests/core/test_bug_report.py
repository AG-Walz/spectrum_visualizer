"""Unit tests for core.bug_report — the Streamlit-free issue-URL helpers."""
import urllib.parse
from pathlib import Path

import yaml

import core.bug_report as bug_report
from core.bug_report import build_issue_url, build_snapshot


def _query(url):
    return dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))


def test_build_snapshot_carries_workspace_page_and_error():
    snap = build_snapshot("proca_run", "Load & Annotate", error="Traceback...")
    assert snap == {
        "workspace": "proca_run",
        "page": "Load & Annotate",
        "error": "Traceback...",
    }
    assert "error" not in build_snapshot("default", "Setup")  # omitted when absent


def test_build_issue_url_inlines_snapshot_fenced():
    url, inlined = build_issue_url("o/r", "bug_report.yml", '{"a": 1}')
    assert inlined is True
    q = _query(url)
    assert q["template"] == "bug_report.yml"
    context = q["context"]
    assert '{"a": 1}' in context  # snapshot present, fenced as a json block
    assert context.startswith("```json") and context.endswith("```")


def test_build_issue_url_drops_oversized_snapshot():
    big = '{"x": "' + "z" * 9000 + '"}'  # blows past the URL byte budget
    url, inlined = build_issue_url("o/r", "bug_report.yml", big)
    assert inlined is False
    assert "context" not in _query(url)
    assert len(url.encode("utf-8")) <= bug_report.MAX_PREFILL_URL_BYTES


_TEMPLATES = Path(__file__).resolve().parents[2] / ".github" / "ISSUE_TEMPLATE"


def _field_ids(name):
    form = yaml.safe_load((_TEMPLATES / name).read_text(encoding="utf-8"))
    return {entry["id"]: entry for entry in form["body"] if "id" in entry}


def test_snapshot_field_exists_and_is_not_rendered():
    # The snapshot lands in `context`; both forms must declare that field and
    # leave it un-rendered, or GitHub silently drops the prefilled value.
    for name in ("bug_report.yml", "feature_request.yml"):
        fields = _field_ids(name)
        assert "context" in fields, f"{name} missing the prefilled `context` field"
        assert "render" not in fields["context"]["attributes"]
