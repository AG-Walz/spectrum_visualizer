"""Unit tests for core.bug_report — the Streamlit-free issue-URL helpers."""
import urllib.parse
from pathlib import Path

import yaml

import core.bug_report as bug_report
from core.bug_report import build_issue_url, build_snapshot


def _query(url):
    return dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))


def test_build_snapshot_carries_workspace_page_error_and_terminal_log():
    snap = build_snapshot(
        "proca_run", "Load & Annotate",
        error="Traceback...", terminal_log=["ERROR core.x: boom"],
    )
    assert snap == {
        "workspace": "proca_run",
        "page": "Load & Annotate",
        "error": "Traceback...",
        "terminal_log": ["ERROR core.x: boom"],
    }
    bare = build_snapshot("default", "Setup")  # both optional fields omitted when absent
    assert "error" not in bare and "terminal_log" not in bare


def test_build_issue_url_inlines_snapshot_fenced():
    url, complete = build_issue_url("o/r", "bug_report.yml", {"workspace": "ws", "page": "Setup"})
    assert complete is True
    q = _query(url)
    assert q["template"] == "bug_report.yml"
    context = q["context"]
    assert '"workspace": "ws"' in context  # snapshot present, fenced as a json block
    assert context.startswith("```json") and context.endswith("```")


def test_build_issue_url_drops_terminal_log_first_to_keep_workspace_pointer():
    # A huge terminal_log blows the budget; it is shed first so the workspace
    # pointer (the whole point of the report) still inlines.
    snapshot = {"workspace": "ws", "page": "Setup", "terminal_log": ["z" * 9000]}
    url, complete = build_issue_url("o/r", "bug_report.yml", snapshot)
    assert complete is False
    assert len(url.encode("utf-8")) <= bug_report.MAX_PREFILL_URL_BYTES
    context = _query(url)["context"]
    assert '"workspace": "ws"' in context and "terminal_log" not in context


def test_build_issue_url_drops_snapshot_entirely_when_even_pointer_overflows():
    huge = {"workspace": "z" * 9000, "page": "Setup"}  # not droppable, can't fit
    url, complete = build_issue_url("o/r", "bug_report.yml", huge)
    assert complete is False
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
