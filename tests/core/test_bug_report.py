"""Unit tests for core.bug_report — the Streamlit-free bug-report snapshot helpers."""
import json
import urllib.parse
from pathlib import Path

import pandas as pd

import core.bug_report as bug_report
from core.bug_report import (
    MAX_CONTAINER_ITEMS,
    build_bug_report_url,
    build_session_snapshot,
    collect_environment,
    create_snapshot_gist,
    format_system_info,
    sanitize_for_json,
)


def test_sanitize_passes_through_json_natives():
    value = {"a": 1, "b": 1.5, "c": "text", "d": True, "e": None, "f": [1, 2]}
    assert sanitize_for_json(value) == value


def test_sanitize_coerces_paths_and_sets():
    result = sanitize_for_json({"p": Path("/tmp/ws"), "s": {3, 1, 2}})
    assert result == {"p": "/tmp/ws", "s": [1, 2, 3]}


def test_sanitize_summarizes_dataframe_with_shape():
    assert sanitize_for_json(pd.DataFrame({"x": [1, 2, 3]})) == "<DataFrame shape=(3, 1)>"


def test_sanitize_summarizes_large_container():
    big = list(range(MAX_CONTAINER_ITEMS + 1))
    assert sanitize_for_json(big) == f"<list len={MAX_CONTAINER_ITEMS + 1}>"


def test_sanitize_summarizes_opaque_object():
    class Opaque:
        pass

    assert sanitize_for_json(Opaque()) == "<Opaque>"


def test_snapshot_is_json_serializable_with_mixed_state():
    state = {
        "idxml_paths": ["/a/b.idXML"],
        "workspace": Path("/tmp/ws"),
        "mz_tolerance": 0.02,
        "mapping_df": pd.DataFrame({"x": [1, 2, 3]}),
        "nested": {"inner_path": Path("/x"), "values": (1, 2)},
    }
    environment = collect_environment("MSpecViz", "1.0.0")
    snapshot = build_session_snapshot(state, environment)

    dumped = json.loads(json.dumps(snapshot))  # must not raise
    assert dumped["environment"]["app_name"] == "MSpecViz"
    assert dumped["session_state"]["workspace"] == "/tmp/ws"
    assert dumped["session_state"]["mapping_df"] == "<DataFrame shape=(3, 1)>"
    assert dumped["session_state"]["nested"]["inner_path"] == "/x"


def test_collect_environment_reports_known_packages():
    environment = collect_environment("MSpecViz", "1.2.3")
    assert environment["app_version"] == "1.2.3"
    assert "streamlit" in environment["packages"]
    assert environment["python_version"].count(".") == 2


def test_format_system_info_lists_versions():
    environment = {
        "app_name": "MSpecViz",
        "app_version": "1.0.0",
        "python_version": "3.11.14",
        "platform": "Linux",
        "packages": {"streamlit": "1.53.1"},
    }
    text = format_system_info(environment)
    assert "MSpecViz version: 1.0.0" in text
    assert "streamlit: 1.53.1" in text


def _parse_query(url):
    return dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))


def test_build_url_targets_template_form_without_gist():
    url = build_bug_report_url("AG-Walz/spectrum_visualizer", "- sys")
    assert url.startswith("https://github.com/AG-Walz/spectrum_visualizer/issues/new?")
    q = _parse_query(url)
    assert q["template"] == "bug_report.yml"
    assert q["labels"] == "bug"
    assert q["system"] == "- sys"
    assert "drag" in q["files"].lower()  # manual-attach instruction


def test_build_url_embeds_gist_link_when_provided():
    gist = "https://gist.github.com/octocat/abc123"
    q = _parse_query(build_bug_report_url("o/r", "- sys", gist))
    assert gist in q["files"]


def test_create_snapshot_gist_posts_secret_gist(monkeypatch):
    captured = {}

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {"html_url": "https://gist.github.com/octocat/deadbeef"}

    def fake_post(url, headers=None, json=None, timeout=None):
        captured["url"] = url
        captured["headers"] = headers
        captured["json"] = json
        return FakeResponse()

    monkeypatch.setattr(bug_report.requests, "post", fake_post)

    result = create_snapshot_gist('{"a": 1}', token="tok123")

    assert result == "https://gist.github.com/octocat/deadbeef"
    assert captured["url"] == bug_report.GIST_API_URL
    assert captured["headers"]["Authorization"] == "Bearer tok123"
    assert captured["json"]["public"] is False  # secret gist
    assert "mspecviz_bug_report.json" in captured["json"]["files"]


def test_create_snapshot_gist_raises_on_http_error(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            raise RuntimeError("401 Unauthorized")

        def json(self):
            return {}

    monkeypatch.setattr(bug_report.requests, "post", lambda *a, **k: FakeResponse())

    try:
        create_snapshot_gist("{}", token="bad")
    except RuntimeError as exc:
        assert "401" in str(exc)
    else:
        raise AssertionError("expected RuntimeError on HTTP error")
