"""Unit tests for core.bug_report — the Streamlit-free bug-report snapshot helpers."""
import json
import urllib.parse
from pathlib import Path

import pandas as pd
import yaml

import core.bug_report as bug_report
from core.bug_report import (
    MAX_CONTAINER_ITEMS,
    SNAPSHOT_SCHEMA,
    build_bug_report_url,
    build_compact_snapshot,
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


def test_build_url_targets_template_form_without_snapshot():
    url, inlined = build_bug_report_url("AG-Walz/spectrum_visualizer", "- sys")
    assert inlined is False  # nothing to inline
    assert url.startswith("https://github.com/AG-Walz/spectrum_visualizer/issues/new?")
    q = _parse_query(url)
    assert q["template"] == "bug_report.yml"
    assert "labels" not in q  # template already declares the label
    assert q["system"] == "- sys"
    # description is deliberately left empty: the required field's placeholder
    # prompts the user, so a stub value can't be submitted by accident.
    assert "description" not in q
    assert "drag" in q["files"].lower()  # manual-attach instruction


def test_build_url_inlines_small_snapshot_into_command_used():
    url, inlined = build_bug_report_url("o/r", "- sys", snapshot_json='{"a": 1}')
    assert inlined is True
    command_used = _parse_query(url)["command_used"]
    assert '{"a": 1}' in command_used  # snapshot is present...
    assert command_used.startswith("```json")  # ...fenced so it renders monospace
    assert command_used.endswith("```")


def test_build_url_drops_snapshot_when_over_budget():
    big = '{"x": "' + "z" * 9000 + '"}'  # blows past the URL byte budget
    url, inlined = build_bug_report_url("o/r", "- sys", snapshot_json=big)
    assert inlined is False
    assert "command_used" not in _parse_query(url)
    assert len(url.encode("utf-8")) <= bug_report.MAX_PREFILL_URL_BYTES


def test_build_url_embeds_gist_link_when_provided():
    gist = "https://gist.github.com/octocat/abc123"
    url, inlined = build_bug_report_url("o/r", "- sys", gist_url=gist)
    assert inlined is False
    assert gist in _parse_query(url)["files"]


def test_compact_snapshot_collapses_paths_and_peptides():
    environment = collect_environment("MSpecViz", "1.0.0")
    settings = {
        "selected_idxmls": ["/data/runs/a.idXML", "/data/runs/b.idXML"],
        "selected_mzmls": [f"/data/{i}.mzML" for i in range(15)],
        "target_peptides": "PEPTIDEK\nELVISLIVES\n",
        "mz_tolerance": 0.02,
        "ion_types": ["b", "y"],
    }
    snapshot = build_compact_snapshot(environment, settings)

    assert snapshot["schema"] == SNAPSHOT_SCHEMA
    assert snapshot["settings"]["selected_idxmls"] == {
        "count": 2,
        "names": ["a.idXML", "b.idXML"],  # basenames only, no leaked dirs
    }
    mzmls = snapshot["settings"]["selected_mzmls"]
    assert mzmls["count"] == 15 and "...+5 more" in mzmls["names"][-1]  # capped sample
    assert snapshot["settings"]["target_peptides"] == {
        "count": 2,
        "sample": ["PEPTIDEK", "ELVISLIVES"],
    }
    assert snapshot["settings"]["mz_tolerance"] == 0.02
    assert "status" not in snapshot and "error" not in snapshot  # omitted when absent


def test_compact_snapshot_includes_status_and_error_when_present():
    environment = collect_environment("MSpecViz", "1.0.0")
    snapshot = build_compact_snapshot(
        environment,
        settings={"workers": 4},
        status={"data_ready": True, "cache_version": 3},
        error={"type": "KeyError", "message": "peptide_mirror.psm_lf"},
    )
    assert snapshot["status"]["data_ready"] is True
    assert snapshot["error"]["type"] == "KeyError"


def test_compact_snapshot_is_small_and_json_serializable():
    environment = collect_environment("MSpecViz", "1.0.0")
    settings = {"selected_idxmls": [f"/data/{i}.idXML" for i in range(50)]}
    snapshot = build_compact_snapshot(environment, settings)
    dumped = json.dumps(snapshot)  # must not raise
    assert len(dumped.encode("utf-8")) < bug_report.MAX_PREFILL_URL_BYTES


# The "issue opens empty" class of bug is a field-id <-> query-param mismatch,
# or a snapshot routed into a field GitHub refuses to prefill. These two tests
# pin the contract between build_bug_report_url and the actual issue-form YAML.
_BUG_REPORT_TEMPLATE = (
    Path(__file__).resolve().parents[2] / ".github" / "ISSUE_TEMPLATE" / "bug_report.yml"
)


def _template_fields() -> dict[str, dict]:
    """Map each issue-form field id -> its attributes block from bug_report.yml."""
    form = yaml.safe_load(_BUG_REPORT_TEMPLATE.read_text(encoding="utf-8"))
    return {entry["id"]: entry for entry in form["body"] if "id" in entry}


def test_prefill_params_match_template_field_ids():
    fields = _template_fields()
    url, _ = build_bug_report_url("o/r", "- sys", snapshot_json='{"a": 1}')
    params = _parse_query(url)
    for key in params:
        if key == "template":  # not a field id; selects the form
            continue
        assert key in fields, f"prefill param {key!r} has no matching field id in bug_report.yml"


def test_snapshot_target_field_is_not_rendered():
    # GitHub will not reliably prefill a textarea that declares `render:`; the
    # snapshot must land in a plain textarea or it silently disappears.
    command_used = _template_fields()["command_used"]
    assert "render" not in command_used["attributes"]


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
