"""Unit tests for core.bug_report — the Streamlit-free bug-report snapshot helpers."""
import json
import urllib.parse
from pathlib import Path

import yaml

import core.bug_report as bug_report
from core.bug_report import (
    SNAPSHOT_SCHEMA,
    build_bug_report_url,
    build_compact_snapshot,
    collect_environment,
    format_system_info,
)


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
    assert "status" not in snapshot  # omitted when absent


def test_compact_snapshot_includes_status_when_present():
    environment = collect_environment("MSpecViz", "1.0.0")
    snapshot = build_compact_snapshot(
        environment,
        settings={"workers": 4},
        status={"data_ready": True, "cache_version": 3},
    )
    assert snapshot["status"]["data_ready"] is True
    assert snapshot["settings"]["workers"] == 4


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
