"""Tests for the per-peptide plot axis and the retention times carried into exported plots/tables.

Issue #114: every plot of one peptidoform must span the same 0 -> max(m/z of the group) x-axis,
so the exported PNGs of a peptide can be flipped through without the scale shifting underneath.
Issue #113: each task carries both spectra's RTs so the worker's summary row reports them.
"""
import math

import numpy as np
from psm_utils import Peptidoform, PSM

from core.peak_annotation import prepare_spectrum_for_plot
from core.pipeline import build_plot_tasks, peptidoform_mz_axis_maxima
from core.plotting.mirror import _MZ_PADDING_DA
from core.plotting.plot_worker import _build_summary_row


def _psm(spectrum_id: str, peptidoform: str = "PEPTIDEK/2", run: str = "sample_run",
         retention_time: float = 10.0) -> PSM:
    return PSM(peptidoform=Peptidoform(peptidoform), spectrum_id=spectrum_id,
               run=run, retention_time=retention_time)


def _frame(top_mz: float, sequence: str = "PEPTIDEK"):
    """An annotated frame whose highest peak sits at top_mz."""
    mz = np.array([227.1026, 276.1554, top_mz])
    intensity = np.array([1000.0, 800.0, 200.0])
    return prepare_spectrum_for_plot(mz, intensity, sequence, charge=2)


class TestPeptidoformAxisMaxima:
    def test_bound_covers_queries_and_references_of_the_group(self):
        short_query, long_query = _psm("scan=1"), _psm("scan=2")
        annotated_queries = [(1, short_query, _frame(500.0)), (2, long_query, _frame(1200.0))]
        reference_store = {"PEPTIDEK": [(_psm("scan=99", run="reference_run"), _frame(900.0))]}

        maxima = peptidoform_mz_axis_maxima(
            annotated_queries, reference_store, {1: "PEPTIDEK", 2: "PEPTIDEK"})

        assert maxima["PEPTIDEK"] == 1200.0 + _MZ_PADDING_DA

    def test_reference_can_set_the_bound_alone(self):
        # A reference reaching further than any query still widens the whole group's axis.
        annotated_queries = [(1, _psm("scan=1"), _frame(500.0))]
        reference_store = {"PEPTIDEK": [(_psm("scan=99", run="reference_run"), _frame(1400.0))]}

        maxima = peptidoform_mz_axis_maxima(annotated_queries, reference_store, {1: "PEPTIDEK"})

        assert maxima["PEPTIDEK"] == 1400.0 + _MZ_PADDING_DA

    def test_groups_are_scoped_per_peptidoform(self):
        annotated_queries = [
            (1, _psm("scan=1"), _frame(500.0)),
            (2, _psm("scan=2", peptidoform="PEPTIDEKK/2"), _frame(1200.0, "PEPTIDEKK")),
        ]
        maxima = peptidoform_mz_axis_maxima(
            annotated_queries, {}, {1: "PEPTIDEK", 2: "PEPTIDEKK"})

        assert maxima["PEPTIDEK"] == 500.0 + _MZ_PADDING_DA
        assert maxima["PEPTIDEKK"] == 1200.0 + _MZ_PADDING_DA


class TestExportedTasksShareTheGroupAxis:
    def _tasks(self, **kwargs):
        short_query, long_query = _psm("scan=1"), _psm("scan=2")
        reference_psm = _psm("scan=99", run="reference_run", retention_time=11.0)
        return build_plot_tasks(
            [(1, short_query, _frame(500.0)), (2, long_query, _frame(1200.0))],
            {"PEPTIDEK": [(reference_psm, _frame(900.0))]},
            {"PEPTIDEK": [short_query, long_query]},
            sanitize_fn=lambda name: name,
            format_dirs={"png": "/tmp/png"},
            **kwargs,
        )

    def test_all_tasks_of_a_peptide_carry_one_bound(self):
        bounds = {task.mz_axis_max for task in self._tasks()}
        assert bounds == {1200.0 + _MZ_PADDING_DA}

    def test_unpaired_single_tasks_get_the_bound_too(self):
        query = _psm("scan=1")
        tasks = build_plot_tasks(
            [(1, query, _frame(500.0))], {}, {"PEPTIDEK": [query]},
            sanitize_fn=lambda name: name, format_dirs={"png": "/tmp/png"},
            plot_unpaired_queries=True,
        )
        assert [task.mz_axis_max for task in tasks] == [500.0 + _MZ_PADDING_DA]

    def test_tasks_carry_both_retention_times(self):
        task = self._tasks()[0]
        assert task.query_rt == 10.0
        assert task.reference_rt == 11.0


class TestWorkerSummaryRowRts:
    """The export worker's summary row reports the RTs it was handed (#113)."""

    def _task(self, **kwargs):
        tasks = build_plot_tasks(
            [(1, _psm("scan=1"), _frame(500.0))],
            {"PEPTIDEK": [(_psm("scan=99", run="reference_run", **kwargs), _frame(900.0))]},
            {"PEPTIDEK": [_psm("scan=1")]},
            sanitize_fn=lambda name: name, format_dirs={"png": "/tmp/png"},
        )
        return tasks[0]

    def test_row_reports_rts_and_their_difference(self):
        row = _build_summary_row(self._task(retention_time=8.0), {"Cos": 0.9})
        assert row["query_rt"] == 10.0
        assert row["ref_rt"] == 8.0
        assert row["rt_delta"] == 2.0

    def test_predicted_reference_leaves_the_delta_blank(self):
        row = _build_summary_row(self._task(retention_time=0.0), {"Cos": 0.9})
        assert math.isnan(row["rt_delta"])
