"""Regression test for issue #88 point 4: mirror-plot filenames must be unique per PSM.

Before the fix the filename was built from run names only (`{query_run}_vs_{reference_run}`),
so multiple PSMs of the same peptide in the same file pair collapsed onto one PNG and
overwrote each other — leaving fewer PNGs than boxplot points. The scan/spectrum ids in
the name now keep one file per (query PSM, reference PSM) pair.
"""
import numpy as np
from psm_utils import Peptidoform, PSM

from core.peak_annotation import prepare_spectrum_for_plot
from core.pipeline import build_plot_tasks


def _query_psm(spectrum_id: str) -> PSM:
    """A query carrier PSM: same peptidoform + run, distinct spectrum id (i.e. a different scan)."""
    return PSM(
        peptidoform=Peptidoform("PEPTIDEK/2"),
        spectrum_id=spectrum_id,
        run="sample_run",
        retention_time=10.0,
    )


def _frame():
    """A real annotated frame so format_metrics_text has something to score."""
    mz = np.array([227.1026, 276.1554, 500.0, 600.0])
    intensity = np.array([1000.0, 800.0, 200.0, 150.0])
    return prepare_spectrum_for_plot(mz, intensity, "PEPTIDEK", charge=2)


def test_two_psms_same_file_pair_get_distinct_filenames():
    peptidoform = "PEPTIDEK"
    query_a, query_b = _query_psm("scan=1"), _query_psm("scan=2")
    frame = _frame()
    reference_psm = PSM(
        peptidoform=Peptidoform("PEPTIDEK/2"),
        spectrum_id="scan=99", run="reference_run", retention_time=11.0,
    )

    annotated_queries = [(1, query_a, frame), (2, query_b, frame)]
    reference_store = {peptidoform: [(reference_psm, frame)]}
    peptidoform_psms = {peptidoform: [query_a, query_b]}

    tasks = build_plot_tasks(
        annotated_queries, reference_store, peptidoform_psms,
        sanitize_fn=lambda name: name,  # identity so the raw name is visible
        format_dirs={"png": "/tmp/png"},
    )

    filenames = [task.safe_filename for task in tasks]
    assert len(tasks) == 2
    assert len(set(filenames)) == 2, f"filenames collided: {filenames}"
    # Each spectrum id shows up in exactly one file name.
    assert any("scan=1" in name for name in filenames)
    assert any("scan=2" in name for name in filenames)
