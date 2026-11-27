"""Rendered as the Load & Annotate tab — one page, three sections: the inputs from Setup, a
live preview of the target peptides found in the idXML files, and one "Run analysis" button.

Heavy lifting lives in load_annotate_actions.py; this module handles layout and validation."""
from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_app.utils.job_registry import CANCELLED, DONE, FAILED, get_job
from streamlit_app.utils.state import (
    IDXML_INDEX_DIR, STATE_PREFIX, get_param, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed, show_stale_results_warning
from core.defaults import DEFAULTS
from core.helpers import format_time, parse_peptide_list
from core.idxml_index import ensure_idxml_index, missing_from_index
from core.peptide_lookup import summarize_lookup_results

from streamlit_app.views.load_annotate.load_annotate_actions import (
    add_missing_mzmls,
    analysis_job_key,
    index_idxmls,
    lookup_from_index,
    show_messages,
    start_analysis,
)

params = page_setup("Load & Annotate")

# Display names for the lookup summary columns
_LOOKUP_COLUMNS = {
    "base_sequence": "Peptide", "modification": "Modification",
    "psm_count": "PSMs", "idxml_files": "idXML Files",
    "mzml_resolved": "mzML Resolved", "unresolved": "Unresolved",
}


@st.cache_data(show_spinner="Looking up peptides...", max_entries=4)
def _lookup_preview(index_files: tuple[str, ...], peptides: tuple[str, ...],
                    mzmls: tuple[str, ...]) -> tuple[pd.DataFrame, list[str]]:
    """Per-peptidoform lookup summary for the preview table, plus the run names of the PSMs
    that no selected mzML matches. Cached on its inputs; an index file name changes
    whenever its idXML changes.
    """
    peptidoform_psms = lookup_from_index([Path(path) for path in index_files], list(peptides), list(mzmls))
    unresolved_runs = sorted({
        psm.run for psm_list in peptidoform_psms.values() for psm in psm_list
        if psm.run and not psm.metadata.get("mzml_resolved_path")
    })
    summary = summarize_lookup_results(peptidoform_psms).collect().to_pandas()
    # `sources` (distinct source-folder count) feeds the source-aware idXML counting but
    # isn't surfaced as its own column — drop it from the view.
    summary = summary.drop(columns=["id", "display_label", "sources"], errors="ignore").rename(
        columns=_LOOKUP_COLUMNS)
    return summary, unresolved_runs


def _render_inputs(idxmls: list[str], mzmls: list[str], peptides: list[str]) -> None:
    """Section 1: what the run will use, so the user can check it before pressing Run."""
    st.subheader("1. Inputs")
    idxml_col, mzml_col, peptide_col = st.columns(3)
    idxml_col.metric("idXML files", len(idxmls))
    mzml_col.metric("mzML files", len(mzmls))
    peptide_col.metric("Target peptides", len(peptides))

    # Comparison mode: reference/query split, Koina fallback, or annotation only
    reference_substring = get_param(params, "reference_substring", "")
    if reference_substring:
        needle = reference_substring.lower()
        reference_count = sum(1 for path in mzmls if needle in path.lower())
        st.caption(
            f"Comparison: {reference_count} reference mzML file(s) (path contains "
            f"'{reference_substring}') and {len(mzmls) - reference_count} query file(s). "
            f"Each query spectrum is compared with the reference spectra of its peptide."
        )
    elif get_param(params, "koina_enabled", DEFAULTS.koina_enabled):
        st.caption("Comparison: no reference identifier, so query spectra are compared with "
                   "Koina-predicted spectra.")
    else:
        st.caption("Comparison: none. No reference identifier is set, so spectra are only annotated.")

    ion_types = get_param(params, "ion_types", list(DEFAULTS.ion_types))
    losses = "on" if get_param(params, "add_neutral_losses", DEFAULTS.add_neutral_losses) else "off"
    st.caption(
        f"Annotation: {', '.join(ion_types)} ions, "
        f"{float(get_param(params, 'mz_tolerance', DEFAULTS.mz_tolerance)):.3f} Da tolerance, "
        f"neutral losses {losses}, "
        f"{get_param(params, 'intensity_mode', DEFAULTS.intensity_mode).lower()} intensity. "
        f"Change these on the Setup page."
    )


def _render_lookup_preview(idxmls: list[str], peptides: list[str], mzmls: list[str]) -> bool:
    """Section 2: which target peptides the idXMLs contain. Returns False only when the
    run is known to find nothing (all files indexed and no target peptide present)."""
    st.subheader("2. Peptides found")
    missing = missing_from_index(idxmls, IDXML_INDEX_DIR)
    if missing:
        st.info(
            f"{len(missing)} of {len(idxmls)} idXML file(s) are new to the app. Index them now to "
            f"preview the peptides, or press **Run analysis**, which indexes them first. "
            f"Each file is indexed only once."
        )
        if st.button(f"Index {len(missing)} file(s) now", key=f"{STATE_PREFIX}_index_btn",
                     disabled=_job_running()):
            _index_files, errors = index_idxmls(idxmls, get_param(params, "workers", DEFAULTS.workers))
            if not errors:
                st.rerun()
            show_messages(errors, "idXML warnings")
        return True

    index_files = ensure_idxml_index(idxmls, IDXML_INDEX_DIR)
    summary, unresolved_runs = _lookup_preview(
        tuple(str(path) for path in index_files), tuple(peptides), tuple(mzmls))
    if summary.empty:
        st.warning("None of the target peptides is in the selected idXML files.")
        return False

    st.caption(f"{int(summary['PSMs'].sum())} PSMs across {len(summary)} peptidoform(s).")
    st.dataframe(summary, width="stretch", hide_index=True)
    unresolved = int(summary["Unresolved"].sum())
    if unresolved:
        st.warning(f"{unresolved} PSM(s) from {len(unresolved_runs)} run(s) have no matching file "
                   f"among the selected mzML files and are skipped.")
        if st.button(f"Find the mzML files of these {len(unresolved_runs)} run(s)",
                     key=f"{STATE_PREFIX}_find_mzml_btn", icon=":material/search:",
                     disabled=_job_running()):
            with st.spinner("Searching for mzML files..."):
                st.session_state[get_state_key("mzml_search_report")] = add_missing_mzmls(
                    params, unresolved_runs, idxmls)
            st.rerun()
    # Report of the last search, shown once after the rerun it triggered
    report = st.session_state.pop(get_state_key("mzml_search_report"), None)
    if report:
        st.info(report, icon=":material/search:")
    not_found = sorted(set(peptides) - set(summary["Peptide"]))
    if not_found:
        st.caption(f"Not found in the idXML files ({len(not_found)}): {', '.join(not_found)}")
    return True


def _render_job_status(was_running: bool) -> None:
    """Status of this workspace's background run. Rendered as a fragment that reruns every
    second while the run is active; when the run ends, one full rerun loads the new results.
    was_running tells whether the run was active at the last full rerun.
    """
    job = get_job(analysis_job_key())
    if job is None:
        return
    if job.is_running:
        counts = f": {job.current}/{job.total} {job.units}" if job.total else ""
        st.progress(min(job.current / job.total, 1.0) if job.total else 0.0,
                    text=f"{job.phase}{counts}")
        st.caption(f"{format_time(job.elapsed)} elapsed. The run continues if you switch pages.")
        if job.cancel_requested:
            st.caption("Cancelling after the current step...")
        elif st.button("Cancel run", key=f"{STATE_PREFIX}_cancel_btn", icon=":material/stop:"):
            job.cancel_requested = True
        return
    if was_running:
        st.rerun()  # full rerun: restore_analysis_if_needed picks up the new cache

    if job.status == DONE:
        st.caption(f"Last run: {job.summary}")
    elif job.status == FAILED:
        st.error(f"The last run failed: {job.error}")
        if job.error_traceback:
            st.session_state[get_state_key("last_error")] = job.error_traceback  # for the bug report
    elif job.status == CANCELLED:
        st.info("The last run was cancelled. The previous results (if any) are unchanged.")
    show_messages(job.warnings, "warnings from the last run")


def _job_running() -> bool:
    """True while this workspace's background run is active."""
    job = get_job(analysis_job_key())
    return job is not None and job.is_running


def _render_job_section() -> None:
    """Job status panel; it polls only while a run is active, so an idle page never reruns."""
    running = _job_running()
    st.fragment(run_every=1.0 if running else None)(_render_job_status)(running)


def _render_last_run() -> None:
    """Section 4: the loaded results and an outdated-settings notice."""
    if not st.session_state.get(get_state_key("data_ready")):
        return
    annotated_queries = st.session_state.get(get_state_key("annotated_queries")) or {}
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
    comparison_count = 0 if comparison_frame is None else len(comparison_frame)
    st.success(
        f"Results loaded: {len(annotated_queries)} spectra and {comparison_count} "
        f"query/reference comparisons. Open the **Explore** page to view them."
    )
    show_stale_results_warning(params, "Press **Run analysis** to update them.")


def render():
    """Page entry point: inputs, peptide preview, Run button, then the loaded results."""
    restore_analysis_if_needed()
    st.title("Load & Annotate")
    st.caption("Find the target peptides in the idXML files, then load, annotate and compare "
               "their spectra.")

    # Guard: every input comes from the Setup page
    idxmls = get_param(params, "selected_idxmls", [])
    mzmls = get_param(params, "selected_mzmls", [])
    peptides = parse_peptide_list(get_param(params, "target_peptides", ""))
    # mzML files are not required to preview: the preview can find them (see below)
    missing_inputs = [label for label, value in
                      (("idXML files", idxmls), ("target peptides", peptides))
                      if not value]
    if missing_inputs:
        st.warning(f"Missing on the Setup page: {', '.join(missing_inputs)}.")
        _render_job_section()
        _render_last_run()
        return

    _render_inputs(idxmls, mzmls, peptides)
    can_run = _render_lookup_preview(idxmls, peptides, mzmls)

    st.subheader("3. Run")
    if not mzmls:
        st.caption("Select mzML files on the Setup page, or let the button above find them.")
    if st.button("Run analysis", key=f"{STATE_PREFIX}_run_btn", type="primary",
                 icon=":material/play_arrow:", disabled=not can_run or not mzmls or _job_running()):
        if start_analysis(params, idxmls, peptides):
            st.rerun()  # show the live status panel at once
        st.warning("A run is already in progress for this workspace.")
    _render_job_section()
    _render_last_run()

    save_params(params)


render()
