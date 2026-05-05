"""
Peptide input widget — paste-list and pick-from-idXML tabs.

The "pick from idXML" flow scans selected idXML files in parallel for PSM counts
so the picker shows how many hits each peptide has. Uses
``mspecviz.common.parallel.run_tasks`` so the sequential/parallel dispatch is
consistent with the rest of the codebase.
"""
import time as _time
from collections import Counter
from pathlib import Path

import streamlit as st

from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.parallel import run_tasks
from mspecviz.common.text_utils import parse_peptide_list
from mspecviz.search import count_psms_per_peptide
from streamlit_app.utils.state import STATE_PREFIX, get_state_key


def _scan_idxmls_for_psm_counts(idxmls: list[str], workers: int) -> tuple[Counter, list[str]]:
    """Run count_psms_per_peptide across files in a st.status block."""
    merged: Counter = Counter()
    errors: list[str] = []
    n = len(idxmls)

    with st.status(f"Scanning {n} idXML files ({workers} workers)...", expanded=True) as status:
        progress = st.progress(0)
        info = st.empty()
        t0 = _time.time()

        def _on_result(task, result):
            merged.update(result)

        def _on_error(task, exc):
            errors.append(f"{Path(task).name}: {exc}")

        def _on_progress(done, total):
            elapsed = _time.time() - t0
            remaining = (elapsed / done) * (total - done) if done < total else 0
            progress.progress(
                done / total,
                f"{done}/{total} files — {elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining",
            )

        def _on_start(task):
            info.caption(f"Scanning {Path(task).name}")

        run_tasks(
            tasks=idxmls,  # single-arg: auto-wrapped by run_tasks
            worker=count_psms_per_peptide,
            workers=min(workers, n),
            on_start=_on_start,
            on_result=_on_result,
            on_error=_on_error,
            on_progress=_on_progress,
        )

        progress.progress(1.0, f"Done — {len(merged)} peptides from {n} files in {_time.time()-t0:.1f}s")
        if errors:
            status.update(label=f"Completed with {len(errors)} errors", state="error")
            for e in errors[:5]:
                st.caption(e)
        else:
            status.update(label=f"{len(merged)} peptides from {n} files", state="complete")

    return merged, errors


def _render_paste_tab(initial_value: str, persisted_key: str) -> None:
    """Paste-list tab: textarea + live count."""
    col1, col2 = st.columns([2, 1])
    with col1:
        peptide_text = st.text_area(
            "Enter peptide sequences (one per line)",
            value=initial_value, height=150, key="target_peptides",
            placeholder="PEPTIDESEQ\nANOTHERPEPTIDE\n...",
        )
    # Persist in a non-widget key that survives page navigation
    st.session_state[persisted_key] = peptide_text
    with col2:
        st.markdown("**Tips:**\n- One peptide per line\n"
                    "- Without modifications\n"
                    "- Case insensitive")
        peptides = parse_peptide_list(peptide_text)
        if peptides:
            st.success(f"{len(peptides)} peptide(s) entered")


def _ensure_psm_counts_for(idxmls: list[str], workers: int) -> dict:
    """Return cached PSM counts for the given idXML files, scanning if needed."""
    scan_key = '_idxml_psm_counts'
    scan_files_key = '_idxml_scan_files'
    psm_counts = st.session_state.get(scan_key, {})
    scanned_files = st.session_state.get(scan_files_key)

    if psm_counts and scanned_files == idxmls:
        return psm_counts

    n = len(idxmls)
    st.info(f"{n} idXML file{'s' if n != 1 else ''} selected")
    if st.button(f"Scan {n} files for peptides",
                 key=f"{STATE_PREFIX}_scan_idxmls", type="secondary"):
        merged, _errors = _scan_idxmls_for_psm_counts(idxmls, workers)
        st.session_state[scan_key] = dict(sorted(merged.items()))
        st.session_state[scan_files_key] = idxmls
        st.rerun()
    return psm_counts


def _render_pick_filters(psm_counts: dict) -> tuple[str, int]:
    """Render the filter input + min-PSM number_input; returns ``(text, min_psm)``."""
    f_col1, f_col2 = st.columns([2, 1])
    with f_col1:
        filter_text = st.text_input(
            "Filter", placeholder="Type to filter by substring...",
            key=f"{STATE_PREFIX}_pick_filter", label_visibility="collapsed",
        )
    with f_col2:
        max_count = max(psm_counts.values())
        min_psm = st.number_input(
            "Min PSMs", min_value=1, max_value=max_count, value=1, step=1,
            key=f"{STATE_PREFIX}_pick_min_psm",
            help="Show only peptides with at least this many PSMs",
        )
    return filter_text, min_psm


def _render_pick_multiselect(psm_counts: dict, filter_text: str, min_psm: int) -> list[str]:
    """Render the peptide multiselect; returns the list of selected sequences."""
    shown_items = [
        (seq, cnt) for seq, cnt in psm_counts.items()
        if cnt >= min_psm and (not filter_text or filter_text.upper() in seq)
    ]
    labels = [f"{seq} ({cnt})" for seq, cnt in shown_items]
    label_to_seq = dict(zip(labels, (seq for seq, _ in shown_items)))
    current = set(parse_peptide_list(st.session_state.get('target_peptides', '')))
    defaults = [lbl for lbl, seq in label_to_seq.items() if seq in current]
    selected_labels = st.multiselect(
        "Select peptides", options=labels, default=defaults,
        key=f"{STATE_PREFIX}_pick_multi",
    )
    return [label_to_seq[lbl] for lbl in selected_labels]


def _commit_picked_peptides(selected: list[str], pending_key: str, mode: str) -> None:
    """Write *selected* into ``pending_key`` either replacing or appending."""
    if mode == "append":
        existing = parse_peptide_list(st.session_state.get('target_peptides', ''))
        peps = list(dict.fromkeys(existing + selected))
    else:
        peps = selected
    st.session_state[pending_key] = "\n".join(peps)
    st.rerun()


def _render_pick_actions(selected: list[str], pending_key: str) -> None:
    """Render the 'Set as target' / 'Add to existing' buttons."""
    c1, c2 = st.columns(2)
    with c1:
        if st.button("Set as target peptides", key=f"{STATE_PREFIX}_pick_set",
                     type="primary"):
            _commit_picked_peptides(selected, pending_key, "replace")
    with c2:
        if st.button("Add to existing", key=f"{STATE_PREFIX}_pick_add"):
            _commit_picked_peptides(selected, pending_key, "append")


def _render_pick_tab(idxmls: list[str], workers: int, pending_key: str) -> None:
    """Pick-from-idXML tab: scans files, renders filterable multiselect."""
    if not idxmls:
        st.info("Select idXML files above first.")
        return

    psm_counts = _ensure_psm_counts_for(idxmls, workers)
    if not psm_counts:
        return

    st.caption(f"{len(psm_counts)} unique peptides found")
    filter_text, min_psm = _render_pick_filters(psm_counts)
    selected = _render_pick_multiselect(psm_counts, filter_text, min_psm)
    _render_pick_actions(selected, pending_key)


def render_peptide_input(params: dict) -> None:
    """Render the Target Peptides widget with paste + pick tabs."""
    pending_key = "_pending_target_peptides"
    persisted_key = get_state_key('target_peptides_persisted')
    if pending_key in st.session_state:
        st.session_state['target_peptides'] = st.session_state.pop(pending_key)

    st.subheader("Target Peptides")

    idxmls = st.session_state.get('selected_idxmls', params.get('selected_idxmls', []))

    # Initial value: persisted session_state (survives page nav) > workspace params
    initial_value = st.session_state.get(persisted_key, params.get('target_peptides', ''))

    paste_tab, pick_tab = st.tabs(["Paste list", "Pick from idXML"])
    with paste_tab:
        _render_paste_tab(initial_value, persisted_key)
    with pick_tab:
        _render_pick_tab(idxmls, params.get('workers', DEFAULTS.workers), pending_key)
