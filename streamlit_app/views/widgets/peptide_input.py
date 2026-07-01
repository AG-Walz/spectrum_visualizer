"""Peptide-input widget for the Setup page: paste-list tab + pick-from-idXML tab (PSM counts
via core.peptide_lookup.load_idxmls + count_psms_per_peptide; the loaded PSMList is cached for Load & Annotate)."""
from collections import Counter

import streamlit as st

from core.defaults import DEFAULTS
from core.helpers import parse_peptide_list
from core.peptide_lookup import count_psms_per_peptide, load_idxmls
from streamlit_app.utils.progress_tracker import progress_status
from streamlit_app.utils.state import (
    STATE_PREFIX, cache_idxml_psmlist, get_state_key,
)


def _scan_idxmls_for_psm_counts(idxmls: list[str], workers: int) -> tuple[Counter, list[str]]:
    """Load idxmls once into a PSMList (cached for reuse by Load & Annotate) with a live
    progress bar, then count PSMs per peptide from memory. Returns (counts, errors).
    """
    # Setup: error accumulator + progress-status context
    errors: list[str] = []
    n = len(idxmls)
    with progress_status(f"Scanning {n} idXML files ({workers} workers)...") as tracker:
        # Parse every idXML once; cache the PSMList so target lookup can reuse it
        psm_list = load_idxmls(
            idxmls,
            workers=min(workers, n),
            progress_callback=lambda done, total, _msg: tracker.update(done, total, units="files"),
            error_callback=errors.append,
        )
        cache_idxml_psmlist(idxmls, psm_list)
        # Derive peptide counts from the loaded PSMList (no file re-read)
        psm_counts = count_psms_per_peptide(psm_list)
        # Finalise: surface error summary or success count
        if errors:
            tracker.finalize(f"Completed with {len(errors)} errors", state="error")
            for error in errors[:5]:
                st.caption(error)
        else:
            tracker.finalize(f"{len(psm_counts)} peptides from {n} files")
    return psm_counts, errors


def _render_paste_tab(initial_value: str, persisted_key: str) -> None:
    """Render the paste-list tab: textarea + live peptide count."""
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
            # Showcase the actual sequences so a loaded/precomputed workspace reveals
            # *which* peptides are in play, not just the count (issue #90).
            with st.expander("Show peptides"):
                st.code("\n".join(peptides))


def _ensure_psm_counts_for(idxmls: list[str], workers: int) -> dict:
    """Return cached PSM counts for idxmls; render a Scan button when no cache exists."""
    # Cache lookup: keys for the counts dict and the list of files last scanned
    scan_key = '_idxml_psm_counts'
    scan_files_key = '_idxml_scan_files'
    psm_counts = st.session_state.get(scan_key, {})
    scanned_files = st.session_state.get(scan_files_key)

    # Fast path: return cached counts when the file list is unchanged
    if psm_counts and scanned_files == idxmls:
        return psm_counts

    # Render scan trigger: button kicks off the parallel scan + caches results
    n = len(idxmls)
    st.info(f"{n} idXML file{'s' if n != 1 else ''} selected")
    if st.button(f"Scan {n} files for peptides",
                 key=f"{STATE_PREFIX}_scan_idxmls", type="secondary"):
        scanned_psm_counts, _errors = _scan_idxmls_for_psm_counts(idxmls, workers)
        st.session_state[scan_key] = dict(sorted(scanned_psm_counts.items()))
        st.session_state[scan_files_key] = idxmls
        st.rerun()
    return psm_counts


def _render_pick_filters(psm_counts: dict) -> tuple[str, int]:
    """Render the substring filter + min-PSM threshold; return (filter_text, min_psm)."""
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
    """Render the peptide multiselect and return the bare sequences the user picked."""
    # Filter rows: respect min-PSM threshold and the optional substring filter
    shown_items = [
        (sequence, count) for sequence, count in psm_counts.items()
        if count >= min_psm and (not filter_text or filter_text.upper() in sequence)
    ]
    # Build label↔sequence mapping and pre-select sequences already in the textarea
    labels = [f"{sequence} ({count})" for sequence, count in shown_items]
    label_to_sequence = dict(zip(labels, (sequence for sequence, _ in shown_items)))
    current = set(parse_peptide_list(st.session_state.get('target_peptides', '')))
    defaults = [label for label, sequence in label_to_sequence.items() if sequence in current]
    # Render multiselect and map chosen labels back to bare sequences
    selected_labels = st.multiselect(
        "Select peptides", options=labels, default=defaults,
        key=f"{STATE_PREFIX}_pick_multi",
    )
    return [label_to_sequence[label] for label in selected_labels]


def _commit_picked_peptides(selected: list[str], pending_key: str, mode: str) -> None:
    """Stage selected into pending_key (replace or append to the textarea) and rerun."""
    if mode == "append":
        existing = parse_peptide_list(st.session_state.get('target_peptides', ''))
        target_peptides = list(dict.fromkeys(existing + selected))
    else:
        target_peptides = selected
    st.session_state[pending_key] = "\n".join(target_peptides)
    st.rerun()


def _render_pick_actions(selected: list[str], pending_key: str) -> None:
    """Render the 'Set as target' and 'Add to existing' commit buttons."""
    c1, c2 = st.columns(2)
    with c1:
        if st.button("Set as target peptides", key=f"{STATE_PREFIX}_pick_set",
                     type="primary"):
            _commit_picked_peptides(selected, pending_key, "replace")
    with c2:
        if st.button("Add to existing", key=f"{STATE_PREFIX}_pick_add"):
            _commit_picked_peptides(selected, pending_key, "append")


def _render_pick_tab(idxmls: list[str], workers: int, pending_key: str) -> None:
    """Render the pick-from-idXML tab: scan files, then filter and multiselect peptides."""
    # Guard: require idXML selections before showing any pick UI
    if not idxmls:
        st.info("Select idXML files above first.")
        return

    # Ensure cached PSM counts; bail until the user runs a scan
    psm_counts = _ensure_psm_counts_for(idxmls, workers)
    if not psm_counts:
        return

    # Render filters → multiselect → commit-action buttons
    st.caption(f"{len(psm_counts)} unique peptides found")
    filter_text, min_psm = _render_pick_filters(psm_counts)
    selected = _render_pick_multiselect(psm_counts, filter_text, min_psm)
    _render_pick_actions(selected, pending_key)


def render_peptide_input(params: dict) -> None:
    """Render the Target Peptides widget (paste-list + pick-from-idXML tabs) on the Setup page."""
    # Commit any pending peptide list flushed by the pick-tab on the previous run
    pending_key = "_pending_target_peptides"
    persisted_key = get_state_key('target_peptides_persisted')
    if pending_key in st.session_state:
        st.session_state['target_peptides'] = st.session_state.pop(pending_key)

    # Resolve current inputs: idXMLs from state, peptide text from persisted/params fallback
    st.subheader("Target Peptides")

    idxmls = st.session_state.get('selected_idxmls', params.get('selected_idxmls', []))

    # Initial value: persisted session_state (survives page nav) > workspace params
    initial_value = st.session_state.get(persisted_key, params.get('target_peptides', ''))

    # Render the two-tab UI: paste textarea + pick-from-idXML browser
    paste_tab, pick_tab = st.tabs(["Paste list", "Pick from idXML"])
    with paste_tab:
        _render_paste_tab(initial_value, persisted_key)
    with pick_tab:
        _render_pick_tab(idxmls, params.get('workers', DEFAULTS.workers), pending_key)
