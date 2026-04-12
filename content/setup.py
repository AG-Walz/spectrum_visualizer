"""
Setup Page — File selection, peptide input, reference config, annotation settings.
"""
import time as _time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import streamlit as st
from psm_utils.io.idxml import IdXMLReader

from src.common.streamlit_setup import page_setup, save_params
from content.page_state import (
    STATE_PREFIX, get_state_key, _get_settings, restore_analysis_if_needed,
    list_subdirectories, list_files,
)
from src.common.text_utils import get_nested as get_settings_value, parse_peptide_list

try:
    from src.koina_prediction import (
        KOINA_AVAILABLE, KOINA_MODELS, INSTRUMENT_TYPES,
        DEFAULT_MODEL, DEFAULT_COLLISION_ENERGY, DEFAULT_INSTRUMENT_TYPE,
        get_koina_settings,
    )
except ImportError:
    KOINA_AVAILABLE = False
    KOINA_MODELS = {}
    INSTRUMENT_TYPES = []
    DEFAULT_MODEL = ""
    DEFAULT_COLLISION_ENERGY = 35
    DEFAULT_INSTRUMENT_TYPE = "TIMSTOF"

    def get_koina_settings(s):
        return {"enabled": False, "model": DEFAULT_MODEL,
                "collision_energy": DEFAULT_COLLISION_ENERGY,
                "instrument_type": DEFAULT_INSTRUMENT_TYPE}

params = page_setup()


def _scan_idxml_psm_counts(idxml_path: str) -> dict:
    """Worker: count PSMs per peptide in one idXML file. Module-level for pickling."""
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            seq = psm.peptidoform.sequence.upper()
            if seq:
                counts[seq] += 1
    except Exception:
        pass
    return dict(counts)


def render_directory_browser(key_prefix: str, extensions: list, state_key: str, default_dir: str = "/"):
    """Render a directory browser with file selection."""
    dir_key = f"{key_prefix}_dir"
    nav_key = f"{key_prefix}_nav_to"
    select_all_key = f"{key_prefix}_select_all"
    recursive_key = f"{key_prefix}_recursive"

    if state_key not in st.session_state:
        st.session_state[state_key] = params.get(state_key, [])

    if dir_key not in st.session_state:
        st.session_state[dir_key] = default_dir

    if nav_key in st.session_state and st.session_state[nav_key]:
        st.session_state[dir_key] = st.session_state[nav_key]
        st.session_state[nav_key] = None
        if select_all_key in st.session_state:
            del st.session_state[select_all_key]
        if recursive_key in st.session_state:
            del st.session_state[recursive_key]

    current_dir = st.session_state[dir_key]

    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path", value=current_dir,
            key=f"{key_prefix}_path_input_{abs(hash(current_dir))}",
            label_visibility="collapsed", placeholder="Enter path..."
        )
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[dir_key] = new_path
            if select_all_key in st.session_state:
                del st.session_state[select_all_key]
            if recursive_key in st.session_state:
                del st.session_state[recursive_key]
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("up", key=f"{key_prefix}_up_btn", help="Parent directory",
                      disabled=(parent == current_dir)):
            st.session_state[nav_key] = parent
            st.rerun()

    subdir_paths = list_subdirectories(current_dir)
    subdirs = [Path(d) for d in subdir_paths]
    if subdirs:
        st.markdown(
            "<style>.folder-nav button{font-size:0.8rem !important;"
            "padding:0.15rem 0.4rem !important;min-height:0 !important}</style>",
            unsafe_allow_html=True,
        )
        cols_per_row = 4
        max_folders = 16
        for i in range(0, min(len(subdirs), max_folders), cols_per_row):
            cols = st.columns(cols_per_row)
            for j, col in enumerate(cols):
                idx = i + j
                if idx < len(subdirs) and idx < max_folders:
                    folder = subdirs[idx]
                    with col:
                        if st.button(f"{folder.name}", key=f"{key_prefix}_folder_{idx}",
                                     help=str(folder)):
                            st.session_state[nav_key] = str(folder)
                            st.rerun()
        if len(subdirs) > max_folders:
            st.caption(f"... and {len(subdirs) - max_folders} more folders")

    recursive = st.checkbox("Include subdirectories", value=False,
                            key=f"{key_prefix}_recursive",
                            help="Scan all subdirectories for matching files")

    all_files = [Path(f) for f in list_files(current_dir, tuple(extensions), recursive)]
    current_path = Path(current_dir)
    file_labels = ([str(f.relative_to(current_path)) for f in all_files]
                   if recursive else [f.name for f in all_files])
    label_to_path = {lbl: str(f) for lbl, f in zip(file_labels, all_files)}

    filter_text = st.text_input("Filter files", key=f"{key_prefix}_filter",
                                placeholder="Type to filter filenames...")
    if filter_text:
        filter_lower = filter_text.lower()
        filtered_labels = [lbl for lbl in file_labels if filter_lower in lbl.lower()]
    else:
        filtered_labels = file_labels

    if filtered_labels:
        all_selected_paths = set(st.session_state.get(state_key, []))
        current_scope_paths = {label_to_path[lbl] for lbl in file_labels}
        if recursive:
            other_dir_selections = [p for p in st.session_state.get(state_key, [])
                                    if p not in current_scope_paths]
        else:
            other_dir_selections = [p for p in st.session_state.get(state_key, [])
                                    if Path(p).parent != current_path]

        valid_defaults = [lbl for lbl in filtered_labels if label_to_path[lbl] in all_selected_paths]
        all_in_scope_selected = set(valid_defaults) == set(filtered_labels) and len(valid_defaults) > 0

        select_all = st.checkbox(
            f"Select all {len(filtered_labels)} files" +
            (f" (filtered from {len(file_labels)})" if filter_text else ""),
            value=all_in_scope_selected, key=select_all_key
        )

        if select_all:
            selected_labels = filtered_labels
        else:
            selected_labels = st.multiselect(
                "Select files" + (" (recursive):" if recursive else ":"),
                options=filtered_labels,
                default=valid_defaults if not all_in_scope_selected else [],
                key=f"{key_prefix}_files"
            )

        current_selected_paths = [label_to_path[lbl] for lbl in selected_labels]
        merged_selections = list(dict.fromkeys(other_dir_selections + current_selected_paths))
        st.session_state[state_key] = merged_selections

    else:
        if filter_text:
            st.info(f"No files match filter '{filter_text}'")
        else:
            st.info(f"No {extensions[0]} files found")

    all_selected = st.session_state.get(state_key, [])
    if all_selected:
        if st.button(f"Clear all selections ({len(all_selected)} files)",
                     key=f"{key_prefix}_clear_all", width="stretch"):
            st.session_state[state_key] = []
            if select_all_key in st.session_state:
                del st.session_state[select_all_key]
            st.rerun()


def file_selection_fragment():
    """File browser for idXML and mzML files."""
    _settings = _get_settings()
    data_dir_default = get_settings_value(_settings, "data_directories", "default", default="/")

    idxml_dir = params.get('idxml-browser-dir', data_dir_default)
    mzml_dir = params.get('mzml-browser-dir', data_dir_default)

    st.header("File Selection")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_idxml",
                                 extensions=['.idXML', '.idxml'],
                                 state_key='selected-idxmls',
                                 default_dir=idxml_dir)
    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_mzml",
                                 extensions=['.mzML', '.mzml'],
                                 state_key='selected-mzmls',
                                 default_dir=mzml_dir)

    for dir_param, dir_ss_key, fallback in [
        ('idxml-browser-dir', f"{STATE_PREFIX}_idxml_dir", idxml_dir),
        ('mzml-browser-dir', f"{STATE_PREFIX}_mzml_dir", mzml_dir),
    ]:
        val = st.session_state.get(dir_ss_key, fallback)
        params[dir_param] = val
        st.session_state[dir_param] = val

    save_params(params)


def peptide_input_fragment():
    """Peptide input fragment with paste and idXML picker tabs."""
    _pending_key = "_pending_target_peptides"
    _persisted_key = get_state_key('target_peptides_persisted')
    if _pending_key in st.session_state:
        st.session_state['target-peptides'] = st.session_state.pop(_pending_key)

    st.subheader("Target Peptides")

    idxmls = st.session_state.get('selected-idxmls', params.get('selected-idxmls', []))

    # Resolve initial value: persisted session_state (survives page nav) > workspace params
    initial_value = st.session_state.get(_persisted_key, params.get('target-peptides', ''))

    paste_tab, pick_tab = st.tabs(["Paste list", "Pick from idXML"])

    with paste_tab:
        col1, col2 = st.columns([2, 1])
        with col1:
            peptide_text = st.text_area(
                "Enter peptide sequences (one per line)",
                value=initial_value,
                height=150, key="target-peptides",
                placeholder="PEPTIDESEQ\nANOTHERPEPTIDE\n..."
            )
        # Persist in a non-widget key that survives page navigation
        st.session_state[_persisted_key] = peptide_text
        with col2:
            st.markdown("**Tips:**\n- One peptide per line\n"
                         "- Without modifications\n"
                         "- Case insensitive")
            peptides = parse_peptide_list(peptide_text)
            if peptides:
                st.success(f"{len(peptides)} peptide(s) entered")

    with pick_tab:
        if not idxmls:
            st.info("Select idXML files above first.")
        else:
            scan_key = '_idxml_psm_counts'
            scan_files_key = '_idxml_scan_files'
            psm_counts = st.session_state.get(scan_key, {})
            scanned_files = st.session_state.get(scan_files_key)

            if not psm_counts or scanned_files != idxmls:
                n = len(idxmls)
                st.info(f"{n} idXML file{'s' if n != 1 else ''} selected")
                workers = params.get('workers', 4)
                if st.button(f"Scan {n} files for peptides",
                             key=f"{STATE_PREFIX}_scan_idxmls", type="secondary"):
                    merged = Counter()
                    errors = []
                    with st.status(f"Scanning {n} idXML files ({workers} workers)...",
                                   expanded=True) as status:
                        progress = st.progress(0)
                        info = st.empty()
                        t0 = _time.time()

                        with ProcessPoolExecutor(max_workers=min(workers, n)) as pool:
                            futures = {
                                pool.submit(_scan_idxml_psm_counts, p): p for p in idxmls
                            }
                            done = 0
                            for future in as_completed(futures):
                                done += 1
                                path = futures[future]
                                try:
                                    merged.update(future.result())
                                except Exception as e:
                                    errors.append(f"{Path(path).name}: {e}")
                                elapsed = _time.time() - t0
                                remaining = (elapsed / done) * (n - done) if done < n else 0
                                progress.progress(done / n,
                                    f"{done}/{n} files — {elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining")
                                info.caption(f"Scanned {Path(path).name}")

                        progress.progress(1.0, f"Done — {len(merged)} peptides from {n} files in {_time.time()-t0:.1f}s")
                        if errors:
                            status.update(label=f"Completed with {len(errors)} errors", state="error")
                            for e in errors[:5]:
                                st.caption(e)
                        else:
                            status.update(label=f"{len(merged)} peptides from {n} files", state="complete")

                    psm_counts = dict(sorted(merged.items()))
                    st.session_state[scan_key] = psm_counts
                    st.session_state[scan_files_key] = idxmls
                    st.rerun()

            if psm_counts:
                st.caption(f"{len(psm_counts)} unique peptides found")

                f_col1, f_col2 = st.columns([2, 1])
                with f_col1:
                    filter_text = st.text_input(
                        "Filter", placeholder="Type to filter by substring...",
                        key=f"{STATE_PREFIX}_pick_filter",
                        label_visibility="collapsed",
                    )
                with f_col2:
                    max_count = max(psm_counts.values())
                    min_psm = st.number_input(
                        "Min PSMs", min_value=1, max_value=max_count,
                        value=1, step=1,
                        key=f"{STATE_PREFIX}_pick_min_psm",
                        help="Show only peptides with at least this many PSMs",
                    )

                shown_items = [
                    (seq, cnt) for seq, cnt in psm_counts.items()
                    if cnt >= min_psm
                    and (not filter_text or filter_text.upper() in seq)
                ]
                labels = [f"{seq} ({cnt})" for seq, cnt in shown_items]
                label_to_seq = {f"{seq} ({cnt})": seq for seq, cnt in shown_items}

                current = set(parse_peptide_list(
                    st.session_state.get('target-peptides', '')
                ))
                defaults = [lbl for lbl, seq in label_to_seq.items() if seq in current]

                selected_labels = st.multiselect(
                    "Select peptides",
                    options=labels,
                    default=defaults,
                    key=f"{STATE_PREFIX}_pick_multi",
                )
                selected = [label_to_seq[lbl] for lbl in selected_labels]

                c1, c2 = st.columns(2)
                with c1:
                    if st.button("Set as target peptides", key=f"{STATE_PREFIX}_pick_set",
                                 type="primary"):
                        st.session_state[_pending_key] = "\n".join(selected)
                        st.rerun()
                with c2:
                    if st.button("Add to existing", key=f"{STATE_PREFIX}_pick_add"):
                        merged = list(dict.fromkeys(
                            parse_peptide_list(
                                st.session_state.get('target-peptides', '')
                            ) + selected
                        ))
                        st.session_state[_pending_key] = "\n".join(merged)
                        st.rerun()


def reference_selection_fragment():
    """Reference substring configuration (not a fragment — must update when file selection changes)."""
    st.subheader("Reference Identification")
    mzmls = st.session_state.get('selected-mzmls', params.get('selected-mzmls', []))
    if not mzmls:
        st.warning("Please select mzML files first")
        return

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")
    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=params.get('reference-substring', ''),
        key="reference-substring", placeholder="e.g., Synpep"
    )

    if ref_substring:
        reference_files = [Path(m).name for m in mzmls if ref_substring.lower() in Path(m).name.lower()]
        query_files = [Path(m).name for m in mzmls if ref_substring.lower() not in Path(m).name.lower()]

        with st.expander(f"Reference ({len(reference_files)}) / Query ({len(query_files)}) files", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Reference files ({len(reference_files)})**")
                for f in reference_files[:5]:
                    st.caption(f)
            with col2:
                st.markdown(f"**Query files ({len(query_files)})**")
                for f in query_files[:5]:
                    st.caption(f)

        if reference_files and query_files:
            st.success("Reference/Query separation configured")
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")
    else:
        st.info("No reference identifier -- individual annotated spectrum plots only.")

    save_params(params)


def settings_fragment():
    """Annotation, display, performance, and Koina settings fragment."""
    st.subheader("Settings")
    available_ions = get_settings_value(_get_settings(), "ion_types", "available",
                                        default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Ion Types**")
        current_ions = params.get('ion-types', ['b', 'y'])
        st.multiselect("Select ion types", options=available_ions,
                        default=[i for i in current_ions if i in available_ions],
                        key="ion-types")

    with col2:
        st.markdown("**m/z Tolerance**")
        st.number_input("Tolerance (Da)", min_value=0.001, max_value=1.0,
                        value=float(params.get('mz-tolerance', 0.02)),
                        step=0.005, format="%.3f", key="mz-tolerance")

    st.markdown("**Annotation Settings**")
    ann_col1, ann_col2, ann_col3 = st.columns(3)

    with ann_col1:
        st.checkbox("Add neutral losses (H2O, NH3)",
                    value=params.get('add-losses', False),
                    key="add-losses")

    with ann_col2:
        st.checkbox("Add first prefix ion",
                    value=params.get('add-first-prefix-ion', True),
                    key="add-first-prefix-ion")

    with ann_col3:
        st.checkbox("Remove precursor peaks",
                    value=params.get('remove-precursor-peaks', False),
                    key="remove-precursor-peaks")

    st.radio("**Intensity Display Mode**",
             options=["Relative", "Raw"],
             captions=["Relative to the highest intensity", "Raw intensity"],
             index=0 if params.get('intensity-mode', 'Relative') == 'Relative' else 1,
             horizontal=True, key="intensity-mode",
             label_visibility="visible")

    st.markdown("**Performance**")
    st.slider("Parallel workers", min_value=1, max_value=32,
              value=params.get('workers', 4),
              key="workers",
              help="Number of parallel processes for plot generation (~30 MB per worker)")

    st.markdown("---")
    st.markdown("**Koina AI Prediction**")
    st.caption(
        "When enabled, spectra without an experimental reference are compared "
        "against AI-predicted fragmentation patterns (intensity + m/z)."
    )

    koina_defaults = get_koina_settings(_get_settings())

    koina_enabled = st.checkbox(
        "Enable Koina fallback reference",
        value=params.get("koina-enabled", koina_defaults["enabled"]),
        key="koina-enabled",
        help="Use AI-predicted spectrum as reference when no experimental reference mzML exists",
        disabled=not KOINA_AVAILABLE,
    )

    if not KOINA_AVAILABLE:
        st.warning("koinapy is not installed. Install it to enable Koina predictions.")

    if koina_enabled and KOINA_AVAILABLE:
        k_col1, k_col2, k_col3 = st.columns(3)

        with k_col1:
            model_names = list(KOINA_MODELS.keys())
            model_ids = list(KOINA_MODELS.values())
            current_model_id = params.get("koina-model", koina_defaults["model"])
            default_idx = model_ids.index(current_model_id) if current_model_id in model_ids else 0

            selected_model_name = st.selectbox(
                "Prediction Model", options=model_names, index=default_idx,
                key=f"{STATE_PREFIX}_koina_model_select",
                help="Deep-learning model for fragment ion prediction",
            )
            st.session_state["koina-model"] = KOINA_MODELS[selected_model_name]

        with k_col2:
            st.number_input(
                "Collision Energy (eV)", min_value=1, max_value=100,
                value=params.get("koina-ce", koina_defaults["collision_energy"]),
                key="koina-ce",
            )

        with k_col3:
            current_inst = params.get("koina-instrument", koina_defaults["instrument_type"])
            inst_idx = INSTRUMENT_TYPES.index(current_inst) if current_inst in INSTRUMENT_TYPES else 0

            st.selectbox(
                "Instrument Type", options=INSTRUMENT_TYPES, index=inst_idx,
                key="koina-instrument",
            )

    save_params(params)


def render():
    restore_analysis_if_needed()
    st.title("Setup")
    st.caption("Configure files, peptides, reference, and annotation settings")

    file_selection_fragment()
    st.markdown("---")
    peptide_input_fragment()
    st.markdown("---")
    reference_selection_fragment()
    st.markdown("---")
    settings_fragment()

    save_params(params)


render()
