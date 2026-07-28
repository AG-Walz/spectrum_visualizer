"""Pure-function tests for streamlit_app.views.explore.mirror_plots (picker frames + eligibility)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from psm_utils import PSMList

from streamlit_app.views.explore.mirror_plots import (
    _best_pair_by_metric,
    _build_peptidoform_index,
    _build_picker_frames,
    _download_basename,
    _eligible_peptidoforms,
    _modification_descriptor,
    _query_charges,
    _ref_row_id,
    _static_mirror_bytes,
)
from tests.conftest import make_peptide_match


def _spec(source_file: str, spectrum_id: str, charge: int = 2,
          rt: float = 12.5, peptidoform: str = "PEPTIDEK",
          match=None, is_ref: bool = False, score: float = 0.95):
    """Build a carrier PSM (two-peak spectrum on psm.spectrum) for picker tests.

    run is derived from source_file so source_label(psm) matches; is_ref is kept
    for call-site compatibility but the carrier no longer flags it."""
    if match is not None:
        return match
    psm = make_peptide_match(spectrum_id=spectrum_id, peptidoform=peptidoform,
                             charge=charge, score=score, mzml_filename=source_file)
    psm.retention_time = rt
    psm.spectrum = {"mz": np.array([100.0, 200.0]), "intensity": np.array([1000.0, 2000.0])}
    return psm


def _empty_plist() -> PSMList:
    """Return an empty PSMList — picker tests only care about the dict-key shape."""
    return PSMList(psm_list=[])


class TestEligiblePeptidoforms:
    def test_keeps_peptidoforms_with_both_query_and_reference(self):
        peptidoform_psms = {"PEPTIDEK": _empty_plist(), "OTHER": _empty_plist()}
        reference_store = {"PEPTIDEK": [(_spec("ref.mzML", "scan=1", is_ref=True), None)]}
        psm_id_to_peptidoform = {0: "PEPTIDEK"}
        assert _eligible_peptidoforms(peptidoform_psms, reference_store, psm_id_to_peptidoform) == ["PEPTIDEK"]

    def test_drops_peptidoform_with_no_reference(self):
        peptidoform_psms = {"PEPTIDEK": _empty_plist()}
        reference_store: dict = {}
        psm_id_to_peptidoform = {0: "PEPTIDEK"}
        assert _eligible_peptidoforms(peptidoform_psms, reference_store, psm_id_to_peptidoform) == []

    def test_drops_peptidoform_with_no_query(self):
        peptidoform_psms = {"PEPTIDEK": _empty_plist()}
        reference_store = {"PEPTIDEK": [(_spec("ref.mzML", "scan=1", is_ref=True), None)]}
        psm_id_to_peptidoform: dict = {}
        assert _eligible_peptidoforms(peptidoform_psms, reference_store, psm_id_to_peptidoform) == []

    def test_preserves_peptidoform_psms_order(self):
        peptidoform_psms = {"AAA": _empty_plist(), "BBB": _empty_plist(), "CCC": _empty_plist()}
        reference_store = {vk: [(_spec("ref.mzML", "scan=1", is_ref=True), None)]
                           for vk in peptidoform_psms}
        psm_id_to_peptidoform = {0: "AAA", 1: "BBB", 2: "CCC"}
        assert _eligible_peptidoforms(peptidoform_psms, reference_store, psm_id_to_peptidoform) == [
            "AAA", "BBB", "CCC",
        ]


class TestBuildPickerFrames:
    def test_columns_and_keys(self):
        annotated_queries = {7: (_spec("q1.mzML", "scan=42", score=0.875), None)}
        reference_store = {
            # score 0.5 is exact in Float32 (the picker's score dtype)
            "PEPTIDEK": [(_spec("r1.mzML", "scan=10", is_ref=True, score=0.5), None)]
        }
        psm_id_to_peptidoform = {7: "PEPTIDEK"}

        ref_lazy, query_lazy = _build_picker_frames(
            "PEPTIDEK", annotated_queries, reference_store, psm_id_to_peptidoform,
        )
        ref = ref_lazy.collect()
        query = query_lazy.collect()

        assert ref.columns == ["id", "mzml", "scan", "rt", "charge", "score"]
        assert query.columns == ["psm_id", "mzml", "scan", "rt", "charge", "score"]
        assert ref.row(0, named=True) == {
            "id": _ref_row_id("PEPTIDEK", 0),
            # Experimental references are PSM-backed, so they now carry a score.
            "mzml": "r1", "scan": 10, "rt": 12.5, "charge": 2, "score": 0.5,
        }
        assert query.row(0, named=True) == {
            "psm_id": 7, "mzml": "q1", "scan": 42, "rt": 12.5, "charge": 2, "score": 0.875,
        }

    def test_ref_row_ids_differ_across_peptides(self):
        # The whole point of _ref_row_id: same positional index across
        # different peptidoforms must produce DIFFERENT ids — otherwise the
        # openms_insight render-cache key collides on the auto-filled
        # selection and the reference table never refreshes when the
        # peptide changes.
        ids_for_a = [_ref_row_id("PEPTIDOFORM_A|", i) for i in range(3)]
        ids_for_b = [_ref_row_id("PEPTIDOFORM_B|", i) for i in range(3)]
        assert set(ids_for_a).isdisjoint(ids_for_b)
        # And the function must be deterministic across calls.
        assert _ref_row_id("PEPTIDOFORM_A|", 0) == ids_for_a[0]

    def test_koina_reference_renders_as_blank_scan(self):
        # Koina spectra use spectrum_id="koina_prediction" (no digits) and have
        # no search score ⇒ scan + score must be None so the cell renders blank.
        annotated_queries = {0: (_spec("q.mzML", "scan=1"), None)}
        reference_store = {
            "PEPTIDEK": [(_spec("Koina:ms2pip", "koina_prediction", rt=0.0, is_ref=True, score=None), None)]
        }
        psm_id_to_peptidoform = {0: "PEPTIDEK"}

        ref_lazy, _ = _build_picker_frames(
            "PEPTIDEK", annotated_queries, reference_store, psm_id_to_peptidoform,
        )
        row = ref_lazy.collect().row(0, named=True)
        # source_label returns the PSM run ("Koina:ms2pip") for Koina references.
        assert row["mzml"] == "Koina:ms2pip"
        assert row["scan"] is None
        assert row["score"] is None

    def test_filters_queries_by_peptidoform(self):
        # Two queries, only one matches the chosen peptidoform — the other is dropped.
        annotated_queries = {
            5: (_spec("q.mzML", "scan=5", peptidoform="A"), None),
            6: (_spec("q.mzML", "scan=6", peptidoform="B"), None),
        }
        reference_store = {
            "A": [(_spec("r.mzML", "scan=1", is_ref=True), None)],
            "B": [(_spec("r.mzML", "scan=2", is_ref=True), None)],
        }
        psm_id_to_peptidoform = {5: "A", 6: "B"}

        _, query_lazy = _build_picker_frames(
            "A", annotated_queries, reference_store, psm_id_to_peptidoform,
        )
        df = query_lazy.collect()
        assert df.height == 1
        assert df.row(0, named=True)["psm_id"] == 5

    def test_empty_inputs_return_empty_frames_with_correct_schema(self):
        ref_lazy, query_lazy = _build_picker_frames("MISSING", {}, {}, {})
        ref = ref_lazy.collect()
        query = query_lazy.collect()
        assert ref.height == 0
        assert query.height == 0
        assert ref.columns == ["id", "mzml", "scan", "rt", "charge", "score"]
        assert query.columns == ["psm_id", "mzml", "scan", "rt", "charge", "score"]

    def test_nan_score_becomes_none(self):
        annotated_queries = {0: (_spec("q.mzML", "scan=1", score=float("nan")), None)}
        reference_store = {
            "PEPTIDEK": [(_spec("r.mzML", "scan=1", is_ref=True), None)]
        }
        psm_id_to_peptidoform = {0: "PEPTIDEK"}
        _, query_lazy = _build_picker_frames(
            "PEPTIDEK", annotated_queries, reference_store, psm_id_to_peptidoform,
        )
        # pl.Float32 cannot represent NaN-as-null distinctly without the
        # explicit conversion in _safe_score — verify it lands as null.
        assert query_lazy.collect().row(0, named=True)["score"] is None

    def test_charge_filter_drops_nonmatching_query(self):
        annotated_queries = {
            1: (_spec("q.mzML", "scan=1", charge=2), None),
            2: (_spec("q.mzML", "scan=2", charge=3), None),
        }
        reference_store = {"PEP": [(_spec("r.mzML", "scan=1", is_ref=True), None)]}
        psm_id_to_peptidoform = {1: "PEP", 2: "PEP"}
        _, query_lazy = _build_picker_frames(
            "PEP", annotated_queries, reference_store, psm_id_to_peptidoform, [2],
        )
        rows = query_lazy.collect()
        assert rows.height == 1 and rows.row(0, named=True)["charge"] == 2


class TestPeptidoformIndex:
    def test_splits_and_sorts_by_sequence_then_modification(self):
        # #106: sequences sort alphabetically; #62: variants group under their base sequence.
        eligible = ["VVM[Oxidation]EYENL", "AAIVGKQVL", "VVMEYENL"]
        sequence_to_keys, descriptor_of = _build_peptidoform_index(eligible)
        assert sorted(sequence_to_keys) == ["AAIVGKQVL", "VVMEYENL"]
        # "(unmodified)" sorts before "M3(Oxidation)".
        assert sequence_to_keys["VVMEYENL"] == ["VVMEYENL", "VVM[Oxidation]EYENL"]
        assert descriptor_of["VVMEYENL"] == "(unmodified)"

    def test_modification_descriptor(self):
        from psm_utils import Peptidoform
        assert _modification_descriptor(Peptidoform("AAIVGKQVL")) == "(unmodified)"
        assert _modification_descriptor(Peptidoform("VVM[Oxidation]EYENL")) == "M3(Oxidation)"
        assert _modification_descriptor(Peptidoform("[Acetyl]-PEPTIDEK")) == "N-term(Acetyl)"


class TestQueryChargesAndBestPair:
    def test_query_charges_sorted_distinct_for_peptidoform(self):
        annotated_queries = {
            1: (_spec("q", "s1", charge=3), None),
            2: (_spec("q", "s2", charge=2), None),
            3: (_spec("q", "s3", charge=2), None),
        }
        psm_id_to_peptidoform = {1: "PEP", 2: "PEP", 3: "OTHER"}
        assert _query_charges("PEP", annotated_queries, psm_id_to_peptidoform) == [2, 3]

    def test_best_pair_picks_highest_metric(self):
        # Reads precomputed metric from comparison_frame and resolves the row back to a store index.
        reference_a = _spec("rA.mzML", "scan=1", is_ref=True)
        reference_b = _spec("rB.mzML", "scan=2", is_ref=True)
        reference_store = {"PEP": [(reference_a, None), (reference_b, None)]}
        comparison_frame = pd.DataFrame([
            {"psm_id": 9, "Charge": 2, "Reference": "rA", "ref_spectrum_id": reference_a.spectrum_id, "R² match": 0.99},
            {"psm_id": 9, "Charge": 2, "Reference": "rB", "ref_spectrum_id": reference_b.spectrum_id, "R² match": 0.40},
        ])
        r2_value, psm_id, reference_index = _best_pair_by_metric(
            comparison_frame, {9}, reference_store, "PEP", [], "R² match",
        )
        assert psm_id == 9 and reference_index == 0 and r2_value > 0.98

    def test_best_pair_ranks_by_any_similarity_metric(self):
        # The auto-pick generalises beyond R²: ranking by Cosine picks a different winner here.
        reference_a = _spec("rA.mzML", "scan=1", is_ref=True)
        reference_b = _spec("rB.mzML", "scan=2", is_ref=True)
        reference_store = {"PEP": [(reference_a, None), (reference_b, None)]}
        comparison_frame = pd.DataFrame([
            {"psm_id": 9, "Charge": 2, "Reference": "rA", "ref_spectrum_id": reference_a.spectrum_id,
             "R² match": 0.99, "Cosine": 0.20},
            {"psm_id": 9, "Charge": 2, "Reference": "rB", "ref_spectrum_id": reference_b.spectrum_id,
             "R² match": 0.40, "Cosine": 0.95},
        ])
        _, _, reference_index = _best_pair_by_metric(
            comparison_frame, {9}, reference_store, "PEP", [], "Cosine",
        )
        assert reference_index == 1  # rB wins on Cosine, though rA wins on R²

    def test_best_pair_disambiguates_colliding_spectrum_ids_by_source(self):
        # spectrum_id is not unique across files: two references share "scan=1" but different mzML.
        # The winning row must resolve to the reference from its own source file, not the last one.
        reference_a = _spec("rA.mzML", "scan=1", is_ref=True)
        reference_b = _spec("rB.mzML", "scan=1", is_ref=True)
        reference_store = {"PEP": [(reference_a, None), (reference_b, None)]}
        comparison_frame = pd.DataFrame([
            {"psm_id": 9, "Charge": 2, "Reference": "rA", "ref_spectrum_id": "scan=1", "R² match": 0.95},
            {"psm_id": 9, "Charge": 2, "Reference": "rB", "ref_spectrum_id": "scan=1", "R² match": 0.10},
        ])
        _, _, reference_index = _best_pair_by_metric(
            comparison_frame, {9}, reference_store, "PEP", [], "R² match",
        )
        assert reference_index == 0  # rA, not the last-written rB

    def test_best_pair_none_when_charge_filter_excludes_all(self):
        reference = _spec("r.mzML", "scan=1", is_ref=True)
        reference_store = {"PEP": [(reference, None)]}
        comparison_frame = pd.DataFrame([
            {"psm_id": 9, "Charge": 2, "Reference": "r", "ref_spectrum_id": reference.spectrum_id, "R² match": 0.9},
        ])
        assert _best_pair_by_metric(comparison_frame, {9}, reference_store, "PEP", [3], "R² match") is None


class TestMirrorDownloads:
    """The Explore mirror plot can be saved straight from the page (issue #107)."""

    def _frame(self) -> pd.DataFrame:
        return pd.DataFrame({
            "mz": [100.0, 200.0], "intensity": [1000.0, 2000.0],
            "annotation": ["b1", "y1"], "ion_type": ["b", "y"],
        })

    def test_renders_png_and_svg_payloads(self):
        # The undecorated renderer — the module caches it separately for the download buttons.
        png_bytes, svg_bytes = _static_mirror_bytes(
            self._frame(), self._frame(), title="PEPTIDEK +2",
            query_label="Query", reference_label="Ref", metrics_text="",
            max_peak_labels=0, export_resolution=100,
        )
        assert png_bytes[:4] == b"\x89PNG"
        assert b"<svg" in svg_bytes[:1000]

    def test_basename_identifies_both_spectra(self):
        query = _spec("query.mzML", "scan=7")
        reference = _spec("ref.mzML", "scan=9", is_ref=True)
        basename = _download_basename(query, reference)
        assert "PEPTIDEK" in basename and "z2" in basename
        assert "query" in basename and "ref" in basename
        # sanitize_filename strips the characters an OS would reject in a filename.
        assert not set(basename) & set('<>:"/\\|?*[]')
