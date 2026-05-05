"""Pure-function tests for the Mirror-Plots picker (``mirror_plots.py``).

The Streamlit fragment + StateManager wiring is covered by the AppTest case
in ``test_gui.py``; this file only exercises the data-shaping helpers so a
broken eligibility filter or picker frame fails fast.
"""
from __future__ import annotations

import numpy as np

from mspecviz.search.models import VariantGroup
from mspecviz.spectra.loader_models import LoadedSpectrum
from streamlit_app.views.explore.mirror_plots import (
    _build_picker_frames,
    _eligible_variants,
    _ref_row_id,
)
from tests.conftest import make_peptide_match


def _spec(source_file: str, spectrum_id: str, charge: int = 2,
          rt: float = 12.5, peptidoform: str = "PEPTIDEK",
          match=None, is_ref: bool = False) -> LoadedSpectrum:
    return LoadedSpectrum(
        spectrum_id=spectrum_id,
        mz=np.array([100.0, 200.0]),
        intensity=np.array([1000.0, 2000.0]),
        precursor_mz=400.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence="PEPTIDEK",
        source_file=source_file,
        is_reference=is_ref,
        rt=rt,
        match=match,
    )


def _vg(label: str = "PEPTIDEK") -> VariantGroup:
    return VariantGroup(base_sequence="PEPTIDEK", modification_key="", display_label=label)


class TestEligibleVariants:
    def test_keeps_variants_with_both_query_and_reference(self):
        variant_groups = {"PEPTIDEK|": _vg(), "OTHER|": _vg()}
        reference_store = {"PEPTIDEK|": [(_spec("ref.mzML", "scan=1", is_ref=True), None)]}
        psm_to_vk = {0: "PEPTIDEK|"}
        assert _eligible_variants(variant_groups, reference_store, psm_to_vk) == ["PEPTIDEK|"]

    def test_drops_variant_with_no_reference(self):
        variant_groups = {"PEPTIDEK|": _vg()}
        reference_store: dict = {}
        psm_to_vk = {0: "PEPTIDEK|"}
        assert _eligible_variants(variant_groups, reference_store, psm_to_vk) == []

    def test_drops_variant_with_no_query(self):
        variant_groups = {"PEPTIDEK|": _vg()}
        reference_store = {"PEPTIDEK|": [(_spec("ref.mzML", "scan=1", is_ref=True), None)]}
        psm_to_vk: dict = {}
        assert _eligible_variants(variant_groups, reference_store, psm_to_vk) == []

    def test_preserves_variant_groups_order(self):
        variant_groups = {"AAA|": _vg("AAA"), "BBB|": _vg("BBB"), "CCC|": _vg("CCC")}
        reference_store = {vk: [(_spec("ref.mzML", "scan=1", is_ref=True), None)]
                           for vk in variant_groups}
        psm_to_vk = {0: "AAA|", 1: "BBB|", 2: "CCC|"}
        assert _eligible_variants(variant_groups, reference_store, psm_to_vk) == [
            "AAA|", "BBB|", "CCC|",
        ]


class TestBuildPickerFrames:
    def test_columns_and_keys(self):
        match = make_peptide_match(score=0.875)
        annotated_queries = {7: (_spec("q1.mzML", "scan=42", match=match), None)}
        reference_store = {
            "PEPTIDEK|": [(_spec("r1.mzML", "scan=10", is_ref=True), None)]
        }
        psm_to_vk = {7: "PEPTIDEK|"}

        ref_lazy, query_lazy = _build_picker_frames(
            "PEPTIDEK|", annotated_queries, reference_store, psm_to_vk,
        )
        ref = ref_lazy.collect()
        query = query_lazy.collect()

        assert ref.columns == ["id", "mzml", "scan", "rt", "charge", "score"]
        assert query.columns == ["psm_id", "mzml", "scan", "rt", "charge", "score"]
        assert ref.row(0, named=True) == {
            "id": _ref_row_id("PEPTIDEK|", 0),
            "mzml": "r1", "scan": 10, "rt": 12.5, "charge": 2, "score": None,
        }
        assert query.row(0, named=True) == {
            "psm_id": 7, "mzml": "q1", "scan": 42, "rt": 12.5, "charge": 2, "score": 0.875,
        }

    def test_ref_row_ids_differ_across_peptides(self):
        # The whole point of _ref_row_id: same positional index, different
        # variant_keys must produce different ids — otherwise the openms_insight
        # render-cache key collides on the auto-filled selection and the
        # reference table never refreshes when the peptide changes.
        ids_for_a = [_ref_row_id("VARIANT_A|", i) for i in range(3)]
        ids_for_b = [_ref_row_id("VARIANT_B|", i) for i in range(3)]
        assert set(ids_for_a).isdisjoint(ids_for_b)
        # And the function is deterministic across calls.
        assert _ref_row_id("VARIANT_A|", 0) == ids_for_a[0]

    def test_koina_reference_renders_as_blank_scan(self):
        # Koina spectra use spectrum_id="koina_prediction" which has no digits
        # and no `match`, so scan + score must be None (blank in Tabulator).
        annotated_queries = {0: (_spec("q.mzML", "scan=1", match=make_peptide_match()), None)}
        reference_store = {
            "PEPTIDEK|": [(_spec("Koina:ms2pip", "koina_prediction", rt=0.0, is_ref=True), None)]
        }
        psm_to_vk = {0: "PEPTIDEK|"}

        ref_lazy, _ = _build_picker_frames(
            "PEPTIDEK|", annotated_queries, reference_store, psm_to_vk,
        )
        row = ref_lazy.collect().row(0, named=True)
        # LoadedSpectrum.source_label returns the source_file stem
        # ("Koina:ms2pip") when set, only falling back to bare "Koina" when
        # source_file is empty.
        assert row["mzml"] == "Koina:ms2pip"
        assert row["scan"] is None
        assert row["score"] is None

    def test_filters_queries_by_variant_key(self):
        # Two queries, only one matches the chosen variant.
        annotated_queries = {
            5: (_spec("q.mzML", "scan=5", match=make_peptide_match(peptidoform="A")), None),
            6: (_spec("q.mzML", "scan=6", match=make_peptide_match(peptidoform="B")), None),
        }
        reference_store = {
            "A|": [(_spec("r.mzML", "scan=1", is_ref=True), None)],
            "B|": [(_spec("r.mzML", "scan=2", is_ref=True), None)],
        }
        psm_to_vk = {5: "A|", 6: "B|"}

        _, query_lazy = _build_picker_frames(
            "A|", annotated_queries, reference_store, psm_to_vk,
        )
        df = query_lazy.collect()
        assert df.height == 1
        assert df.row(0, named=True)["psm_id"] == 5

    def test_empty_inputs_return_empty_frames_with_correct_schema(self):
        ref_lazy, query_lazy = _build_picker_frames("MISSING|", {}, {}, {})
        ref = ref_lazy.collect()
        query = query_lazy.collect()
        assert ref.height == 0
        assert query.height == 0
        assert ref.columns == ["id", "mzml", "scan", "rt", "charge", "score"]
        assert query.columns == ["psm_id", "mzml", "scan", "rt", "charge", "score"]

    def test_nan_score_becomes_none(self):
        match = make_peptide_match(score=float("nan"))
        annotated_queries = {0: (_spec("q.mzML", "scan=1", match=match), None)}
        reference_store = {
            "PEPTIDEK|": [(_spec("r.mzML", "scan=1", is_ref=True), None)]
        }
        psm_to_vk = {0: "PEPTIDEK|"}
        _, query_lazy = _build_picker_frames(
            "PEPTIDEK|", annotated_queries, reference_store, psm_to_vk,
        )
        # pl.Float32 cannot represent NaN-as-null distinctly without the explicit
        # conversion in _safe_score; verify it lands as null.
        assert query_lazy.collect().row(0, named=True)["score"] is None
