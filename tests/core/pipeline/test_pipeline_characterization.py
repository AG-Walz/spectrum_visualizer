"""End-to-end characterization of the CLI default (single-spectrum) pipeline.

This test pins the *observable* output of the full lookup -> stream+annotate ->
build-frames pipeline on the small_dataset fixture, mirroring what
`python -m cli -c cli/config/config_example.yaml` runs (reference_substring="",
Koina disabled). It is the safety net for the PSMList-native migration: every
refactor step that reshapes how spectra and PSMs are carried through the pipeline
must keep these invariants green. Assertions use range checks for peak counts
(matching the spectrum-loader integration test) so harmless library upgrades
don't read as regressions, but pin schemas and row counts exactly.
"""
import copy

import pytest

from core.dataframes import (
    HEATMAP_SCHEMA,
    PEAKS_SCHEMA,
    PSM_TABLE_SCHEMA,
    SEQUENCE_SCHEMA,
    build_psm_table,
    build_visualization_frames,
)
from core.peptide_lookup import lookup_peptides_in_idxmls
from core.pipeline import (
    build_annotation_kwargs,
    build_koina_kwargs,
    run_streaming_annotation,
)


@pytest.mark.integration
class TestSingleSpectrumPipelineCharacterization:
    @pytest.fixture(scope="class")
    def pipeline_output(self, idxml_path, mzml_paths, target_peptides):
        """Run the whole CLI default pipeline once; share results across tests."""
        variant_groups = lookup_peptides_in_idxmls(
            [idxml_path], target_peptides, mzml_paths=mzml_paths,
        )
        # The CLI deep-copies before annotating because streaming stamps PSMs.
        variant_groups = copy.deepcopy(variant_groups)
        psm_frame, psm_lookup = build_psm_table(variant_groups)

        # getter returning the default keeps Koina off and uses DEFAULTS settings.
        getter = lambda key, default: default
        annotated_queries, reference_store, comparison_rows = run_streaming_annotation(
            variant_groups=variant_groups,
            mzml_paths=mzml_paths,
            reference_substring="",
            psm_lookup=psm_lookup,
            annotation_kwargs=build_annotation_kwargs(getter),
            koina_kwargs=build_koina_kwargs(getter),
            workers=1,
        )
        frames = build_visualization_frames(annotated_queries, psm_frame)
        return {
            "variant_groups": variant_groups,
            "psm_frame": psm_frame.collect(),
            "annotated_queries": annotated_queries,
            "reference_store": reference_store,
            "comparison_rows": comparison_rows,
            "frames": {name: lf.collect() for name, lf in frames.items()},
        }

    # --- Step 1: lookup ----------------------------------------------------

    def test_lookup_finds_two_psms_in_two_variants(self, pipeline_output):
        variant_groups = pipeline_output["variant_groups"]
        assert set(variant_groups) == {"ASDDFDQASY", "RTEFNLNQY"}
        assert sum(len(group) for group in variant_groups.values()) == 2

    # --- PSM selector table ------------------------------------------------

    def test_psm_table_schema_and_rows(self, pipeline_output):
        psm_table = pipeline_output["psm_frame"]
        assert psm_table.columns == list(PSM_TABLE_SCHEMA)
        assert psm_table.height == 2
        assert set(psm_table["id"].to_list()) == {0, 1}
        assert set(psm_table["base_sequence"].to_list()) == {"ASDDFDQASY", "RTEFNLNQY"}
        assert all(charge > 0 for charge in psm_table["charge"].to_list())

    # --- Step 2: stream + annotate ----------------------------------------

    def test_two_query_spectra_annotated(self, pipeline_output):
        annotated_queries = pipeline_output["annotated_queries"]
        assert len(annotated_queries) == 2
        for _psm_id, _spectrum, query_frame in annotated_queries:
            assert {"mz", "intensity", "annotation", "ion_type", "is_annotated"} <= set(query_frame.columns)
            # Every PSM should match at least one theoretical fragment.
            assert query_frame["is_annotated"].sum() > 0

    def test_single_spectrum_mode_has_no_references(self, pipeline_output):
        # reference_substring="" + Koina disabled -> no experimental or predicted refs.
        assert pipeline_output["reference_store"] == {}
        assert pipeline_output["comparison_rows"] == []

    # --- Visualization frames ---------------------------------------------

    def test_frame_bundle_keys(self, pipeline_output):
        assert set(pipeline_output["frames"]) == {"peaks", "sequence", "heatmap", "psm"}

    def test_peaks_frame_schema_and_count(self, pipeline_output):
        peaks = pipeline_output["frames"]["peaks"]
        assert peaks.columns == list(PEAKS_SCHEMA)
        # Total peaks across both spectra — range check (~373), like the loader test.
        assert 350 < peaks.height < 400
        assert peaks["highlight"].sum() > 0

    def test_sequence_frame_schema_and_rows(self, pipeline_output):
        sequence = pipeline_output["frames"]["sequence"]
        assert sequence.columns == list(SEQUENCE_SCHEMA)
        assert sequence.height == 2
        assert all(charge > 0 for charge in sequence["precursor_charge"].to_list())

    def test_heatmap_frame_schema_and_count(self, pipeline_output):
        heatmap = pipeline_output["frames"]["heatmap"]
        assert heatmap.columns == list(HEATMAP_SCHEMA)
        # Heatmap has one row per peak, same source as the peaks frame.
        assert 350 < heatmap.height < 400
