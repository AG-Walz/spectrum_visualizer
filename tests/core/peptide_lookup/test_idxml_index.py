"""Tests for core.idxml_index: the index gives the same counts and PSMs as a full parse,
is reused while the idXML is unchanged, and keeps reporting parse errors."""
import os

import pytest

from core.idxml_index import (
    count_indexed_peptides,
    ensure_idxml_index,
    index_file_for,
    load_indexed_psms,
    missing_from_index,
)
from core.peptide_lookup import count_psms_per_peptide, load_idxmls
from tests.conftest import SMALL_DATASET

# Smallest example idXML (19 MB, 2194 PSMs) keeps the integration test to a few seconds
SMALL_IDXML = str(SMALL_DATASET / "idXML" / "ProT002_1_pout_filtered.idXML")


def _sorted_proteins(psm):
    """psm_utils builds protein_list from a set, so its order differs between processes."""
    return psm.model_copy(update={"protein_list": sorted(psm.protein_list or [])})


@pytest.mark.integration
class TestIndexMatchesFullParse:
    """Counts and target PSMs read from the index equal those of core.peptide_lookup."""

    def test_counts_and_psms_equal_full_parse(self, tmp_path):
        index_files = ensure_idxml_index([SMALL_IDXML], tmp_path)
        parsed = load_idxmls([SMALL_IDXML])
        assert count_indexed_peptides(index_files) == count_psms_per_peptide(parsed)

        # Pick the three most frequent peptides so the comparison covers several PSMs
        targets = [sequence for sequence, _ in count_psms_per_peptide(parsed).most_common(3)]
        expected = [psm for psm in parsed if psm.peptidoform.sequence.upper() in targets]
        indexed = load_indexed_psms(index_files, targets)
        assert len(indexed) == len(expected) > 0
        assert [_sorted_proteins(p) for p in indexed] == [_sorted_proteins(p) for p in expected]

    def test_unchanged_file_is_not_parsed_again(self, tmp_path):
        ensure_idxml_index([SMALL_IDXML], tmp_path)
        assert missing_from_index([SMALL_IDXML], tmp_path) == []
        progress = []
        ensure_idxml_index([SMALL_IDXML], tmp_path,
                           progress_callback=lambda *args: progress.append(args))
        assert progress == []


class TestIndexKeysAndErrors:
    """Cache keys follow file changes; unreadable and broken files are reported every time."""

    def test_modified_file_gets_a_new_index_file(self, tmp_path):
        idxml = tmp_path / "a.idXML"
        idxml.write_text("x")
        before = index_file_for(str(idxml), tmp_path)
        os.utime(idxml, ns=(1, 1))
        assert index_file_for(str(idxml), tmp_path) != before

    def test_rebuild_after_a_change_removes_the_old_index(self, tmp_path):
        broken = tmp_path / "changing.idXML"
        broken.write_text("<not an idXML>")
        cache = tmp_path / "cache"
        first = ensure_idxml_index([str(broken)], cache)[0]
        os.utime(broken, ns=(1, 1))
        second = ensure_idxml_index([str(broken)], cache)[0]
        assert second != first and second.exists()
        assert sorted(path.name for path in cache.iterdir()) == sorted(
            [second.name, second.with_suffix(".errors.json").name])

    def test_parse_errors_are_reported_on_every_use(self, tmp_path):
        broken = tmp_path / "broken.idXML"
        broken.write_text("<not an idXML>")
        cache = tmp_path / "cache"
        for _ in range(2):
            errors = []
            index_files = ensure_idxml_index([str(broken)], cache, error_callback=errors.append)
            assert len(index_files) == 1
            assert errors and "broken.idXML" in errors[0]
        assert count_indexed_peptides(index_files) == {}
        assert len(load_indexed_psms(index_files, ["PEPTIDE"])) == 0

    def test_missing_file_is_skipped_with_a_message(self, tmp_path):
        errors = []
        index_files = ensure_idxml_index([str(tmp_path / "gone.idXML")], tmp_path,
                                         error_callback=errors.append)
        assert index_files == []
        assert errors and "gone.idXML" in errors[0]
        assert missing_from_index([str(tmp_path / "gone.idXML")], tmp_path) == [
            str(tmp_path / "gone.idXML")]

    def test_no_targets_or_no_files_returns_empty(self, tmp_path):
        assert count_indexed_peptides([]) == {}
        assert len(load_indexed_psms([], ["PEPTIDE"])) == 0
