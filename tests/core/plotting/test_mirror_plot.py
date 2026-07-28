"""Unit + integration tests for core.plotting.mirror (pyopenms-viz wrapper)."""
from io import BytesIO

import matplotlib
matplotlib.use("Agg")  # MUST come before importing plotting modules (headless backend).

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from matplotlib.figure import Figure

from core.peak_annotation import prepare_spectrum_for_plot
from core.plotting import create_spectrum_plot, group_mz_axis_max
from core.plotting.mirror import _MZ_PADDING_DA, _prep


class TestMaxPeakLabels:
    """The max_peak_labels cap keeps only the N most intense annotated peaks labelled (issue #93)."""

    def _annotated_frame(self, n: int) -> pd.DataFrame:
        # Ascending intensity so the top-N by intensity are the last N rows.
        return pd.DataFrame({
            "mz": np.linspace(100, 1000, n),
            "intensity": np.arange(1, n + 1, dtype=float),
            "annotation": [f"b{i}" for i in range(n)],
            "ion_type": ["b"] * n,
        })

    def test_zero_keeps_all_labels(self):
        prepped = _prep(self._annotated_frame(10), max_peak_labels=0)
        assert (prepped["_label"] != "").sum() == 10

    def test_cap_keeps_only_top_n_by_intensity(self):
        prepped = _prep(self._annotated_frame(10), max_peak_labels=3)
        kept = set(prepped.loc[prepped["_label"] != "", "annotation"])
        # Highest three intensities are the last three rows (b7, b8, b9).
        assert kept == {"b7", "b8", "b9"}

    def test_cap_larger_than_labelled_keeps_all(self):
        prepped = _prep(self._annotated_frame(4), max_peak_labels=99)
        assert (prepped["_label"] != "").sum() == 4

    def test_cap_never_drops_peaks_only_labels(self):
        # Every peak stem must survive; only the text labels are thinned.
        prepped = _prep(self._annotated_frame(10), max_peak_labels=2)
        assert len(prepped) == 10

    def test_empty_frame_is_safe(self):
        prepped = _prep(pd.DataFrame({"mz": [], "intensity": []}), max_peak_labels=5)
        assert prepped.empty


class TestPlotCreation:
    """Return-type + empty-frame safety across ms_matplotlib / ms_plotly and single/mirror modes."""

    @pytest.mark.parametrize("backend, expected_type", [
        ("ms_matplotlib", Figure),
        ("ms_plotly", go.Figure),
    ])
    @pytest.mark.parametrize("with_reference", [True, False])
    def test_returns_correct_figure_type(self, spectrum_frame, backend, expected_type, with_reference):
        reference = spectrum_frame.copy() if with_reference else None
        fig = create_spectrum_plot(spectrum_frame, reference, backend=backend)
        assert isinstance(fig, expected_type)

    @pytest.mark.parametrize("backend", ["ms_matplotlib", "ms_plotly"])
    @pytest.mark.parametrize("with_reference", [True, False])
    def test_empty_df_does_not_crash(self, empty_spectrum_frame, backend, with_reference):
        reference = empty_spectrum_frame.copy() if with_reference else None
        fig = create_spectrum_plot(empty_spectrum_frame, reference, backend=backend)
        assert isinstance(fig, (Figure, go.Figure))


class TestXAxisStartsAtZero:
    """Every plot's m/z axis is anchored at 0 so PSMs share one origin (issue #114)."""

    def _frame(self, first_mz: float) -> pd.DataFrame:
        return pd.DataFrame({
            "mz": [first_mz, first_mz + 300.0],
            "intensity": [100.0, 50.0],
            "annotation": ["b1", "y1"],
            "ion_type": ["b", "y"],
        })

    @pytest.mark.parametrize("first_mz", [150.0, 800.0])
    @pytest.mark.parametrize("with_reference", [True, False])
    def test_matplotlib_xlim_starts_at_zero(self, first_mz, with_reference):
        frame = self._frame(first_mz)
        fig = create_spectrum_plot(
            frame, frame.copy() if with_reference else None, backend="ms_matplotlib",
        )
        left, right = fig.axes[0].get_xlim()
        assert left == 0.0
        # The upper bound still follows the spectrum, padded so outer peaks aren't clipped.
        assert right > first_mz + 300.0

    @pytest.mark.parametrize("first_mz", [150.0, 800.0])
    @pytest.mark.parametrize("with_reference", [True, False])
    def test_plotly_x_range_starts_at_zero(self, first_mz, with_reference):
        frame = self._frame(first_mz)
        fig = create_spectrum_plot(
            frame, frame.copy() if with_reference else None, backend="ms_plotly",
        )
        low, high = fig.layout.xaxis.range
        assert low == 0.0
        assert high > first_mz + 300.0

    def test_two_spectra_share_the_same_origin(self):
        # The regression: different first peaks used to give different axis starts.
        low_fig = create_spectrum_plot(self._frame(150.0), backend="ms_matplotlib")
        high_fig = create_spectrum_plot(self._frame(800.0), backend="ms_matplotlib")
        assert low_fig.axes[0].get_xlim()[0] == high_fig.axes[0].get_xlim()[0]


class TestGroupWideXAxis:
    """A peptide group's plots all end at the group's highest peak (issue #114)."""

    def _frame(self, last_mz: float) -> pd.DataFrame:
        return pd.DataFrame({
            "mz": [120.0, last_mz], "intensity": [100.0, 50.0],
            "annotation": ["b1", "y1"], "ion_type": ["b", "y"],
        })

    def test_group_max_pads_the_highest_peak_across_frames(self):
        bound = group_mz_axis_max([self._frame(400.0), self._frame(1500.0), None])
        assert bound == 1500.0 + _MZ_PADDING_DA

    def test_group_max_ignores_empty_frames(self):
        empty = pd.DataFrame({"mz": [], "intensity": []})
        assert group_mz_axis_max([empty, self._frame(900.0)]) == 900.0 + _MZ_PADDING_DA

    def test_group_max_of_nothing_is_zero(self):
        # 0.0 is the "auto" sentinel create_spectrum_plot falls back on.
        assert group_mz_axis_max([]) == 0.0

    @pytest.mark.parametrize("backend", ["ms_matplotlib", "ms_plotly"])
    def test_every_plot_in_a_group_gets_the_same_axis(self, backend):
        # The short spectrum must be drawn on the long one's axis, not its own.
        frames = [self._frame(400.0), self._frame(1500.0)]
        bound = group_mz_axis_max(frames)
        ranges = []
        for frame in frames:
            fig = create_spectrum_plot(frame, backend=backend, mz_axis_max=bound)
            ranges.append(tuple(fig.axes[0].get_xlim()) if backend == "ms_matplotlib"
                          else tuple(fig.layout.xaxis.range))
        assert ranges[0] == ranges[1] == (0.0, bound)

    def test_zero_falls_back_to_the_plots_own_peaks(self):
        fig = create_spectrum_plot(self._frame(400.0), backend="ms_matplotlib", mz_axis_max=0.0)
        assert fig.axes[0].get_xlim() == (0.0, 400.0 + _MZ_PADDING_DA)


class TestMatplotlibExtras:
    def test_metrics_text_accepted(self, spectrum_frame):
        metrics_text = "R² match: 0.92\nCosine: 0.95"
        fig = create_spectrum_plot(
            spectrum_frame, spectrum_frame.copy(),
            backend="ms_matplotlib", metrics_text=metrics_text,
        )
        assert isinstance(fig, Figure)


class TestFigureSaveable:
    def test_save_as_png(self, spectrum_frame):
        fig = create_spectrum_plot(
            spectrum_frame, spectrum_frame.copy(), backend="ms_matplotlib",
        )
        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        content = buffer.getvalue()
        assert len(content) > 0
        assert content[:4] == b"\x89PNG"


@pytest.mark.integration
class TestMirrorPlotWithRealAnnotation:
    """End-to-end: annotate realistic spectra then render via each backend."""

    def test_plotly_mirror_with_annotation(self):
        mz = np.array([175.12, 290.14, 405.17, 520.20], dtype=np.float64)
        intensity = np.array([100, 200, 150, 300], dtype=np.float64)

        query = prepare_spectrum_for_plot(
            mz, intensity, "ASDDFDQASY", charge=1,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
        )
        ref = prepare_spectrum_for_plot(
            mz * 1.001, intensity * 0.9, "ASDDFDQASY", charge=1,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
        )

        fig = create_spectrum_plot(
            query, ref, backend="ms_plotly", title="Test Mirror",
            query_label="Query", reference_label="Reference",
        )
        assert isinstance(fig, go.Figure)
        assert len(fig.data) > 0

    def test_matplotlib_mirror_with_annotation(self, tmp_path):
        mz = np.array([175.12, 290.14, 405.17], dtype=np.float64)
        intensity = np.array([100, 200, 150], dtype=np.float64)

        df = prepare_spectrum_for_plot(
            mz, intensity, "RTEFNLNQY", charge=2,
            ion_types=["b", "y"], mz_tolerance=0.02, normalize=True,
        )

        fig = create_spectrum_plot(df, df, backend="ms_matplotlib", title="Test")
        out = tmp_path / "test.png"
        fig.savefig(str(out), format="png", bbox_inches="tight", facecolor="white")
        assert out.exists()
        # >1 kB ⇒ Matplotlib actually rendered something (not a blank page)
        assert out.stat().st_size > 1000
