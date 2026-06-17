"""Unit + integration tests for core.plotting.mirror (pyopenms-viz wrapper)."""
from io import BytesIO

import matplotlib
matplotlib.use("Agg")  # MUST come before importing plotting modules (headless backend).

import numpy as np
import plotly.graph_objects as go
import pytest
from matplotlib.figure import Figure

from core.peak_annotation import prepare_spectrum_for_plot
from core.plotting import create_spectrum_plot


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
