"""
Unit tests for mspecviz/mirror_plot.py

Covers:
- create_plotly_mirror_plot / create_plotly_single_plot
- figure_to_bytes
"""
from io import BytesIO

import matplotlib
matplotlib.use("Agg")  # Must be before importing mirror_plot

import numpy as np
import plotly.graph_objects as go
import pytest
from matplotlib.figure import Figure

from mspecviz.annotation import prepare_spectrum_for_plot
from mspecviz.plotting import (
    create_mirror_plot,
    create_plotly_mirror_plot,
    create_plotly_single_plot,
    create_single_spectrum_plot,
    figure_to_bytes,
)


# ---------------------------------------------------------------------------
# TestPlotCreation (all 4 plot functions: returns correct type, handles empty)
# ---------------------------------------------------------------------------

class TestPlotCreation:
    """Parametrized tests for all plot creation functions."""

    @pytest.mark.parametrize("plot_fn, is_mirror", [
        (create_mirror_plot, True),
        (create_single_spectrum_plot, False),
        (create_plotly_mirror_plot, True),
        (create_plotly_single_plot, False),
    ])
    def test_returns_correct_figure_type(self, spectrum_frame, plot_fn, is_mirror):
        if is_mirror:
            fig = plot_fn(spectrum_frame, spectrum_frame.copy())
        else:
            fig = plot_fn(spectrum_frame)
        assert isinstance(fig, (Figure, go.Figure))

    @pytest.mark.parametrize("plot_fn, is_mirror", [
        (create_mirror_plot, True),
        (create_single_spectrum_plot, False),
        (create_plotly_mirror_plot, True),
        (create_plotly_single_plot, False),
    ])
    def test_empty_df_does_not_crash(self, empty_spectrum_frame, plot_fn, is_mirror):
        if is_mirror:
            fig = plot_fn(empty_spectrum_frame, empty_spectrum_frame.copy())
        else:
            fig = plot_fn(empty_spectrum_frame)
        assert isinstance(fig, (Figure, go.Figure))


# ---------------------------------------------------------------------------
# TestMirrorPlotExtras (Matplotlib-specific features)
# ---------------------------------------------------------------------------

class TestMirrorPlotExtras:
    def test_figsize_respected(self, spectrum_frame):
        fig = create_mirror_plot(spectrum_frame, spectrum_frame.copy(), figure_size=(10, 6))
        size = fig.get_size_inches()
        assert abs(size[0] - 10.0) < 0.1
        assert abs(size[1] - 6.0) < 0.1

    def test_metrics_text_accepted(self, spectrum_frame):
        metrics_text = "R\u00b2 match: 0.92\nCosine: 0.95"
        fig = create_mirror_plot(spectrum_frame, spectrum_frame.copy(), metrics_text=metrics_text)
        assert isinstance(fig, Figure)


# ---------------------------------------------------------------------------
# TestFigureSaveable
# ---------------------------------------------------------------------------

class TestFigureSaveable:
    def test_save_as_png(self, spectrum_frame):
        fig = create_mirror_plot(spectrum_frame, spectrum_frame.copy())
        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        buffer.seek(0)
        content = buffer.read()
        assert len(content) > 0
        # PNG files start with the PNG magic bytes
        assert content[:4] == b"\x89PNG"


# ---------------------------------------------------------------------------
# Integration tests with annotated spectra
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestMirrorPlotWithRealAnnotation:
    """Generate plots from annotated real-like spectra."""

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

        fig = create_plotly_mirror_plot(
            query, ref, title="Test Mirror",
            upper_label="Query", lower_label="Reference",
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

        fig = create_mirror_plot(df, df, title="Test")
        png_bytes = figure_to_bytes(fig, format="png")
        assert len(png_bytes) > 1000  # reasonable PNG size

        out = tmp_path / "test.png"
        out.write_bytes(png_bytes)
        assert out.exists()
