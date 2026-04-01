"""
Unit tests for src/mirror_plot.py

MSInsight changes from MSpecViz:
- create_interactive_mirror_plot -> create_plotly_mirror_plot
- create_interactive_single_plot -> create_plotly_single_plot
- Removed: create_legend_figure, add_color_column, figure_to_base64
- Added: figure_to_bytes
"""
import sys
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Must be before importing mirror_plot

import numpy as np
import pytest
from matplotlib.figure import Figure

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.mirror_plot import (
    ION_COLORS,
    create_mirror_plot,
    create_single_spectrum_plot,
    create_plotly_mirror_plot,
    create_plotly_single_plot,
    get_ion_color,
)


# ---------------------------------------------------------------------------
# TestGetIonColor
# ---------------------------------------------------------------------------

class TestGetIonColor:
    def test_unannotated_returns_gray(self):
        color = get_ion_color("unannotated")
        assert color == ION_COLORS["unannotated"]

    def test_b_ion_returns_color(self):
        color = get_ion_color("b3")
        assert isinstance(color, str)
        assert color.startswith("#")

    def test_y_ion_returns_color(self):
        color = get_ion_color("y2")
        assert isinstance(color, str)
        assert color.startswith("#")

    def test_unknown_returns_gray_fallback(self):
        # "hello_world" does not start with any ion letter (a,b,c,x,y,z) nor [M
        color = get_ion_color("hello_world")
        assert color == ION_COLORS["unannotated"]


# ---------------------------------------------------------------------------
# TestCreateMirrorPlot (Matplotlib static)
# ---------------------------------------------------------------------------

class TestCreateMirrorPlot:
    def test_returns_figure(self, spectrum_df):
        fig = create_mirror_plot(spectrum_df, spectrum_df.copy())
        assert isinstance(fig, Figure)

    def test_empty_dfs_do_not_crash(self, empty_spectrum_df):
        fig = create_mirror_plot(empty_spectrum_df, empty_spectrum_df.copy())
        assert isinstance(fig, Figure)

    def test_figsize_respected(self, spectrum_df):
        fig = create_mirror_plot(spectrum_df, spectrum_df.copy(), figsize=(10, 6))
        size = fig.get_size_inches()
        assert abs(size[0] - 10.0) < 0.1
        assert abs(size[1] - 6.0) < 0.1

    def test_metrics_dict_accepted(self, spectrum_df):
        metrics_text = "R\u00b2 match: 0.92\nCosine: 0.95"
        fig = create_mirror_plot(spectrum_df, spectrum_df.copy(), metrics_text=metrics_text)
        assert isinstance(fig, Figure)


# ---------------------------------------------------------------------------
# TestCreateSingleSpectrumPlot (Matplotlib static)
# ---------------------------------------------------------------------------

class TestCreateSingleSpectrumPlot:
    def test_returns_figure(self, spectrum_df):
        fig = create_single_spectrum_plot(spectrum_df)
        assert isinstance(fig, Figure)

    def test_empty_df_does_not_crash(self, empty_spectrum_df):
        fig = create_single_spectrum_plot(empty_spectrum_df)
        assert isinstance(fig, Figure)


# ---------------------------------------------------------------------------
# TestPlotlyMirrorPlot
# ---------------------------------------------------------------------------

class TestPlotlyMirrorPlot:
    def test_returns_plotly_figure(self, spectrum_df):
        import plotly.graph_objects as go
        fig = create_plotly_mirror_plot(spectrum_df, spectrum_df.copy())
        assert isinstance(fig, go.Figure)

    def test_empty_dfs_do_not_crash(self, empty_spectrum_df):
        import plotly.graph_objects as go
        fig = create_plotly_mirror_plot(empty_spectrum_df, empty_spectrum_df.copy())
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# TestPlotlySinglePlot
# ---------------------------------------------------------------------------

class TestPlotlySinglePlot:
    def test_returns_plotly_figure(self, spectrum_df):
        import plotly.graph_objects as go
        fig = create_plotly_single_plot(spectrum_df)
        assert isinstance(fig, go.Figure)

    def test_empty_df_does_not_crash(self, empty_spectrum_df):
        import plotly.graph_objects as go
        fig = create_plotly_single_plot(empty_spectrum_df)
        assert isinstance(fig, go.Figure)


# ---------------------------------------------------------------------------
# TestFigureSaveable
# ---------------------------------------------------------------------------

class TestFigureSaveable:
    def test_save_as_png(self, spectrum_df):
        fig = create_mirror_plot(spectrum_df, spectrum_df.copy())
        buf = BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        content = buf.read()
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
        from src.peak_annotation import prepare_spectrum_for_plot
        import plotly.graph_objects as go

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
        from src.peak_annotation import prepare_spectrum_for_plot
        from src.mirror_plot import figure_to_bytes

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
