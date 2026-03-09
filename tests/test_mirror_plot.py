"""
Unit tests for src/mirror_plot.py
"""
import sys
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Must be before importing mirror_plot

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.mirror_plot import (
    ION_COLORS,
    create_mirror_plot,
    create_single_spectrum_plot,
    get_ion_color,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_spectrum_df(n_peaks=5):
    """Build a minimal annotated spectrum DataFrame."""
    mz = np.linspace(100.0, 600.0, n_peaks)
    intensity = np.linspace(0.2, 1.0, n_peaks)
    annotations = ["b1", "y2", "b3", "y4", "unannotated"][:n_peaks]
    is_annotated = [a != "unannotated" for a in annotations]
    return pd.DataFrame({
        "mz": mz,
        "intensity": intensity,
        "annotation": annotations,
        "is_annotated": is_annotated,
    })


def _empty_spectrum_df():
    return pd.DataFrame({
        "mz": [],
        "intensity": [],
        "annotation": [],
        "is_annotated": [],
    })


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
# TestCreateMirrorPlot
# ---------------------------------------------------------------------------

class TestCreateMirrorPlot:
    def test_returns_figure(self):
        upper = _make_spectrum_df()
        lower = _make_spectrum_df()
        fig = create_mirror_plot(upper, lower)
        assert isinstance(fig, Figure)

    def test_empty_dfs_do_not_crash(self):
        upper = _empty_spectrum_df()
        lower = _empty_spectrum_df()
        fig = create_mirror_plot(upper, lower)
        assert isinstance(fig, Figure)

    def test_figsize_respected(self):
        upper = _make_spectrum_df()
        lower = _make_spectrum_df()
        fig = create_mirror_plot(upper, lower, figsize=(10, 6))
        size = fig.get_size_inches()
        assert abs(size[0] - 10.0) < 0.1
        assert abs(size[1] - 6.0) < 0.1

    def test_metrics_dict_accepted(self):
        upper = _make_spectrum_df()
        lower = _make_spectrum_df()
        metrics_text = "R² match: 0.92\nCosine: 0.95"
        fig = create_mirror_plot(upper, lower, metrics_text=metrics_text)
        assert isinstance(fig, Figure)


# ---------------------------------------------------------------------------
# TestCreateSingleSpectrumPlot
# ---------------------------------------------------------------------------

class TestCreateSingleSpectrumPlot:
    def test_returns_figure(self):
        spectrum = _make_spectrum_df()
        fig = create_single_spectrum_plot(spectrum)
        assert isinstance(fig, Figure)

    def test_empty_df_does_not_crash(self):
        spectrum = _empty_spectrum_df()
        fig = create_single_spectrum_plot(spectrum)
        assert isinstance(fig, Figure)


# ---------------------------------------------------------------------------
# TestFigureSaveable
# ---------------------------------------------------------------------------

class TestFigureSaveable:
    def test_save_as_png(self):
        upper = _make_spectrum_df()
        lower = _make_spectrum_df()
        fig = create_mirror_plot(upper, lower)
        buf = BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        content = buf.read()
        assert len(content) > 0
        # PNG files start with the PNG magic bytes
        assert content[:4] == b"\x89PNG"
