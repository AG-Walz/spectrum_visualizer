"""Regression tests for gallery._show_plot image routing (mocked Streamlit)."""
from unittest.mock import patch

from streamlit_app.views.export.gallery import _show_plot


def test_show_plot_passes_svg_as_utf8_markup_not_path(tmp_path):
    # matplotlib SVGs carry non-ASCII glyphs (·, R², Å). Handing st.image the .svg PATH makes
    # Streamlit re-read the file with the platform-default (ascii) codec and crash with
    # UnicodeDecodeError. _show_plot must instead pass the decoded UTF-8 markup string so
    # Streamlit base64-embeds it directly, bypassing that file read.
    svg = tmp_path / "plot.svg"
    svg.write_text('<?xml version="1.0"?>\n<svg><text>R² · Å</text></svg>',
                   encoding="utf-8")
    assert any(byte > 127 for byte in svg.read_bytes())  # genuinely non-ASCII

    with patch("streamlit_app.views.export.gallery.st") as st_mock:
        _show_plot(svg)

    (passed,), _ = st_mock.image.call_args
    assert passed.lstrip().startswith(("<?xml", "<svg")), "SVG must be passed as markup, not a path"
    assert "R²" in passed, "SVG markup must be decoded as UTF-8"


def test_show_plot_passes_png_as_path(tmp_path):
    # Raster files are read as binary by Streamlit, so the path form is fine (and cheaper).
    png = tmp_path / "plot.png"
    png.write_bytes(b"\x89PNG\r\n\x1a\n")

    with patch("streamlit_app.views.export.gallery.st") as st_mock:
        _show_plot(png)

    (passed,), _ = st_mock.image.call_args
    assert passed == str(png), "PNG should be passed by path, not read into a string"
