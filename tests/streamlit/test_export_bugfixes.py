"""
Regression tests for ``streamlit_app.views.export.page.export_controls_fragment``.

Each test here locks in a specific past bug so it can't come back. The main
export-pipeline tests live in ``test_export_pipeline.py``.
"""
import inspect

from streamlit_app.views.export.page import export_controls_fragment


class TestExportControlsBugFix:
    """Verify the NameError fix: export_dir is no longer referenced as undefined."""

    def test_export_controls_no_name_error(self):
        """The function source should reference get_param for export-dir, not bare export_dir."""
        source = inspect.getsource(export_controls_fragment)
        # The bug was: _run_export(export_dir, ...) where export_dir was undefined.
        # The fix uses get_param(params, 'export-dir', '') or export_fallback.
        assert "get_param" in source or "export_fallback" in source
        # Make sure the bare "export_dir," pattern on _run_export call is gone
        lines = source.split("\n")
        for line in lines:
            if "_run_export(" in line and "export_dir," in line:
                # Should NOT have bare export_dir as first arg
                assert "get_param" in line or "export_fallback" in line, \
                    f"export_dir should not be referenced directly: {line.strip()}"
