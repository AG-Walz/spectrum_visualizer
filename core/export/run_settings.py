"""Single source of truth for run_settings.json (built by cli.main and streamlit_app.views.export).
Each frontend supplies a (key, default) -> value getter; the field list lives here so frontends cannot drift.
"""
from datetime import datetime
from typing import Any, Callable

from core.defaults import DEFAULTS

# A frontend's parameter source, adapted to a uniform (key, default) -> value
# lookup — the CLI passes getattr(config, ...), Streamlit reads params.
ParamGetter = Callable[[str, Any], Any]


# Each entry: (run_settings_key, default). Frontends map their internal storage
# (CLI SimpleNamespace attr, Streamlit params widget key) onto these keys
# via their getter — see cli.main and streamlit_app.utils.page_state.
RUN_SETTINGS_FIELDS: list[tuple[str, object]] = [
    ("idxml_files",            []),
    ("mzml_files",             []),
    ("reference_substring",    ""),
    ("ion_types",              list(DEFAULTS.ion_types)),
    ("mz_tolerance",           DEFAULTS.mz_tolerance),
    ("add_neutral_losses",     DEFAULTS.add_neutral_losses),
    ("intensity_mode",         DEFAULTS.intensity_mode),
    ("koina_enabled",          DEFAULTS.koina_enabled),
    ("koina_model",            DEFAULTS.koina_model),
    ("koina_collision_energy", DEFAULTS.koina_collision_energy),
    ("koina_instrument_type",  DEFAULTS.koina_instrument_type),
    ("export_formats",         list(DEFAULTS.export_formats)),
    ("export_resolution",      DEFAULTS.export_resolution),
    ("figure_size",            list(DEFAULTS.figure_size)),
    ("summary_metric",         DEFAULTS.summary_metric),
    ("summary_plot_types",     list(DEFAULTS.summary_plot_types)),
]


def build_run_settings(getter: ParamGetter, peptides: list[str]) -> dict:
    """Assemble the run_settings.json dict from any frontend's parameter source.
    peptides is passed separately because frontends parse it differently (YAML list vs newline string).
    The returned dict holds every RUN_SETTINGS_FIELDS key plus peptides and the
    underscore-prefixed timestamps _created_at / _last_run (both set to now on creation).
    """
    # Pull each declared field through the frontend getter, falling back to its default.
    settings: dict = {key: getter(key, default) for key, default in RUN_SETTINGS_FIELDS}
    settings["peptides"] = peptides
    # Both timestamps start equal; _last_run is later refreshed when an export reuses this folder.
    now = datetime.now().isoformat()
    settings["_created_at"] = now
    settings["_last_run"] = now
    return settings
