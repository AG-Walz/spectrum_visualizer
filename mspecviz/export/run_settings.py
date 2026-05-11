"""Single source of truth for ``run_settings.json`` field list and assembly.

Both the CLI (``cli/main.py``) and the Streamlit Export page need to capture
the user's parameters into ``run_settings.json`` for reproducibility. The
field list lives here so the two frontends cannot drift; each frontend supplies
a *getter* callable ``(key, default) -> value`` that knows how to read from its
local config representation.
"""
from datetime import datetime
from typing import Any, Callable

from mspecviz.common.defaults import DEFAULTS

ParamGetter = Callable[[str, Any], Any]


# Each entry: (run_settings_key, default). Frontends map their internal storage
# (CLI ``SimpleNamespace`` attribute, Streamlit ``params`` widget key) onto these
# run_settings keys via their getter — see ``cli/main.py`` and
# ``streamlit_app/utils/page_state.py``.
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
    """Assemble the ``run_settings.json`` dict from any frontend's parameter source.

    *peptides* is supplied separately because the two frontends parse it
    differently (CLI: list from YAML / file; Streamlit: newline-joined string
    in session ``target_peptides``). Timestamps are filled in at call time.
    """
    settings: dict = {key: getter(key, default) for key, default in RUN_SETTINGS_FIELDS}
    settings["peptides"] = peptides
    now = datetime.now().isoformat()
    settings["_created_at"] = now
    settings["_last_run"] = now
    return settings
