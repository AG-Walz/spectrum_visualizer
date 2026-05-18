"""Pack annotation + Koina settings into kwargs dicts (called by ``cli.main`` and ``streamlit_app.views.search_load``)."""
from typing import Any, Callable

from core.common.defaults import DEFAULTS

ParamGetter = Callable[[str, Any], Any]


def build_annotation_kwargs(getter: ParamGetter) -> dict:
    """Read ion types, tolerance, intensity mode, and neutral-loss flag via ``getter`` for ``prepare_spectrum_for_plot``."""
    return dict(
        ion_types=getter("ion_types", DEFAULTS.ion_types),
        mz_tolerance=getter("mz_tolerance", DEFAULTS.mz_tolerance),
        normalize=getter("intensity_mode", DEFAULTS.intensity_mode) == "Relative",
        add_neutral_losses=getter("add_neutral_losses", DEFAULTS.add_neutral_losses),
    )


def build_koina_kwargs(getter: ParamGetter) -> dict:
    """Read Koina enable flag, model, collision energy, and instrument via ``getter`` for ``process_query_result``."""
    return dict(
        koina_enabled=getter("koina_enabled", DEFAULTS.koina_enabled),
        koina_model=getter("koina_model", DEFAULTS.koina_model),
        koina_collision_energy=getter("koina_collision_energy", DEFAULTS.koina_collision_energy),
        koina_instrument_type=getter("koina_instrument_type", DEFAULTS.koina_instrument_type),
    )
