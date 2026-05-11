"""Frontend-agnostic kwargs builders for the streaming-annotation pipeline.

Both the CLI (`cli/main.py`) and the Streamlit search-and-load view need to
pack annotation + Koina settings into the kwargs dicts consumed by
``process_query_result``. The two frontends store the same values under
different shapes (CLI uses a ``SimpleNamespace``; Streamlit uses session
``params``). The functions here take a *getter* callable
``(key, default) -> value`` so each frontend can supply its own lookup.
"""
from typing import Any, Callable

from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.helpers import is_relative_intensity

ParamGetter = Callable[[str, Any], Any]


def build_annotation_kwargs(getter: ParamGetter) -> dict:
    """Pack annotation params into the kwargs dict for ``prepare_spectrum_for_plot``."""
    return dict(
        ion_types=getter("ion_types", DEFAULTS.ion_types),
        mz_tolerance=getter("mz_tolerance", DEFAULTS.mz_tolerance),
        normalize=is_relative_intensity(getter("intensity_mode", DEFAULTS.intensity_mode)),
        add_neutral_losses=getter("add_neutral_losses", DEFAULTS.add_neutral_losses),
    )


def build_koina_kwargs(getter: ParamGetter) -> dict:
    """Pack Koina settings into the kwargs dict for ``process_query_result``."""
    return dict(
        koina_enabled=getter("koina_enabled", DEFAULTS.koina_enabled),
        koina_model=getter("koina_model", DEFAULTS.koina_model),
        koina_collision_energy=getter("koina_collision_energy", DEFAULTS.koina_collision_energy),
        koina_instrument_type=getter("koina_instrument_type", DEFAULTS.koina_instrument_type),
    )
