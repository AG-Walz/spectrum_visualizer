"""Single source of truth for pipeline parameter defaults (``mspecviz``/``cli``/``streamlit_app`` all read here).

Precedence: this file < settings.json/params.json/YAML < CLI flags /
Streamlit widgets. All keys are ``snake_case`` everywhere.
"""
from dataclasses import dataclass, field, fields
from typing import List, Tuple


@dataclass(frozen=True)
class Defaults:
    """Canonical defaults. Instantiated once as ``DEFAULTS`` below."""

    # --- Annotation ---------------------------------------------------------
    ion_types: Tuple[str, ...] = ("b", "y")
    mz_tolerance: float = 0.02                 # Da
    add_neutral_losses: bool = False           # include H2O, NH3 fragment losses
    add_first_prefix_ion: bool = True          # include a1/b1/c1 ions
    remove_precursor_peaks: bool = False       # mask peaks near precursor m/z
    intensity_mode: str = "Relative"           # "Relative" | "Raw"

    # --- Performance --------------------------------------------------------
    workers: int = 4                           # 1 = sequential fallback

    # --- Koina AI prediction ------------------------------------------------
    koina_enabled: bool = False                # off by default: requires network
    koina_model: str = "ms2pip_timsTOF2024"
    koina_collision_energy: int = 35           # eV
    koina_instrument_type: str = "TIMSTOF"

    # --- Export -------------------------------------------------------------
    export_directory: str = "./exports"
    project_name: str = "MSpecViz_analysis"
    export_formats: Tuple[str, ...] = ("png",)  # CLI batch list: any of png/svg/pdf/html
    image_format: str = "png"                  # Streamlit single-image picker default
    export_resolution: int = 300                      # raster resolution (png/pdf)
    figure_size: Tuple[int, int] = (12, 8)     # [width, height] in inches
    plot_unpaired_queries: bool = False        # render queries without a reference

    # --- Summary plots ------------------------------------------------------
    summary_metric: str = "R\u00b2 binned"
    summary_plot_types: Tuple[str, ...] = ("dual_axis", "ridgeline")

    def as_dict(self) -> dict:
        """Return a plain dict — tuples stay tuples; copy lists at use sites."""
        return {f.name: getattr(self, f.name) for f in fields(self)}


DEFAULTS = Defaults()


# Module-level constants kept for convenience at import sites that previously
# declared their own local DEFAULT_* (peak_theory, koina_prediction,
# _summary_shared). Re-exporting here lets those modules import a single
# name without dragging in the whole dataclass.
DEFAULT_ION_TYPES: List[str] = list(DEFAULTS.ion_types)
DEFAULT_KOINA_MODEL: str = DEFAULTS.koina_model
DEFAULT_KOINA_COLLISION_ENERGY: int = DEFAULTS.koina_collision_energy
DEFAULT_KOINA_INSTRUMENT_TYPE: str = DEFAULTS.koina_instrument_type
DEFAULT_SUMMARY_METRIC: str = DEFAULTS.summary_metric
DEFAULT_SUMMARY_PLOT_TYPES: List[str] = list(DEFAULTS.summary_plot_types)
