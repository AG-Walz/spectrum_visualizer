"""Source for all pipeline parameter defaults (Both cli and streamlit app read from here)."""

from dataclasses import dataclass
from typing import Tuple

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


DEFAULTS = Defaults()
