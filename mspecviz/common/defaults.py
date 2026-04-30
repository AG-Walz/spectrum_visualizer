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
    export_resolution: int = 300               # raster resolution (png/pdf)
    figure_size: Tuple[int, int] = (12, 8)     # [width, height] in inches
    plot_unpaired_queries: bool = False        # render queries without a reference

    # --- Summary plots ------------------------------------------------------
    summary_metric: str = "R² binned"
    summary_plot_types: Tuple[str, ...] = ("dual_axis", "ridgeline")

    def __post_init__(self):
        """Fail loud at construction time if a config supplies an unknown enum value."""
        if self.intensity_mode not in _INTENSITY_MODES:
            raise ValueError(
                f"intensity_mode={self.intensity_mode!r} not one of {sorted(_INTENSITY_MODES)}"
            )
        if self.koina_instrument_type not in _KOINA_INSTRUMENTS:
            raise ValueError(
                f"koina_instrument_type={self.koina_instrument_type!r} "
                f"not one of {sorted(_KOINA_INSTRUMENTS)}"
            )


_INTENSITY_MODES = frozenset({"Relative", "Raw"})
_KOINA_INSTRUMENTS = frozenset({"QE", "LUMOS", "TIMSTOF", "ELITE", "SCIEXTOF"})


DEFAULTS = Defaults()
