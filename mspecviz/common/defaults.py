"""Single source of truth for user-tunable pipeline parameters."""

from dataclasses import dataclass


# Parameter option catalogues (imported by the rest of the codebase)

KOINA_MODELS: dict[str, str] = {
    # MS2PIP family
    "MS2PIP TimsTOF 2024": "ms2pip_timsTOF2024",
    "MS2PIP HCD 2021": "ms2pip_HCD2021",
    "MS2PIP Immuno HCD": "ms2pip_Immuno_HCD",
    "MS2PIP TTOF5600": "ms2pip_TTOF5600",
    "MS2PIP CID TMT": "ms2pip_CID_TMT",
    "MS2PIP iTRAQ Phospho": "ms2pip_iTRAQphospho",
    # Prosit family
    "Prosit 2019": "Prosit_2019_intensity",
    "Prosit 2020 HCD": "Prosit_2020_intensity_HCD",
    "Prosit 2020 CID": "Prosit_2020_intensity_CID",
    "Prosit 2020 TMT": "Prosit_2020_intensity_TMT",
    "Prosit 2023 timsTOF": "Prosit_2023_intensity_timsTOF",

    # Other
    "AlphaPept Generic": "AlphaPept_ms2_generic",
    "UniSpec": "UniSpec",
    "PredFull": "PredFull",
}

INSTRUMENT_TYPES: tuple[str, ...] = ("QE", "LUMOS", "TIMSTOF", "SCIEXTOF")

# Extended spectrum_utils' ion palette 
ION_COLORS: dict[str, str] = {
    "a":           "#388E3C",  # green       (Material Green 700)
    "b":           "#1976D2",  # blue        (Material Blue 700)
    "c":           "#00796B",  # teal        (Material Teal 700)
    "x":           "#7B1FA2",  # purple      (Material Purple 700)
    "y":           "#D32F2F",  # red         (Material Red 700)
    "z":           "#F57C00",  # orange      (Material Orange 700)
    "internal":    "#FBC02D",  # yellow      (Material Yellow 700) — spectrum_utils' "m"
    "precursor":   "#512DA8",  # deep purple (Material Deep Purple 700)
    "immonium":    "#455A64",  # blue-grey   (Material Blue Grey 700)
    "unannotated": "#9E9E9E",  # grey        (Material Grey 500)
}

IMAGE_FORMATS: tuple[str, ...] = ("png", "svg", "pdf", "html")

# Defaults for all user-tunable parameters, with validation at construction time

@dataclass(frozen=True)
class Defaults:
    """Canonical defaults. Instantiated once as ``DEFAULTS`` below."""

    # --- Annotation ---------------------------------------------------------
    ion_types: tuple[str, ...] = ("b", "y")
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
    export_formats: tuple[str, ...] = ("png",)  # CLI batch list: any of IMAGE_FORMATS
    image_format: str = "png"                  # Streamlit single-image picker default
    export_resolution: int = 300               # raster resolution (png/pdf)
    figure_size: tuple[int, int] = (12, 8)     # [width, height] in inches
    plot_unpaired_queries: bool = False        # render queries without a reference

    # --- Summary plots ------------------------------------------------------
    summary_metric: str = "R² binned"
    summary_plot_types: tuple[str, ...] = ("dual_axis", "ridgeline")

    def __post_init__(self):
        """Fail loud at construction time if a config supplies an unknown enum value."""
        for field, allowed in (
            ("intensity_mode", _INTENSITY_MODES),
            ("koina_instrument_type", _KOINA_INSTRUMENTS),
            ("koina_model", _KOINA_MODEL_IDS),
            ("image_format", _IMAGE_FORMATS),
        ):
            value = getattr(self, field)
            if value not in allowed:
                raise ValueError(f"{field}={value!r} not one of {sorted(allowed)}")
        bad = set(self.export_formats) - _IMAGE_FORMATS
        if bad:
            raise ValueError(
                f"export_formats contains {sorted(bad)!r}; each must be one of {sorted(_IMAGE_FORMATS)}"
            )

# Dictionaries and sets for validating config values at construction time
_INTENSITY_MODES = frozenset({"Relative", "Raw"})
_KOINA_INSTRUMENTS = frozenset(INSTRUMENT_TYPES)
_KOINA_MODEL_IDS = frozenset(KOINA_MODELS.values())
_IMAGE_FORMATS = frozenset(IMAGE_FORMATS)


DEFAULTS = Defaults()