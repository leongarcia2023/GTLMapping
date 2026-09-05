"""GTLMapping: spatially varying foreground extinction mapping."""

from .background import (
    estimate_box_background,
    estimate_lmf_background,
    estimate_smf_background,
    measure_box_background,
)
from .catalog import find_cloud, read_simon_catalog
from .exceptions import (
    CatalogError,
    GTLMappingError,
    GridMismatchError,
    InsufficientSamplesError,
)
from .extinction import (
    G_PER_CM2_TO_MSUN_PER_PC2,
    compute_extinction,
    convert_surface_density,
    propagate_uncertainty,
    transmission_std,
    unresolved_transmission,
)
from .foreground import (
    cross_validate_foreground,
    detect_saturated_samples,
    estimate_bt12_foreground,
    fit_conservative_foreground,
    fit_liberal_foreground,
    fit_moderate_foreground,
    interpolate_foreground,
)
from .mapper import GTLMapper
from .models import (
    BackgroundResult,
    CloudEllipse,
    ForegroundSamples,
    InterpolationResult,
    MappingResult,
    UncertaintyResult,
)
from .opacity import FilterOpacity, get_filter_opacity, list_filter_opacities

__all__ = [
    "BackgroundResult",
    "CatalogError",
    "CloudEllipse",
    "ForegroundSamples",
    "FilterOpacity",
    "GTLMapper",
    "GTLMappingError",
    "G_PER_CM2_TO_MSUN_PER_PC2",
    "GridMismatchError",
    "InsufficientSamplesError",
    "InterpolationResult",
    "MappingResult",
    "UncertaintyResult",
    "compute_extinction",
    "convert_surface_density",
    "cross_validate_foreground",
    "detect_saturated_samples",
    "estimate_bt12_foreground",
    "estimate_box_background",
    "estimate_lmf_background",
    "estimate_smf_background",
    "find_cloud",
    "fit_conservative_foreground",
    "fit_liberal_foreground",
    "fit_moderate_foreground",
    "get_filter_opacity",
    "interpolate_foreground",
    "list_filter_opacities",
    "measure_box_background",
    "propagate_uncertainty",
    "read_simon_catalog",
    "transmission_std",
    "unresolved_transmission",
]

__version__ = "0.5.0"
