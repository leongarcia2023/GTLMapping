"""Package-specific exceptions."""


class GTLMappingError(Exception):
    """Base exception for GTLMapping."""


class CatalogError(GTLMappingError):
    """Raised when a catalog cannot be parsed or queried."""


class GridMismatchError(GTLMappingError):
    """Raised when arrays do not describe the same celestial grid."""


class InsufficientSamplesError(GTLMappingError):
    """Raised when interpolation has too few independent samples."""
