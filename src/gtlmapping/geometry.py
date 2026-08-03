"""Celestial-grid and cloud-geometry utilities."""

from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from .exceptions import GridMismatchError
from .models import CloudEllipse


def pixel_scale_arcsec(wcs: WCS, *, fallback: float | None = None) -> float:
    """Return the mean projected celestial pixel scale in arcseconds."""

    try:
        celestial = wcs.celestial
        scales = np.abs(proj_plane_pixel_scales(celestial)) * 3600.0
        scale = float(np.mean(scales))
    except Exception as exc:
        if fallback is None:
            raise ValueError("A valid celestial WCS is required.") from exc
        return float(fallback)
    if not np.isfinite(scale) or scale <= 0:
        if fallback is None:
            raise ValueError("WCS has no positive finite pixel scale.")
        return float(fallback)
    return scale


def ellipse_mask(
    shape: tuple[int, int],
    wcs: WCS,
    cloud: CloudEllipse,
) -> np.ndarray:
    """Create a WCS-aware Simon ellipse mask.

    Position angle is interpreted east of north, independent of image-axis
    parity. This avoids the sign/rotation ambiguity of adding 90 degrees in
    pixel coordinates.
    """

    if len(shape) != 2:
        raise ValueError("Ellipse masks require a two-dimensional shape.")
    yy, xx = np.indices(shape, dtype=float)
    sky = wcs.celestial.pixel_to_world(xx, yy)
    center = SkyCoord(
        l=cloud.glon_deg * u.deg,
        b=cloud.glat_deg * u.deg,
        frame="galactic",
    )
    east, north = center.spherical_offsets_to(sky.galactic)
    east_arcmin = east.to_value(u.arcmin)
    north_arcmin = north.to_value(u.arcmin)

    theta = np.deg2rad(cloud.pa_deg)
    along_major = east_arcmin * np.sin(theta) + north_arcmin * np.cos(theta)
    along_minor = east_arcmin * np.cos(theta) - north_arcmin * np.sin(theta)
    semi_major = cloud.major_axis_arcmin / 2.0
    semi_minor = cloud.minor_axis_arcmin / 2.0
    if semi_major <= 0 or semi_minor <= 0:
        raise ValueError("Cloud axes must be positive.")
    return (along_major / semi_major) ** 2 + (along_minor / semi_minor) ** 2 <= 1


def validate_compatible_grids(
    shape_a: tuple[int, int],
    wcs_a: WCS,
    shape_b: tuple[int, int],
    wcs_b: WCS,
    *,
    tolerance_pixels: float = 0.1,
) -> dict[str, float]:
    """Validate that two arrays share shape and celestial pixel centers."""

    if tuple(shape_a) != tuple(shape_b):
        raise GridMismatchError(
            f"Array shapes differ: {tuple(shape_a)} versus {tuple(shape_b)}."
        )
    if tolerance_pixels < 0:
        raise ValueError("tolerance_pixels must be non-negative.")

    ny, nx = shape_a
    x = np.array([0.0, nx - 1.0, 0.0, nx - 1.0, (nx - 1.0) / 2.0])
    y = np.array([0.0, 0.0, ny - 1.0, ny - 1.0, (ny - 1.0) / 2.0])
    sky_a = wcs_a.celestial.pixel_to_world(x, y)
    sky_b = wcs_b.celestial.pixel_to_world(x, y)
    separations = sky_a.separation(sky_b).to_value(u.arcsec)
    scale = min(pixel_scale_arcsec(wcs_a), pixel_scale_arcsec(wcs_b))
    max_pixels = float(np.nanmax(separations) / scale)
    diagnostics = {
        "max_separation_arcsec": float(np.nanmax(separations)),
        "max_separation_pixels": max_pixels,
    }
    if not np.all(np.isfinite(separations)) or max_pixels > tolerance_pixels:
        raise GridMismatchError(
            "Celestial grids are misaligned by up to "
            f"{diagnostics['max_separation_arcsec']:.3f} arcsec "
            f"({max_pixels:.3f} pixels); allowed {tolerance_pixels:.3f} pixels."
        )
    return diagnostics
