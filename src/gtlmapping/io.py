"""FITS input and optional reprojection helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def read_fits_image(
    path: str | Path,
    *,
    hdu: int | str = 0,
) -> tuple[np.ndarray, fits.Header, WCS]:
    """Read a finite-or-NaN two-dimensional FITS image."""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    with fits.open(source, memmap=False) as hdul:
        selected = hdul[hdu]
        if selected.data is None:
            raise ValueError(f"HDU {hdu!r} in {source} has no image data.")
        data = np.asarray(selected.data, dtype=float)
        header = selected.header.copy()
    if data.ndim != 2:
        raise ValueError(
            f"HDU {hdu!r} in {source} is {data.ndim}D; a 2D image is required."
        )
    return data, header, WCS(header)


def reproject_to_grid(
    data: np.ndarray,
    source_wcs: WCS,
    target_wcs: WCS,
    target_shape: tuple[int, int],
) -> np.ndarray:
    """Reproject an array to a target celestial grid using interpolation."""

    try:
        from reproject import reproject_interp
    except ImportError as exc:
        raise ImportError(
            "Automatic alignment requires the optional 'align' extra: "
            "pip install 'GTLMapping[align]'."
        ) from exc
    aligned, footprint = reproject_interp(
        (np.asarray(data, dtype=float), source_wcs.celestial),
        target_wcs.celestial,
        shape_out=target_shape,
    )
    aligned = np.asarray(aligned, dtype=float)
    aligned[np.asarray(footprint) <= 0] = np.nan
    return aligned
