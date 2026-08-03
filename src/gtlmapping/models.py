"""Typed result containers used throughout GTLMapping."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


@dataclass(frozen=True, slots=True)
class CloudEllipse:
    """One cloud or core entry from the Simon et al. (2006) catalog."""

    name: str
    component: str
    glon_deg: float
    glat_deg: float
    major_axis_arcmin: float
    minor_axis_arcmin: float
    pa_deg: float
    area_arcmin2: float | None = None
    peak_contrast: float | None = None
    contrast_snr: float | None = None
    integrated_contrast_arcmin2: float | None = None

    @property
    def is_cloud(self) -> bool:
        """Whether this entry is the cloud rather than a cataloged core."""

        return self.component == "0"


@dataclass(frozen=True, slots=True)
class ForegroundSamples:
    """Validated foreground sample locations and intensities."""

    rows: np.ndarray
    cols: np.ndarray
    values: np.ndarray
    multiplicity: np.ndarray
    accepted_windows: int
    rejected_windows: int
    total_windows: int
    pixel_scale_arcsec: float
    min_separation_arcsec: float

    def __post_init__(self) -> None:
        lengths = {
            len(self.rows),
            len(self.cols),
            len(self.values),
            len(self.multiplicity),
        }
        if len(lengths) != 1:
            raise ValueError("Sample arrays must have matching lengths.")
        if np.any(np.asarray(self.multiplicity) < 1):
            raise ValueError("Sample multiplicities must be positive integers.")

    def __len__(self) -> int:
        return len(self.values)

    @property
    def raw_detection_count(self) -> int:
        """Number of accepted windows before duplicate coordinates are merged."""

        return int(np.sum(self.multiplicity))

    @property
    def points_xy(self) -> np.ndarray:
        """Unique sample coordinates in ``(x, y) = (column, row)`` order."""

        return np.column_stack((self.cols, self.rows))


@dataclass(frozen=True, slots=True)
class InterpolationResult:
    """A modeled foreground and optional model variance."""

    values: np.ndarray
    method: str
    variance: np.ndarray | None = None
    constraint_mask: np.ndarray | None = None
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class BackgroundResult:
    """A background map and information about its construction."""

    values: np.ndarray
    method: str
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class UncertaintyResult:
    """First-order uncertainty products and their surface-density components."""

    optical_depth_std: np.ma.MaskedArray
    surface_density_std: np.ma.MaskedArray
    components: dict[str, np.ma.MaskedArray] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class MappingResult:
    """Extinction products and masks on a shared celestial grid."""

    surface_density: np.ma.MaskedArray
    optical_depth: np.ma.MaskedArray
    foreground: np.ndarray
    background: np.ndarray
    saturated_mask: np.ndarray
    invalid_background_mask: np.ndarray
    bright_mask: np.ndarray
    header: fits.Header
    wcs: WCS
    kappa_cm2_g: float
    foreground_constraint_mask: np.ndarray | None = None
    uncertainty: UncertaintyResult | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def write(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Write a transparent multi-extension FITS product.

        The primary HDU is mass surface density. Additional extensions retain
        optical depth, the foreground/background models, and diagnostic masks.
        """

        output = Path(path)
        header = self.header.copy()
        header["BUNIT"] = ("g cm-2", "Mass surface density")
        header["EXTNAME"] = "SIGMA"
        header["HIERARCH GTL KAPPA"] = (
            float(self.kappa_cm2_g),
            "Opacity in cm2 g-1",
        )
        for key, value in self.metadata.items():
            if isinstance(value, (str, int, float, bool)) and len(str(key)) <= 48:
                header[f"HIERARCH GTL {str(key).upper()}"] = value
        header.add_history("Created by GTLMapping.")

        hdus = [
            fits.PrimaryHDU(self.surface_density.filled(np.nan), header=header),
            fits.ImageHDU(
                self.optical_depth.filled(np.nan),
                header=self._image_header("TAU", ""),
                name="TAU",
            ),
            fits.ImageHDU(
                np.asarray(self.foreground, dtype=float),
                header=self._image_header("FOREGROUND", "MJy sr-1"),
                name="FOREGROUND",
            ),
            fits.ImageHDU(
                np.asarray(self.background, dtype=float),
                header=self._image_header("BACKGROUND", "MJy sr-1"),
                name="BACKGROUND",
            ),
            fits.ImageHDU(
                np.asarray(self.saturated_mask, dtype=np.uint8),
                header=self._image_header("SATURATED", "bool"),
                name="SATURATED",
            ),
            fits.ImageHDU(
                np.asarray(self.invalid_background_mask, dtype=np.uint8),
                header=self._image_header("INVALID_BG", "bool"),
                name="INVALID_BG",
            ),
            fits.ImageHDU(
                np.asarray(self.bright_mask, dtype=np.uint8),
                header=self._image_header("BRIGHT", "bool"),
                name="BRIGHT",
            ),
        ]
        if self.uncertainty is not None:
            hdus.extend(
                [
                    fits.ImageHDU(
                        self.uncertainty.surface_density_std.filled(np.nan),
                        header=self._image_header("SIGMA_ERR", "g cm-2"),
                        name="SIGMA_ERR",
                    ),
                    fits.ImageHDU(
                        self.uncertainty.optical_depth_std.filled(np.nan),
                        header=self._image_header("TAU_ERR", ""),
                        name="TAU_ERR",
                    ),
                ]
            )
        if self.foreground_constraint_mask is not None:
            hdus.append(
                fits.ImageHDU(
                    np.asarray(self.foreground_constraint_mask, dtype=np.uint8),
                    header=self._image_header("FG_CONSTRAINT", "bool"),
                    name="FG_CONSTRAINT",
                )
            )
        fits.HDUList(hdus).writeto(output, overwrite=overwrite)
        return output

    def _image_header(self, name: str, unit: str) -> fits.Header:
        header = self.wcs.to_header()
        header["EXTNAME"] = name
        if unit:
            header["BUNIT"] = unit
        return header
