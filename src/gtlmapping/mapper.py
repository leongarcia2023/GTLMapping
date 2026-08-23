"""High-level stateful interface for an extinction-mapping workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .background import (
    estimate_box_background,
    estimate_lmf_background,
    estimate_smf_background,
)
from .catalog import find_cloud
from .extinction import compute_extinction, propagate_uncertainty
from .foreground import (
    detect_saturated_samples,
    estimate_bt12_foreground,
    fit_conservative_foreground,
    fit_liberal_foreground,
    fit_moderate_foreground,
    interpolate_foreground,
)
from .geometry import ellipse_mask, validate_compatible_grids
from .io import read_fits_image, reproject_to_grid
from .models import (
    BackgroundResult,
    CloudEllipse,
    ForegroundSamples,
    InterpolationResult,
    MappingResult,
)
from .opacity import get_filter_opacity


class GTLMapper:
    """Build a GTL extinction map one inspectable stage at a time."""

    def __init__(
        self,
        observed: np.ndarray,
        *,
        header: fits.Header | None = None,
        wcs: WCS | None = None,
        observed_std: np.ndarray | None = None,
    ) -> None:
        data = np.asarray(observed, dtype=float)
        if data.ndim != 2:
            raise ValueError("observed must be a two-dimensional image.")
        self.observed = data
        self.header = fits.Header() if header is None else header.copy()
        self.wcs = WCS(self.header) if wcs is None else wcs
        if observed_std is None:
            self.observed_std: np.ndarray | None = None
        else:
            uncertainty = np.asarray(observed_std, dtype=float)
            if uncertainty.shape != data.shape:
                raise ValueError("observed_std must match observed.shape.")
            if np.any(np.isfinite(uncertainty) & (uncertainty < 0)):
                raise ValueError("observed_std must be non-negative.")
            self.observed_std = uncertainty
        self.cloud: CloudEllipse | None = None
        self.samples: ForegroundSamples | None = None
        self.foreground_result: InterpolationResult | None = None
        self.background_result: BackgroundResult | None = None
        self.background_std: np.ndarray | None = None

    @classmethod
    def from_fits(
        cls,
        path: str | Path,
        *,
        hdu: int | str = 0,
        uncertainty_hdu: int | str | None = None,
    ) -> "GTLMapper":
        """Create a mapper from a FITS image and optional uncertainty HDU."""

        data, header, wcs = read_fits_image(path, hdu=hdu)
        uncertainty = None
        if uncertainty_hdu is not None:
            uncertainty, _, uncertainty_wcs = read_fits_image(
                path,
                hdu=uncertainty_hdu,
            )
            if uncertainty.shape != data.shape:
                raise ValueError(
                    "The uncertainty image shape does not match the science "
                    f"image: {uncertainty.shape} versus {data.shape}."
                )
            if uncertainty_wcs.has_celestial:
                validate_compatible_grids(
                    data.shape,
                    wcs,
                    uncertainty.shape,
                    uncertainty_wcs,
                )
        return cls(
            data,
            header=header,
            wcs=wcs,
            observed_std=uncertainty,
        )

    def set_observed_uncertainty(self, uncertainty: np.ndarray) -> np.ndarray:
        """Attach an observed-intensity standard-deviation array."""

        values = np.asarray(uncertainty, dtype=float)
        if values.shape != self.observed.shape:
            raise ValueError("uncertainty must match observed.shape.")
        if np.any(np.isfinite(values) & (values < 0)):
            raise ValueError("uncertainty must be non-negative.")
        self.observed_std = values
        return values

    def select_cloud(
        self,
        catalog: str | Path,
        name: str,
        *,
        component: str = "0",
    ) -> CloudEllipse:
        """Select a cloud/core directly from ``catalog.dat``."""

        self.cloud = find_cloud(catalog, name, component=component)
        return self.cloud

    def detect_foreground(
        self,
        *,
        cloud: CloudEllipse | None = None,
        region_mask: np.ndarray | None = None,
        restrict_to_cloud: bool = False,
        **kwargs: Any,
    ) -> ForegroundSamples:
        """Detect and store local GTL foreground samples.

        The notebook's GTL search spans the full image, which is therefore
        the default. Set ``restrict_to_cloud=True`` (or pass ``region_mask``)
        to limit the search spatially.
        """

        selected = cloud or self.cloud
        if region_mask is None and restrict_to_cloud:
            if selected is None:
                raise RuntimeError(
                    "restrict_to_cloud=True requires a selected cloud."
                )
            region_mask = ellipse_mask(self.observed.shape, self.wcs, selected)
        self.samples = detect_saturated_samples(
            self.observed,
            wcs=self.wcs,
            region_mask=region_mask,
            **kwargs,
        )
        return self.samples

    def fit_foreground(
        self,
        *,
        method: str = "conservative",
        samples: ForegroundSamples | None = None,
        cloud: CloudEllipse | None = None,
        region_mask: np.ndarray | None = None,
        **kwargs: Any,
    ) -> InterpolationResult:
        """Fit and store a conservative/liberal GTL, legacy, or BT12 foreground."""

        normalized = method.lower()
        if normalized == "bt12":
            selected = cloud or self.cloud
            if region_mask is None:
                if selected is None:
                    raise RuntimeError(
                        "BT12 foreground estimation requires a selected cloud "
                        "or region_mask."
                    )
                region_mask = ellipse_mask(
                    self.observed.shape,
                    self.wcs,
                    selected,
                )
            self.foreground_result = estimate_bt12_foreground(
                self.observed,
                wcs=self.wcs,
                region_mask=region_mask,
                **kwargs,
            )
            return self.foreground_result

        selected = samples or self.samples
        if selected is None:
            raise RuntimeError("Call detect_foreground() before fit_foreground().")
        if normalized == "conservative":
            selected_cloud = cloud or self.cloud
            if region_mask is None and selected_cloud is not None:
                region_mask = ellipse_mask(
                    self.observed.shape,
                    self.wcs,
                    selected_cloud,
                )
            self.foreground_result = fit_conservative_foreground(
                selected,
                self.observed,
                wcs=self.wcs,
                region_mask=region_mask,
                **kwargs,
            )
            return self.foreground_result
        if normalized == "liberal":
            selected_cloud = cloud or self.cloud
            if region_mask is None and selected_cloud is not None:
                region_mask = ellipse_mask(
                    self.observed.shape,
                    self.wcs,
                    selected_cloud,
                )
            self.foreground_result = fit_liberal_foreground(
                selected,
                self.observed,
                wcs=self.wcs,
                region_mask=region_mask,
                **kwargs,
            )
            return self.foreground_result
        if normalized == "moderate":
            selected_cloud = cloud or self.cloud
            if region_mask is None and selected_cloud is not None:
                region_mask = ellipse_mask(
                    self.observed.shape,
                    self.wcs,
                    selected_cloud,
                )
            self.foreground_result = fit_moderate_foreground(
                selected,
                self.observed,
                wcs=self.wcs,
                region_mask=region_mask,
                **kwargs,
            )
            return self.foreground_result
        self.foreground_result = interpolate_foreground(
            selected,
            self.observed.shape,
            method=normalized,
            **kwargs,
        )
        return self.foreground_result

    def estimate_background(
        self,
        *,
        method: str = "smf",
        cloud: CloudEllipse | None = None,
        **kwargs: Any,
    ) -> BackgroundResult:
        """Estimate and store an LMF, SMF, or adjacent-box background."""

        normalized = method.lower()
        if normalized == "lmf":
            result = estimate_lmf_background(self.observed, self.wcs, **kwargs)
        elif normalized == "smf":
            selected = cloud or self.cloud
            if selected is None:
                raise RuntimeError("SMF estimation requires a selected cloud.")
            result = estimate_smf_background(
                self.observed,
                self.wcs,
                selected,
                **kwargs,
            )
        elif normalized == "boxes":
            result = estimate_box_background(self.observed, **kwargs)
        else:
            raise ValueError(
                "Background method must be 'smf', 'lmf', or 'boxes'."
            )
        self.background_result = result
        self.background_std = None
        return result

    def set_background(
        self,
        background: np.ndarray,
        *,
        method: str = "provided",
        uncertainty: np.ndarray | float | None = None,
    ) -> BackgroundResult:
        """Use an already aligned background array."""

        values = np.asarray(background, dtype=float)
        if values.shape != self.observed.shape:
            raise ValueError(
                f"Background shape {values.shape} does not match "
                f"observed shape {self.observed.shape}."
            )
        if uncertainty is None:
            self.background_std = None
        else:
            standard_deviation = np.broadcast_to(
                np.asarray(uncertainty, dtype=float),
                values.shape,
            ).copy()
            if np.any(
                np.isfinite(standard_deviation) & (standard_deviation < 0)
            ):
                raise ValueError("Background uncertainty must be non-negative.")
            self.background_std = standard_deviation
        self.background_result = BackgroundResult(values=values, method=method)
        return self.background_result

    def set_background_from_fits(
        self,
        path: str | Path,
        *,
        hdu: int | str = 0,
        align: bool = False,
        tolerance_pixels: float = 0.1,
        uncertainty_hdu: int | str | None = None,
    ) -> BackgroundResult:
        """Load a background FITS file, rejecting WCS mismatch by default."""

        values, _, background_wcs = read_fits_image(path, hdu=hdu)
        if align:
            values = reproject_to_grid(
                values,
                background_wcs,
                self.wcs,
                self.observed.shape,
            )
            diagnostics = {"reprojected": True}
        else:
            diagnostics = validate_compatible_grids(
                self.observed.shape,
                self.wcs,
                values.shape,
                background_wcs,
                tolerance_pixels=tolerance_pixels,
            )
            diagnostics["reprojected"] = False
        if uncertainty_hdu is None:
            self.background_std = None
        else:
            uncertainty, _, uncertainty_wcs = read_fits_image(
                path,
                hdu=uncertainty_hdu,
            )
            if align:
                uncertainty = reproject_to_grid(
                    uncertainty,
                    uncertainty_wcs,
                    self.wcs,
                    self.observed.shape,
                )
            else:
                validate_compatible_grids(
                    self.observed.shape,
                    self.wcs,
                    uncertainty.shape,
                    uncertainty_wcs,
                    tolerance_pixels=tolerance_pixels,
                )
            if np.any(np.isfinite(uncertainty) & (uncertainty < 0)):
                raise ValueError("Background uncertainty must be non-negative.")
            self.background_std = np.asarray(uncertainty, dtype=float)
        self.background_result = BackgroundResult(
            values=values,
            method="provided",
            diagnostics=diagnostics,
        )
        return self.background_result

    def constrain_foreground(
        self,
        *,
        minimum_transmitted_intensity: float | np.ndarray,
        minimum_foreground: float | None = None,
    ) -> InterpolationResult:
        """Project the foreground onto the radiatively feasible range.

        The observed off-cloud background contains both foreground and
        transmitted background emission, so physically
        ``I_fore <= I_bg - minimum_transmitted_intensity``. Kriging and a
        separately estimated background can violate this inequality. This
        explicit projection removes those model-conflict NaNs while retaining
        a pixel mask and adjustment diagnostics in the result.
        """

        if self.foreground_result is None:
            raise RuntimeError("No foreground is available.")
        if self.background_result is None:
            raise RuntimeError("No background is available.")
        transmission = np.broadcast_to(
            np.asarray(minimum_transmitted_intensity, dtype=float),
            self.observed.shape,
        )
        if np.any(~np.isfinite(transmission) | (transmission <= 0)):
            raise ValueError(
                "minimum_transmitted_intensity must be positive and finite."
            )
        ceiling = self.background_result.values - transmission
        selected_floor = minimum_foreground
        if (
            selected_floor is None
            and self.foreground_result.diagnostics.get(
                "bt12_floor_enforced",
                False,
            )
        ):
            selected_floor = float(
                self.foreground_result.diagnostics["reference_foreground"]
            )
        if selected_floor is not None:
            floor = float(selected_floor)
            if np.any(np.isfinite(ceiling) & (ceiling < floor)):
                raise ValueError(
                    "The requested transmitted-intensity floor is incompatible "
                    "with the background and minimum_foreground."
                )
        original = np.asarray(self.foreground_result.values, dtype=float)
        constraint_mask = (
            np.isfinite(original)
            & np.isfinite(ceiling)
            & (original > ceiling)
        )
        adjusted = np.minimum(original, ceiling)
        if selected_floor is not None:
            adjusted = np.maximum(adjusted, float(selected_floor))
        difference = original - adjusted
        diagnostics = dict(self.foreground_result.diagnostics)
        diagnostics.update(
            {
                "foreground_constraint": "background_minus_transmission",
                "preserved_foreground_floor": selected_floor,
                "constrained_pixel_count": int(np.count_nonzero(constraint_mask)),
                "constrained_pixel_fraction": float(np.mean(constraint_mask)),
                "maximum_constraint_adjustment": float(
                    np.nanmax(np.where(constraint_mask, difference, 0.0))
                ),
            }
        )
        self.foreground_result = InterpolationResult(
            values=adjusted,
            method=self.foreground_result.method,
            variance=self.foreground_result.variance,
            constraint_mask=constraint_mask,
            diagnostics=diagnostics,
        )
        return self.foreground_result

    def compute(
        self,
        *,
        kappa_cm2_g: float | None = None,
        filter_name: str | None = None,
        dust_model: str = "oh94_thin_ice_coagulated",
        gas_to_dust_ratio: float = 156.0,
        kappa_std_cm2_g: float = 0.0,
        bright_pixel_policy: str = "allow",
        saturation_policy: str = "mask",
        intensity_floor: float | np.ndarray | None = None,
        observed_std: float | np.ndarray | None = None,
        background_std: float | np.ndarray | None = None,
        foreground_std: float | np.ndarray | None = None,
        use_kriging_variance: bool = True,
    ) -> MappingResult:
        """Compute a mapping result from the stored foreground/background."""

        if self.foreground_result is None:
            raise RuntimeError("No foreground is available.")
        if self.background_result is None:
            raise RuntimeError("No background is available.")
        if filter_name is not None and kappa_cm2_g is not None:
            raise ValueError("Specify either filter_name or kappa_cm2_g, not both.")
        if filter_name is not None:
            selected_kappa = get_filter_opacity(
                filter_name,
                dust_model=dust_model,
                gas_to_dust_ratio=gas_to_dust_ratio,
            )
        else:
            selected_kappa = 7.5 if kappa_cm2_g is None else float(kappa_cm2_g)
        tau, sigma, saturated, invalid_background, bright = compute_extinction(
            self.observed,
            self.background_result.values,
            self.foreground_result.values,
            kappa_cm2_g=selected_kappa,
            bright_pixel_policy=bright_pixel_policy,
            saturation_policy=saturation_policy,
            intensity_floor=intensity_floor,
        )
        selected_observed_std = (
            self.observed_std if observed_std is None else observed_std
        )
        selected_background_std = (
            self.background_std if background_std is None else background_std
        )
        selected_foreground_std = foreground_std
        if (
            selected_foreground_std is None
            and use_kriging_variance
            and self.foreground_result.variance is not None
        ):
            selected_foreground_std = np.sqrt(
                np.maximum(self.foreground_result.variance, 0.0)
            )
        has_uncertainty = (
            selected_observed_std is not None
            or selected_background_std is not None
            or selected_foreground_std is not None
            or kappa_std_cm2_g > 0
        )
        uncertainty = None
        if has_uncertainty:
            uncertainty = propagate_uncertainty(
                self.observed,
                self.background_result.values,
                self.foreground_result.values,
                observed_std=(
                    0.0
                    if selected_observed_std is None
                    else selected_observed_std
                ),
                background_std=(
                    0.0
                    if selected_background_std is None
                    else selected_background_std
                ),
                foreground_std=(
                    0.0
                    if selected_foreground_std is None
                    else selected_foreground_std
                ),
                kappa_cm2_g=selected_kappa,
                kappa_std_cm2_g=kappa_std_cm2_g,
                additional_mask=saturated | invalid_background,
            )
        metadata = {
            "fg_method": self.foreground_result.method,
            "bg_method": self.background_result.method,
            "n_samples": len(self.samples) if self.samples is not None else 0,
            "saturation_policy": saturation_policy,
        }
        if filter_name is not None:
            metadata.update(
                {
                    "filter": filter_name,
                    "dust_model": dust_model,
                    "gas_dust": float(gas_to_dust_ratio),
                }
            )
        return MappingResult(
            surface_density=sigma,
            optical_depth=tau,
            foreground=self.foreground_result.values,
            background=self.background_result.values,
            saturated_mask=saturated,
            invalid_background_mask=invalid_background,
            bright_mask=bright,
            header=self.header,
            wcs=self.wcs,
            kappa_cm2_g=selected_kappa,
            foreground_constraint_mask=self.foreground_result.constraint_mask,
            uncertainty=uncertainty,
            metadata=metadata,
        )

    def compute_liberal(
        self,
        *,
        intensity_floor: float | np.ndarray | None = None,
        minimum_foreground: float | None = None,
        **kwargs: Any,
    ) -> MappingResult:
        """Compute a finite liberal map with censored pixels as lower limits.

        The liberal foreground intentionally permits a controlled number of
        pixels with ``I_obs <= I_fg``. This convenience method projects that
        foreground below ``I_bg - intensity_floor`` and computes those pixels
        with ``saturation_policy='lower_limit'``. The ``SATURATED`` extension
        remains the authoritative lower-limit mask; finite values there are
        not ordinary detections.

        If ``intensity_floor`` is omitted, the liberal fit's documented
        recommendation (normally ``2 * noise_sigma``) is used. A physically
        impossible background, non-finite input, or an explicitly conflicting
        keyword still raises instead of silently manufacturing a value.
        """

        if self.foreground_result is None:
            raise RuntimeError("No foreground is available.")
        if self.foreground_result.method != "liberal":
            raise RuntimeError(
                "compute_liberal() requires fit_foreground(method='liberal')."
            )
        if self.background_result is None:
            raise RuntimeError("No background is available.")
        if "saturation_policy" in kwargs:
            raise ValueError(
                "compute_liberal() always uses saturation_policy='lower_limit'."
            )
        if "intensity_floor" in kwargs:
            raise ValueError(
                "Pass intensity_floor directly to compute_liberal(), not in kwargs."
            )

        selected_floor: float | np.ndarray
        if intensity_floor is None:
            recommended = self.foreground_result.diagnostics.get(
                "recommended_intensity_floor"
            )
            if recommended is None:
                raise ValueError(
                    "No recommended intensity floor is stored; pass one explicitly."
                )
            selected_floor = float(recommended)
        else:
            selected_floor = intensity_floor

        self.constrain_foreground(
            minimum_transmitted_intensity=selected_floor,
            minimum_foreground=minimum_foreground,
        )
        return self.compute(
            saturation_policy="lower_limit",
            intensity_floor=selected_floor,
            **kwargs,
        )

    def compute_moderate(
        self,
        *,
        intensity_floor: float | np.ndarray | None = None,
        minimum_foreground: float | None = None,
        **kwargs: Any,
    ) -> MappingResult:
        """Compute a finite moderate map with flagged lower limits.

        This is the moderate counterpart to :meth:`compute_liberal`: it
        projects the foreground below the background by the transmitted-
        intensity floor and keeps censored pixels finite but marked in
        ``SATURATED``.
        """

        if self.foreground_result is None:
            raise RuntimeError("No foreground is available.")
        if self.foreground_result.method != "moderate":
            raise RuntimeError(
                "compute_moderate() requires fit_foreground(method='moderate')."
            )
        if self.background_result is None:
            raise RuntimeError("No background is available.")
        if "saturation_policy" in kwargs:
            raise ValueError(
                "compute_moderate() always uses "
                "saturation_policy='lower_limit'."
            )
        if "intensity_floor" in kwargs:
            raise ValueError(
                "Pass intensity_floor directly to compute_moderate(), not in kwargs."
            )

        selected_floor: float | np.ndarray
        if intensity_floor is None:
            recommended = self.foreground_result.diagnostics.get(
                "recommended_intensity_floor"
            )
            if recommended is None:
                raise ValueError(
                    "No recommended intensity floor is stored; pass one explicitly."
                )
            selected_floor = float(recommended)
        else:
            selected_floor = intensity_floor

        self.constrain_foreground(
            minimum_transmitted_intensity=selected_floor,
            minimum_foreground=minimum_foreground,
        )
        return self.compute(
            saturation_policy="lower_limit",
            intensity_floor=selected_floor,
            **kwargs,
        )
