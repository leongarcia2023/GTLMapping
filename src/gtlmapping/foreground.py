"""Detection and interpolation of spatially varying foreground emission."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np
from astropy.wcs import WCS
from scipy.interpolate import RBFInterpolator, SmoothBivariateSpline, griddata
from scipy.optimize import least_squares

from .exceptions import InsufficientSamplesError
from .geometry import pixel_scale_arcsec
from .models import ForegroundSamples, InterpolationResult


def _scan_starts(
    length: int,
    box_size: int,
    stride: int,
    *,
    cover_edges: bool,
) -> np.ndarray:
    if box_size < 1 or box_size > length:
        raise ValueError("box_size must be between 1 and the image dimension.")
    if stride < 1:
        raise ValueError("stride must be at least one pixel.")
    starts = list(range(0, length - box_size + 1, stride))
    last = length - box_size
    if cover_edges and (not starts or starts[-1] != last):
        starts.append(last)
    return np.asarray(starts, dtype=int)


def detect_saturated_samples(
    image: np.ndarray,
    *,
    wcs: WCS | None = None,
    region_mask: np.ndarray | None = None,
    excluded_mask: np.ndarray | None = None,
    grid_n: int = 5,
    overlap: float = 0.5,
    noise_sigma: float = 0.6,
    beam_fwhm_arcsec: float = 2.0,
    min_separation_arcsec: float | None = None,
    fallback_pixel_scale_arcsec: float = 1.2,
    cover_edges: bool = False,
) -> ForegroundSamples:
    """Find locally saturated minima using overlapping aliasing windows.

    A minimum is accepted when a distinct pixel in the same window lies
    within ``2 * noise_sigma`` and is sufficiently far away. Duplicate
    coordinates from overlapping windows are merged and their multiplicity
    retained. By default, the scan origins exactly match the prototype
    notebook. Set ``cover_edges=True`` to append a final window when the
    stride does not land on the last valid origin.
    """

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be a two-dimensional array.")
    if grid_n < 1:
        raise ValueError("grid_n must be at least one.")
    if not 0 <= overlap < 1:
        raise ValueError("overlap must satisfy 0 <= overlap < 1.")
    if noise_sigma <= 0 or beam_fwhm_arcsec <= 0:
        raise ValueError("Noise and beam FWHM must be positive.")

    mask = (
        np.ones(data.shape, dtype=bool)
        if region_mask is None
        else np.asarray(region_mask, dtype=bool)
    )
    if mask.shape != data.shape:
        raise ValueError("region_mask must match image.shape.")
    if excluded_mask is not None:
        excluded = np.asarray(excluded_mask, dtype=bool)
        if excluded.shape != data.shape:
            raise ValueError("excluded_mask must match image.shape.")
        mask = mask & ~excluded
    if not np.any(mask & np.isfinite(data)):
        raise ValueError("No finite image pixels are inside region_mask.")

    pixscale = (
        pixel_scale_arcsec(wcs, fallback=fallback_pixel_scale_arcsec)
        if wcs is not None
        else float(fallback_pixel_scale_arcsec)
    )
    separation = (
        float(min_separation_arcsec)
        if min_separation_arcsec is not None
        else 4.0 * (beam_fwhm_arcsec / 2.0)
    )
    if separation <= 0:
        raise ValueError("min_separation_arcsec must be positive.")
    separation_pix = separation / pixscale

    ny, nx = data.shape
    box_y = max(1, ny // grid_n)
    box_x = max(1, nx // grid_n)
    stride_y = max(1, int(np.floor(box_y * (1.0 - overlap))))
    stride_x = max(1, int(np.floor(box_x * (1.0 - overlap))))
    y_starts = _scan_starts(
        ny,
        box_y,
        stride_y,
        cover_edges=cover_edges,
    )
    x_starts = _scan_starts(
        nx,
        box_x,
        stride_x,
        cover_edges=cover_edges,
    )

    detections: list[tuple[int, int, float]] = []
    accepted = 0
    rejected = 0
    for y0 in y_starts:
        for x0 in x_starts:
            window = data[y0 : y0 + box_y, x0 : x0 + box_x]
            window_mask = mask[y0 : y0 + box_y, x0 : x0 + box_x]
            eligible = window_mask & np.isfinite(window)
            if not np.any(eligible):
                rejected += 1
                continue
            eligible_values = np.where(eligible, window, np.nan)
            flat_index = int(np.nanargmin(eligible_values))
            row_rel, col_rel = np.unravel_index(flat_index, window.shape)
            min_value = float(window[row_rel, col_rel])

            candidate = eligible & (window <= min_value + 2.0 * noise_sigma)
            rows, cols = np.indices(window.shape)
            distance = np.hypot(rows - row_rel, cols - col_rel)
            has_partner = bool(np.any(candidate & (distance >= separation_pix)))
            if has_partner:
                detections.append((y0 + row_rel, x0 + col_rel, min_value))
                accepted += 1
            else:
                rejected += 1

    grouped: dict[tuple[int, int], list[float]] = defaultdict(list)
    for row, col, value in detections:
        grouped[(row, col)].append(value)
    coordinates = sorted(grouped)
    rows = np.asarray([row for row, _ in coordinates], dtype=int)
    cols = np.asarray([col for _, col in coordinates], dtype=int)
    values = np.asarray(
        [np.mean(grouped[coordinate]) for coordinate in coordinates], dtype=float
    )
    multiplicity = np.asarray(
        [len(grouped[coordinate]) for coordinate in coordinates], dtype=int
    )
    return ForegroundSamples(
        rows=rows,
        cols=cols,
        values=values,
        multiplicity=multiplicity,
        accepted_windows=accepted,
        rejected_windows=rejected,
        total_windows=len(y_starts) * len(x_starts),
        pixel_scale_arcsec=pixscale,
        min_separation_arcsec=separation,
    )


def estimate_bt12_foreground(
    image: np.ndarray,
    *,
    wcs: WCS | None = None,
    region_mask: np.ndarray | None = None,
    noise_sigma: float = 0.6,
    min_separation_arcsec: float = 8.0,
    foreground_margin: float | None = None,
    fallback_pixel_scale_arcsec: float = 1.2,
) -> InterpolationResult:
    """Estimate the spatially constant foreground prescribed by BT12.

    The global minimum is measured inside ``region_mask``. All pixels
    strictly between that minimum and ``minimum + 2 * noise_sigma`` are
    labeled saturated if at least one such pixel is spatially independent
    of the minimum. The foreground is their mean intensity minus
    ``2 * noise_sigma`` (or an explicitly supplied ``foreground_margin``).
    """

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be a two-dimensional array.")
    if noise_sigma <= 0:
        raise ValueError("noise_sigma must be positive.")
    if min_separation_arcsec <= 0:
        raise ValueError("min_separation_arcsec must be positive.")
    margin = (
        2.0 * float(noise_sigma)
        if foreground_margin is None
        else float(foreground_margin)
    )
    if margin < 0:
        raise ValueError("foreground_margin must be non-negative.")

    mask = np.ones(data.shape, dtype=bool) if region_mask is None else np.asarray(
        region_mask, dtype=bool
    )
    if mask.shape != data.shape:
        raise ValueError("region_mask must match image.shape.")
    eligible = mask & np.isfinite(data)
    if not np.any(eligible):
        raise ValueError("No finite image pixels are inside region_mask.")

    pixscale = (
        pixel_scale_arcsec(wcs, fallback=fallback_pixel_scale_arcsec)
        if wcs is not None
        else float(fallback_pixel_scale_arcsec)
    )
    minimum_flat = int(np.nanargmin(np.where(eligible, data, np.nan)))
    minimum_row, minimum_col = np.unravel_index(minimum_flat, data.shape)
    minimum = float(data[minimum_row, minimum_col])
    saturated = (
        eligible
        & (data > minimum)
        & (data < minimum + 2.0 * float(noise_sigma))
    )

    yy, xx = np.indices(data.shape)
    separation_pixels = float(min_separation_arcsec) / pixscale
    distance = np.hypot(yy - minimum_row, xx - minimum_col)
    independent = saturated & (distance >= separation_pixels)
    independent_count = int(np.count_nonzero(independent))
    if independent_count == 0:
        raise InsufficientSamplesError(
            "BT12 foreground estimation found no saturated pixel at least "
            f"{min_separation_arcsec:g} arcsec from the global minimum."
        )

    saturated_values = data[saturated]
    saturated_mean = float(np.mean(saturated_values))
    foreground_value = saturated_mean - margin
    return InterpolationResult(
        values=np.full(data.shape, foreground_value, dtype=float),
        method="bt12",
        diagnostics={
            "global_minimum": minimum,
            "global_minimum_row": int(minimum_row),
            "global_minimum_col": int(minimum_col),
            "saturated_pixel_count": int(saturated_values.size),
            "independent_saturated_pixel_count": independent_count,
            "saturated_mean": saturated_mean,
            "noise_sigma": float(noise_sigma),
            "foreground_margin": margin,
            "pixel_scale_arcsec": pixscale,
            "min_separation_arcsec": float(min_separation_arcsec),
        },
    )


def fit_conservative_foreground(
    samples: ForegroundSamples,
    image: np.ndarray,
    *,
    wcs: WCS | None = None,
    region_mask: np.ndarray | None = None,
    noise_sigma: float = 0.6,
    min_separation_arcsec: float | None = None,
    foreground_margin: float | None = None,
    robust_loss: str = "soft_l1",
    maximum_local_saturation_fraction: float = 0.01,
    maximum_strict_saturation_fraction: float = 0.0,
    blend_steps: int = 101,
    floor: float | None = 0.0,
    fallback_pixel_scale_arcsec: float = 1.2,
) -> InterpolationResult:
    """Fit a robust spatial foreground with a hard BT12 lower bound.

    Local-window minima do not by themselves prove saturation: on structured
    clouds they can encode the cloud morphology. This conservative model fits
    only a broad plane to the local minima, subtracts the same foreground
    margin used by BT12, and then applies the historical GTL rule
    ``foreground = maximum(spatial_foreground, BT12_foreground)``. Thus the
    spatial model can add foreground relative to BT12 but cannot remove it.
    A data-driven blend factor limits both pixels within ``2 * noise_sigma``
    of the foreground and pixels strictly below it.

    The result is deliberately a trend rather than an exact interpolation
    through every minimum. Raw kriging remains available explicitly through
    :func:`interpolate_foreground`, but is not a safe default for extinction
    mapping without independent validation of every anchor. By default, the
    model is not allowed to introduce any new strictly saturated pixels in the
    fitted region. This preserves the pointwise BT12 mass lower bound wherever
    the BT12 radiative-transfer solution is valid.
    """

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be a two-dimensional array.")
    if noise_sigma <= 0 or not np.isfinite(noise_sigma):
        raise ValueError("noise_sigma must be positive and finite.")
    if robust_loss not in {"linear", "soft_l1", "huber", "cauchy", "arctan"}:
        raise ValueError(
            "robust_loss must be linear, soft_l1, huber, cauchy, or arctan."
        )
    for name, fraction in (
        ("maximum_local_saturation_fraction", maximum_local_saturation_fraction),
        ("maximum_strict_saturation_fraction", maximum_strict_saturation_fraction),
    ):
        if not 0 <= fraction <= 1:
            raise ValueError(f"{name} must lie between zero and one.")
    if blend_steps < 2:
        raise ValueError("blend_steps must be at least two.")
    _validate_samples(samples, "kriging")

    mask = (
        np.ones(data.shape, dtype=bool)
        if region_mask is None
        else np.asarray(region_mask, dtype=bool)
    )
    if mask.shape != data.shape:
        raise ValueError("region_mask must match image.shape.")
    eligible = mask & np.isfinite(data)
    if not np.any(eligible):
        raise ValueError("No finite image pixels are inside region_mask.")

    selected_separation = (
        samples.min_separation_arcsec
        if min_separation_arcsec is None
        else float(min_separation_arcsec)
    )
    reference = estimate_bt12_foreground(
        data,
        wcs=wcs,
        region_mask=mask,
        noise_sigma=noise_sigma,
        min_separation_arcsec=selected_separation,
        foreground_margin=foreground_margin,
        fallback_pixel_scale_arcsec=fallback_pixel_scale_arcsec,
    )

    ny, nx = data.shape
    x_center = 0.5 * (nx - 1)
    y_center = 0.5 * (ny - 1)
    x_scale = max(x_center, 1.0)
    y_scale = max(y_center, 1.0)
    sample_design = np.column_stack(
        (
            np.ones(len(samples), dtype=float),
            (samples.cols.astype(float) - x_center) / x_scale,
            (samples.rows.astype(float) - y_center) / y_scale,
        )
    )
    weights = np.sqrt(samples.multiplicity.astype(float))
    initial = np.array(
        [float(np.median(samples.values)), 0.0, 0.0],
        dtype=float,
    )
    fit = least_squares(
        lambda coefficients: (
            (sample_design @ coefficients - samples.values) * weights
        ),
        initial,
        loss=robust_loss,
        f_scale=2.0 * float(noise_sigma),
    )

    yy, xx = np.indices(data.shape, dtype=float)
    full_design = np.column_stack(
        (
            np.ones(data.size, dtype=float),
            (xx.ravel() - x_center) / x_scale,
            (yy.ravel() - y_center) / y_scale,
        )
    )
    plane = (full_design @ fit.x).reshape(data.shape)
    reference_values = np.asarray(reference.values, dtype=float)
    margin = float(reference.diagnostics["foreground_margin"])
    spatial_candidate = plane - margin
    enhancement = np.maximum(spatial_candidate - reference_values, 0.0)

    reference_local = eligible & (
        data <= reference_values + 2.0 * float(noise_sigma)
    )
    reference_strict = eligible & (data <= reference_values)
    eligible_count = int(np.count_nonzero(eligible))
    local_limit = max(
        int(np.count_nonzero(reference_local)),
        int(np.ceil(maximum_local_saturation_fraction * eligible_count)),
    )
    strict_limit = max(
        int(np.count_nonzero(reference_strict)),
        int(np.ceil(maximum_strict_saturation_fraction * eligible_count)),
    )

    selected_blend = 0.0
    selected_local_count = int(np.count_nonzero(reference_local))
    selected_strict_count = int(np.count_nonzero(reference_strict))
    for blend in np.linspace(0.0, 1.0, blend_steps):
        candidate = reference_values + float(blend) * enhancement
        if floor is not None and np.any(candidate[eligible] < float(floor)):
            continue
        local_count = int(
            np.count_nonzero(
                eligible
                & (data <= candidate + 2.0 * float(noise_sigma))
            )
        )
        strict_count = int(
            np.count_nonzero(eligible & (data <= candidate))
        )
        if local_count <= local_limit and strict_count <= strict_limit:
            selected_blend = float(blend)
            selected_local_count = local_count
            selected_strict_count = strict_count

    foreground = reference_values + selected_blend * enhancement
    if floor is not None:
        foreground = np.maximum(foreground, float(floor))

    residuals = samples.values - sample_design @ fit.x
    degrees_of_freedom = max(1, len(samples) - sample_design.shape[1])
    residual_variance = float(
        np.sum(samples.multiplicity * residuals**2) / degrees_of_freedom
    )
    normal_matrix = (
        sample_design.T
        @ (samples.multiplicity[:, None] * sample_design)
    )
    coefficient_covariance = residual_variance * np.linalg.pinv(normal_matrix)
    plane_variance = np.einsum(
        "ij,jk,ik->i",
        full_design,
        coefficient_covariance,
        full_design,
    ).reshape(data.shape)
    active_enhancement = spatial_candidate > reference_values
    variance = (
        selected_blend**2
        * np.maximum(plane_variance, 0.0)
        * active_enhancement
    )

    raw_candidate = reference_values + enhancement
    raw_local_count = int(
        np.count_nonzero(
            eligible
            & (data <= raw_candidate + 2.0 * float(noise_sigma))
        )
    )
    raw_strict_count = int(
        np.count_nonzero(eligible & (data <= raw_candidate))
    )
    return InterpolationResult(
        values=foreground,
        method="conservative",
        variance=variance,
        diagnostics={
            "reference_method": "bt12",
            "reference_foreground": float(reference_values[eligible][0]),
            "reference_saturated_pixel_count": reference.diagnostics[
                "saturated_pixel_count"
            ],
            "reference_independent_saturated_pixel_count": (
                reference.diagnostics["independent_saturated_pixel_count"]
            ),
            "trend_model": "robust_plane",
            "robust_loss": robust_loss,
            "plane_coefficients": fit.x.tolist(),
            "anchor_policy": "hard_bt12_floor",
            "bt12_floor_enforced": True,
            "spatial_candidate_min": float(
                np.min(spatial_candidate[eligible])
            ),
            "spatial_candidate_max": float(
                np.max(spatial_candidate[eligible])
            ),
            "enhancement_min": float(np.min(enhancement[eligible])),
            "enhancement_max": float(np.max(enhancement[eligible])),
            "foreground_below_bt12_count": int(
                np.count_nonzero(
                    eligible & (foreground < reference_values)
                )
            ),
            "blend_factor": selected_blend,
            "blend_steps": int(blend_steps),
            "eligible_pixel_count": eligible_count,
            "local_saturation_limit_count": local_limit,
            "strict_saturation_limit_count": strict_limit,
            "local_saturation_count": selected_local_count,
            "strict_saturation_count": selected_strict_count,
            "unregularized_local_saturation_count": raw_local_count,
            "unregularized_strict_saturation_count": raw_strict_count,
            "sample_residual_rms": float(np.sqrt(residual_variance)),
            "foreground_margin": reference.diagnostics[
                "foreground_margin"
            ],
        },
    )


def fit_liberal_foreground(
    samples: ForegroundSamples,
    image: np.ndarray,
    *,
    wcs: WCS | None = None,
    region_mask: np.ndarray | None = None,
    noise_sigma: float = 0.6,
    min_separation_arcsec: float | None = 8.0,
    foreground_margin: float | None = None,
    robust_loss: str = "soft_l1",
    trend_degree: int = 2,
    target_local_saturation_fraction: float = 0.01,
    maximum_strict_saturation_fraction: float = 0.001,
    bt12_anchor_weight: float = 0.0,
    floor: float | None = 0.0,
    clip_to_sample_range: bool = True,
    fallback_pixel_scale_arcsec: float = 1.2,
    _ordered_floor: InterpolationResult | None = None,
) -> InterpolationResult:
    """Fit a controlled, sample-driven spatial foreground.

    This is the deliberately more permissive counterpart to
    :func:`fit_conservative_foreground`. It fits a robust first- or
    second-degree spatial trend directly to the GTL local-minimum samples and
    uses BT12 only through a one-sided ordering floor: the final liberal
    foreground cannot fall below moderate GTL, which cannot fall below
    conservative GTL or BT12. The absolute level is shifted upward until
    either the requested near-saturation budget or the strict-saturation
    budget is reached. Consequently the model can identify more saturated
    pixels than BT12 without allowing an unconstrained interpolator to imprint
    every local minimum into the extinction map or lower the inferred mass.

    ``target_local_saturation_fraction`` controls pixels satisfying
    ``I_obs <= I_fg + 2 * noise_sigma``. Truly censored pixels satisfy
    ``I_obs <= I_fg`` and are separately capped by
    ``maximum_strict_saturation_fraction``. Such pixels must be computed with
    ``saturation_policy='lower_limit'`` and a positive intensity floor; the
    returned diagnostics recommend ``2 * noise_sigma``. Use
    :meth:`gtlmapping.GTLMapper.compute_liberal` to apply the foreground/
    background feasibility projection and lower-limit calculation together.

    ``bt12_anchor_weight`` controls the two-sided pull on the raw quadratic
    trend. It is separate from the one-sided ordering floor, which is always
    enforced for the named moderate and liberal profiles.
    """

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be a two-dimensional array.")
    if noise_sigma <= 0 or not np.isfinite(noise_sigma):
        raise ValueError("noise_sigma must be positive and finite.")
    if robust_loss not in {"linear", "soft_l1", "huber", "cauchy", "arctan"}:
        raise ValueError(
            "robust_loss must be linear, soft_l1, huber, cauchy, or arctan."
        )
    if trend_degree not in {1, 2}:
        raise ValueError("trend_degree must be one or two.")
    for name, fraction in (
        ("target_local_saturation_fraction", target_local_saturation_fraction),
        (
            "maximum_strict_saturation_fraction",
            maximum_strict_saturation_fraction,
        ),
        ("bt12_anchor_weight", bt12_anchor_weight),
    ):
        if not 0 <= fraction <= 1:
            raise ValueError(f"{name} must lie between zero and one.")
    _validate_samples(samples, "kriging")

    mask = (
        np.ones(data.shape, dtype=bool)
        if region_mask is None
        else np.asarray(region_mask, dtype=bool)
    )
    if mask.shape != data.shape:
        raise ValueError("region_mask must match image.shape.")
    eligible = mask & np.isfinite(data)
    if not np.any(eligible):
        raise ValueError("No finite image pixels are inside region_mask.")

    if _ordered_floor is None:
        _ordered_floor = fit_moderate_foreground(
            samples,
            data,
            wcs=wcs,
            region_mask=mask,
            noise_sigma=noise_sigma,
            min_separation_arcsec=min_separation_arcsec,
            foreground_margin=foreground_margin,
            robust_loss=robust_loss,
            trend_degree=trend_degree,
            floor=floor,
            clip_to_sample_range=clip_to_sample_range,
            fallback_pixel_scale_arcsec=fallback_pixel_scale_arcsec,
        )

    margin = (
        2.0 * float(noise_sigma)
        if foreground_margin is None
        else float(foreground_margin)
    )
    if margin < 0 or not np.isfinite(margin):
        raise ValueError("foreground_margin must be non-negative and finite.")

    ny, nx = data.shape
    x_center = 0.5 * (nx - 1)
    y_center = 0.5 * (ny - 1)
    x_scale = max(x_center, 1.0)
    y_scale = max(y_center, 1.0)

    def design_matrix(cols: np.ndarray, rows: np.ndarray) -> np.ndarray:
        x = (np.asarray(cols, dtype=float) - x_center) / x_scale
        y = (np.asarray(rows, dtype=float) - y_center) / y_scale
        columns = [np.ones_like(x), x, y]
        if trend_degree == 2:
            columns.extend((x**2, x * y, y**2))
        return np.column_stack(columns)

    sample_design = design_matrix(samples.cols, samples.rows)
    weights = np.sqrt(samples.multiplicity.astype(float))
    initial = np.zeros(sample_design.shape[1], dtype=float)
    initial[0] = float(np.median(samples.values))
    fit = least_squares(
        lambda coefficients: (
            (sample_design @ coefficients - samples.values) * weights
        ),
        initial,
        loss=robust_loss,
        f_scale=2.0 * float(noise_sigma),
    )

    yy, xx = np.indices(data.shape, dtype=float)
    full_design = design_matrix(xx.ravel(), yy.ravel())
    trend = (full_design @ fit.x).reshape(data.shape)
    if clip_to_sample_range:
        trend = np.clip(
            trend,
            float(np.min(samples.values)),
            float(np.max(samples.values)),
        )
    spatial_candidate = trend - margin

    selected_separation = (
        samples.min_separation_arcsec
        if min_separation_arcsec is None
        else float(min_separation_arcsec)
    )
    reference: InterpolationResult | None
    try:
        reference = estimate_bt12_foreground(
            data,
            wcs=wcs,
            region_mask=mask,
            noise_sigma=noise_sigma,
            min_separation_arcsec=selected_separation,
            foreground_margin=margin,
            fallback_pixel_scale_arcsec=fallback_pixel_scale_arcsec,
        )
    except InsufficientSamplesError:
        reference = None

    if bt12_anchor_weight > 0:
        if reference is None:
            raise InsufficientSamplesError(
                "bt12_anchor_weight is positive, but BT12 has no independent "
                "saturated reference pixels."
            )
        spatial_candidate = (
            (1.0 - float(bt12_anchor_weight)) * spatial_candidate
            + float(bt12_anchor_weight) * reference.values
        )

    eligible_count = int(np.count_nonzero(eligible))
    if reference is None:
        reference_foreground = None
        reference_local_count = 0
        reference_strict_count = 0
        reference_saturated_count = 0
        reference_independent_count = 0
    else:
        reference_values = np.asarray(reference.values, dtype=float)
        reference_foreground = float(reference_values[eligible][0])
        reference_local_count = int(
            np.count_nonzero(
                eligible
                & (data <= reference_values + 2.0 * float(noise_sigma))
            )
        )
        reference_strict_count = int(
            np.count_nonzero(eligible & (data <= reference_values))
        )
        reference_saturated_count = int(
            reference.diagnostics["saturated_pixel_count"]
        )
        reference_independent_count = int(
            reference.diagnostics["independent_saturated_pixel_count"]
        )

    local_limit = min(
        eligible_count,
        max(
            reference_local_count + 1,
            int(
                np.ceil(
                    float(target_local_saturation_fraction) * eligible_count
                )
            ),
        ),
    )
    strict_limit = min(
        eligible_count,
        int(
            np.ceil(
                float(maximum_strict_saturation_fraction) * eligible_count
            )
        ),
    )

    def shifted_foreground(shift: float) -> np.ndarray:
        values = spatial_candidate + float(shift)
        if floor is not None:
            values = np.maximum(values, float(floor))
        return values

    def saturation_counts(shift: float) -> tuple[int, int]:
        values = shifted_foreground(shift)
        local = int(
            np.count_nonzero(
                eligible & (data <= values + 2.0 * float(noise_sigma))
            )
        )
        strict = int(np.count_nonzero(eligible & (data <= values)))
        return local, strict

    residual = data[eligible] - spatial_candidate[eligible]
    residual_scale = max(float(np.ptp(residual)), 2.0 * float(noise_sigma), 1.0)
    lower_shift = float(np.min(residual) - 2.0 * noise_sigma - residual_scale)
    upper_shift = float(np.max(residual) + residual_scale)
    lower_counts = saturation_counts(lower_shift)
    if lower_counts[0] > local_limit or lower_counts[1] > strict_limit:
        raise ValueError(
            "The requested saturation budgets are incompatible with the "
            "foreground floor. Lower floor or relax the budgets."
        )
    for _ in range(80):
        midpoint = 0.5 * (lower_shift + upper_shift)
        local_count, strict_count = saturation_counts(midpoint)
        if local_count <= local_limit and strict_count <= strict_limit:
            lower_shift = midpoint
        else:
            upper_shift = midpoint

    selected_shift = float(lower_shift)
    foreground = shifted_foreground(selected_shift)
    selected_local_count, selected_strict_count = saturation_counts(
        selected_shift
    )

    residuals = samples.values - sample_design @ fit.x
    degrees_of_freedom = max(1, len(samples) - sample_design.shape[1])
    residual_variance = float(
        np.sum(samples.multiplicity * residuals**2) / degrees_of_freedom
    )
    normal_matrix = (
        sample_design.T
        @ (samples.multiplicity[:, None] * sample_design)
    )
    coefficient_covariance = residual_variance * np.linalg.pinv(normal_matrix)
    trend_variance = np.einsum(
        "ij,jk,ik->i",
        full_design,
        coefficient_covariance,
        full_design,
    ).reshape(data.shape)
    if bt12_anchor_weight > 0:
        trend_variance *= (1.0 - float(bt12_anchor_weight)) ** 2
    if floor is not None:
        trend_variance = np.where(
            spatial_candidate + selected_shift > float(floor),
            trend_variance,
            0.0,
        )

    raw_foreground = foreground
    raw_local_count = selected_local_count
    raw_strict_count = selected_strict_count
    ordered_floor = np.asarray(_ordered_floor.values, dtype=float)
    if ordered_floor.shape != data.shape:
        raise ValueError("The ordered foreground floor must match image.shape.")
    if np.any(eligible & ~np.isfinite(ordered_floor)):
        raise ValueError(
            "The ordered foreground floor is non-finite inside region_mask."
        )
    floor_local_count = int(
        np.count_nonzero(
            eligible
            & (data <= ordered_floor + 2.0 * float(noise_sigma))
        )
    )
    floor_strict_count = int(
        np.count_nonzero(eligible & (data <= ordered_floor))
    )
    local_limit = max(local_limit, floor_local_count)
    strict_limit = max(strict_limit, floor_strict_count)
    ordered_enhancement = np.maximum(raw_foreground - ordered_floor, 0.0)

    def ordered_counts(blend: float) -> tuple[int, int]:
        candidate = ordered_floor + float(blend) * ordered_enhancement
        local_count = int(
            np.count_nonzero(
                eligible
                & (data <= candidate + 2.0 * float(noise_sigma))
            )
        )
        strict_count = int(np.count_nonzero(eligible & (data <= candidate)))
        return local_count, strict_count

    ordering_blend = 1.0
    selected_local_count, selected_strict_count = ordered_counts(
        ordering_blend
    )
    if (
        selected_local_count > local_limit
        or selected_strict_count > strict_limit
    ):
        lower_blend = 0.0
        upper_blend = 1.0
        for _ in range(80):
            midpoint = 0.5 * (lower_blend + upper_blend)
            local_count, strict_count = ordered_counts(midpoint)
            if local_count <= local_limit and strict_count <= strict_limit:
                lower_blend = midpoint
            else:
                upper_blend = midpoint
        ordering_blend = float(lower_blend)
        selected_local_count, selected_strict_count = ordered_counts(
            ordering_blend
        )
    foreground = ordered_floor + ordering_blend * ordered_enhancement

    floor_variance = (
        np.zeros(data.shape, dtype=float)
        if _ordered_floor.variance is None
        else np.maximum(
            np.asarray(_ordered_floor.variance, dtype=float),
            0.0,
        )
    )
    active_enhancement = raw_foreground > ordered_floor
    trend_variance = floor_variance + (
        ordering_blend**2 * trend_variance * active_enhancement
    )

    return InterpolationResult(
        values=foreground,
        method="liberal",
        variance=np.maximum(trend_variance, 0.0),
        diagnostics={
            "reference_method": "bt12_diagnostic_only",
            "reference_available": reference is not None,
            "reference_foreground": reference_foreground,
            "reference_saturated_pixel_count": reference_saturated_count,
            "reference_independent_saturated_pixel_count": (
                reference_independent_count
            ),
            "reference_local_saturation_count": reference_local_count,
            "reference_strict_saturation_count": reference_strict_count,
            "trend_model": (
                "robust_quadratic" if trend_degree == 2 else "robust_plane"
            ),
            "trend_degree": int(trend_degree),
            "robust_loss": robust_loss,
            "trend_coefficients": fit.x.tolist(),
            "anchor_policy": "ordered_one_sided_floor",
            "bt12_anchor_weight": float(bt12_anchor_weight),
            "bt12_floor_enforced": True,
            "ordering_floor_method": _ordered_floor.method,
            "ordering_floor_enforced": True,
            "ordering_blend_factor": ordering_blend,
            "raw_foreground_below_ordering_floor_count": int(
                np.count_nonzero(
                    eligible & (raw_foreground < ordered_floor)
                )
            ),
            "foreground_below_ordering_floor_count": int(
                np.count_nonzero(eligible & (foreground < ordered_floor))
            ),
            "ordering_floor_local_saturation_count": floor_local_count,
            "ordering_floor_strict_saturation_count": floor_strict_count,
            "raw_profile_local_saturation_count": raw_local_count,
            "raw_profile_strict_saturation_count": raw_strict_count,
            "variance_combination": (
                "ordered_floor_plus_scaled_candidate_without_covariance"
            ),
            "foreground_margin": margin,
            "level_shift": selected_shift,
            "clip_to_sample_range": bool(clip_to_sample_range),
            "spatial_candidate_min": float(
                np.min(spatial_candidate[eligible])
            ),
            "spatial_candidate_max": float(
                np.max(spatial_candidate[eligible])
            ),
            "eligible_pixel_count": eligible_count,
            "target_local_saturation_fraction": float(
                target_local_saturation_fraction
            ),
            "local_saturation_limit_count": local_limit,
            "local_saturation_count": selected_local_count,
            "maximum_strict_saturation_fraction": float(
                maximum_strict_saturation_fraction
            ),
            "strict_saturation_limit_count": strict_limit,
            "strict_saturation_count": selected_strict_count,
            "more_local_saturation_than_bt12": (
                selected_local_count > reference_local_count
            ),
            "requires_lower_limit_policy": selected_strict_count > 0,
            "recommended_intensity_floor": 2.0 * float(noise_sigma),
            "sample_residual_rms": float(np.sqrt(residual_variance)),
            "foreground_nonfinite_count": int(
                np.count_nonzero(~np.isfinite(foreground))
            ),
        },
    )


def fit_moderate_foreground(
    samples: ForegroundSamples,
    image: np.ndarray,
    *,
    wcs: WCS | None = None,
    region_mask: np.ndarray | None = None,
    noise_sigma: float = 0.6,
    min_separation_arcsec: float | None = 8.0,
    foreground_margin: float | None = None,
    robust_loss: str = "soft_l1",
    trend_degree: int = 2,
    target_local_saturation_fraction: float = 0.005,
    maximum_strict_saturation_fraction: float = 0.0001,
    bt12_anchor_weight: float = 0.5,
    floor: float | None = 0.0,
    clip_to_sample_range: bool = True,
    fallback_pixel_scale_arcsec: float = 1.2,
) -> InterpolationResult:
    """Fit the intermediate GTL foreground preset.

    Moderate GTL uses the same controlled robust spatial trend as liberal GTL
    but defaults to a 50-percent soft BT12 pull, a 0.5-percent near-saturation
    budget, and a 0.01-percent strict-censoring ceiling. Conservative GTL is a
    one-sided pointwise floor, so moderate GTL cannot lower its foreground or
    inferred surface density on jointly valid, uncensored pixels.

    The parameters remain explicit and can be changed for sensitivity tests.
    Any strictly saturated pixels are still censored measurements and should
    be computed with :meth:`gtlmapping.GTLMapper.compute_moderate`.
    """

    conservative_floor = fit_conservative_foreground(
        samples,
        image,
        wcs=wcs,
        region_mask=region_mask,
        noise_sigma=noise_sigma,
        min_separation_arcsec=min_separation_arcsec,
        foreground_margin=foreground_margin,
        robust_loss=robust_loss,
        floor=floor,
        fallback_pixel_scale_arcsec=fallback_pixel_scale_arcsec,
    )
    result = fit_liberal_foreground(
        samples,
        image,
        wcs=wcs,
        region_mask=region_mask,
        noise_sigma=noise_sigma,
        min_separation_arcsec=min_separation_arcsec,
        foreground_margin=foreground_margin,
        robust_loss=robust_loss,
        trend_degree=trend_degree,
        target_local_saturation_fraction=target_local_saturation_fraction,
        maximum_strict_saturation_fraction=(
            maximum_strict_saturation_fraction
        ),
        bt12_anchor_weight=bt12_anchor_weight,
        floor=floor,
        clip_to_sample_range=clip_to_sample_range,
        fallback_pixel_scale_arcsec=fallback_pixel_scale_arcsec,
        _ordered_floor=conservative_floor,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "profile": "moderate",
            "profile_description": (
                "conservative_floor_with_soft_bt12_spatial_trend"
            ),
        }
    )
    return InterpolationResult(
        values=result.values,
        method="moderate",
        variance=result.variance,
        constraint_mask=result.constraint_mask,
        diagnostics=diagnostics,
    )


def _validate_samples(samples: ForegroundSamples, method: str) -> None:
    minimum = {
        "flat": 1,
        "gaussian": 1,
        "cauchy": 1,
        "rbf": 3,
        "spline": 9,
        "kriging": 3,
    }[method]
    if len(samples) < minimum:
        raise InsufficientSamplesError(
            f"{method} interpolation needs at least {minimum} unique samples; "
            f"received {len(samples)}."
        )
    if method in {"rbf", "spline", "kriging"}:
        centered = samples.points_xy - np.mean(samples.points_xy, axis=0)
        if np.linalg.matrix_rank(centered) < 2:
            raise InsufficientSamplesError(
                f"{method} interpolation needs non-collinear sample coordinates."
            )


def _weighted_interpolation(
    samples: ForegroundSamples,
    shape: tuple[int, int],
    *,
    kernel: str,
    length_scale_pixels: float,
    chunk_size: int,
    use_multiplicity: bool,
) -> np.ndarray:
    if length_scale_pixels <= 0:
        raise ValueError("length_scale_pixels must be positive.")
    points = samples.points_xy.astype(float)
    values = samples.values.astype(float)
    yy, xx = np.indices(shape, dtype=float)
    targets = np.column_stack((xx.ravel(), yy.ravel()))
    output = np.empty(len(targets), dtype=float)
    scale2 = length_scale_pixels**2
    for start in range(0, len(targets), chunk_size):
        stop = min(start + chunk_size, len(targets))
        delta = targets[start:stop, None, :] - points[None, :, :]
        radius2 = np.sum(delta**2, axis=2)
        if kernel == "gaussian":
            weights = np.exp(-0.5 * radius2 / scale2)
        else:
            weights = 1.0 / (1.0 + radius2 / scale2)
        if use_multiplicity:
            weights *= samples.multiplicity[None, :]
        denominator = np.sum(weights, axis=1)
        output[start:stop] = np.divide(
            weights @ values,
            denominator,
            out=np.full(stop - start, np.nan),
            where=denominator > 0,
        )
    return output.reshape(shape)


def interpolate_foreground(
    samples: ForegroundSamples,
    shape: tuple[int, int],
    *,
    method: str = "kriging",
    foreground_margin: float = 1.2,
    floor: float | None = 0.0,
    clip_to_sample_range: bool = True,
    variogram_model: str = "gaussian",
    length_scale_pixels: float = 500.0,
    smoothing: float = 0.0,
    chunk_size: int = 50_000,
    use_multiplicity: bool = True,
    kriging_duplicate_policy: str = "aggregate",
    fill_interpolation_gaps: bool = True,
    **kwargs: Any,
) -> InterpolationResult:
    """Interpolate samples over an image and subtract a conservative margin.

    The stable kriging path aggregates overlapping-window duplicate
    coordinates and uses a pseudo-inverse. Set
    ``kriging_duplicate_policy='repeat'`` to reproduce the prototype
    notebook's repeated-coordinate fit. Any numerical holes in the
    interpolated intensity surface are filled from the nearest accepted
    foreground sample by default and counted in the diagnostics.
    """

    normalized = method.lower()
    if normalized not in {"kriging", "rbf", "spline", "gaussian", "cauchy", "flat"}:
        raise ValueError(
            "method must be kriging, rbf, spline, gaussian, cauchy, or flat."
        )
    if foreground_margin < 0:
        raise ValueError("foreground_margin must be non-negative.")
    if kriging_duplicate_policy not in {"aggregate", "repeat"}:
        raise ValueError(
            "kriging_duplicate_policy must be 'aggregate' or 'repeat'."
        )
    _validate_samples(samples, normalized)

    variance = None
    diagnostics: dict[str, Any] = {
        "sample_count": len(samples),
        "raw_detection_count": samples.raw_detection_count,
        "foreground_margin": float(foreground_margin),
        "clip_to_sample_range": bool(clip_to_sample_range),
        "use_multiplicity": bool(use_multiplicity),
        "fill_interpolation_gaps": bool(fill_interpolation_gaps),
    }
    if normalized == "flat":
        values = (
            np.repeat(samples.values, samples.multiplicity)
            if use_multiplicity
            else samples.values
        )
        prediction = np.full(shape, np.median(values), dtype=float)
    elif normalized in {"gaussian", "cauchy"}:
        prediction = _weighted_interpolation(
            samples,
            shape,
            kernel=normalized,
            length_scale_pixels=length_scale_pixels,
            chunk_size=chunk_size,
            use_multiplicity=use_multiplicity,
        )
        diagnostics["length_scale_pixels"] = float(length_scale_pixels)
    elif normalized == "rbf":
        interpolator = RBFInterpolator(
            samples.points_xy.astype(float),
            samples.values.astype(float),
            kernel=kwargs.pop("kernel", "thin_plate_spline"),
            smoothing=smoothing,
            **kwargs,
        )
        yy, xx = np.indices(shape, dtype=float)
        targets = np.column_stack((xx.ravel(), yy.ravel()))
        flat = np.empty(len(targets), dtype=float)
        for start in range(0, len(targets), chunk_size):
            stop = min(start + chunk_size, len(targets))
            flat[start:stop] = interpolator(targets[start:stop])
        prediction = flat.reshape(shape)
        diagnostics["smoothing"] = float(smoothing)
    elif normalized == "spline":
        spline_smoothing = float(
            kwargs.pop("s", 200.0 if smoothing == 0.0 else smoothing)
        )
        spline = SmoothBivariateSpline(
            samples.cols,
            samples.rows,
            samples.values,
            kx=int(kwargs.pop("kx", 2)),
            ky=int(kwargs.pop("ky", 2)),
            s=spline_smoothing,
            **kwargs,
        )
        yy, xx = np.indices(shape, dtype=float)
        prediction = spline.ev(xx.ravel(), yy.ravel()).reshape(shape)
        diagnostics["smoothing"] = spline_smoothing
    else:
        try:
            from pykrige.ok import OrdinaryKriging
        except ImportError as exc:  # pragma: no cover - dependency is declared
            raise ImportError(
                "Kriging requires PyKrige; install GTLMapping with its dependencies."
            ) from exc
        if use_multiplicity and kriging_duplicate_policy == "repeat":
            repeated = np.repeat(np.arange(len(samples)), samples.multiplicity)
        else:
            repeated = np.arange(len(samples))
        kwargs.setdefault("pseudo_inv", True)
        kwargs.setdefault("pseudo_inv_type", "pinv")
        kriging = OrdinaryKriging(
            samples.cols[repeated].astype(float),
            samples.rows[repeated].astype(float),
            samples.values[repeated].astype(float),
            variogram_model=variogram_model,
            verbose=False,
            enable_plotting=False,
            **kwargs,
        )
        grid_x = np.arange(shape[1], dtype=float)
        grid_y = np.arange(shape[0], dtype=float)
        predicted, kriging_variance = kriging.execute("grid", grid_x, grid_y)
        prediction = np.ma.filled(predicted, np.nan)
        variance = np.maximum(np.ma.filled(kriging_variance, np.nan), 0.0)
        diagnostics["variogram_model"] = variogram_model
        diagnostics["kriging_duplicate_policy"] = kriging_duplicate_policy
        diagnostics["kriging_pseudo_inverse"] = bool(kwargs["pseudo_inv"])

    prediction = np.asarray(prediction, dtype=float)
    interpolation_gaps = ~np.isfinite(prediction)
    gap_count = int(np.count_nonzero(interpolation_gaps))
    diagnostics["interpolation_gap_count"] = gap_count
    if gap_count and fill_interpolation_gaps:
        gap_rows, gap_cols = np.nonzero(interpolation_gaps)
        nearest = griddata(
            samples.points_xy.astype(float),
            samples.values.astype(float),
            (gap_cols.astype(float), gap_rows.astype(float)),
            method="nearest",
        )
        prediction[gap_rows, gap_cols] = nearest
        diagnostics["filled_interpolation_gap_count"] = int(
            np.count_nonzero(np.isfinite(nearest))
        )
    remaining_gaps = int(np.count_nonzero(~np.isfinite(prediction)))
    diagnostics["remaining_interpolation_gap_count"] = remaining_gaps
    if remaining_gaps:
        raise RuntimeError(
            f"{normalized} interpolation left {remaining_gaps} non-finite pixels."
        )
    if variance is not None:
        variance_gaps = ~np.isfinite(variance)
        variance_gap_count = int(np.count_nonzero(variance_gaps))
        diagnostics["variance_gap_count"] = variance_gap_count
        if variance_gap_count:
            finite_variance = variance[np.isfinite(variance)]
            conservative_fill = (
                float(np.max(finite_variance)) if finite_variance.size else 0.0
            )
            variance[variance_gaps] = conservative_fill
            diagnostics["variance_gap_fill"] = conservative_fill
    diagnostics["raw_prediction_min"] = float(np.nanmin(prediction))
    diagnostics["raw_prediction_max"] = float(np.nanmax(prediction))
    if clip_to_sample_range:
        lower = float(np.nanmin(samples.values))
        upper = float(np.nanmax(samples.values))
        clipped = np.clip(prediction, lower, upper)
        diagnostics["clipped_pixel_count"] = int(
            np.count_nonzero(~np.isclose(clipped, prediction, equal_nan=True))
        )
        prediction = clipped
    foreground = prediction - foreground_margin
    if floor is not None:
        foreground = np.maximum(foreground, float(floor))
    return InterpolationResult(
        values=foreground,
        method=normalized,
        variance=variance,
        diagnostics=diagnostics,
    )


def cross_validate_foreground(
    samples: ForegroundSamples,
    *,
    method: str = "gaussian",
    foreground_margin: float = 0.0,
    **kwargs: Any,
) -> dict[str, float]:
    """Leave-one-out diagnostics at the sample locations."""

    if len(samples) < 4:
        raise InsufficientSamplesError(
            "Cross-validation requires at least four unique samples."
        )
    predictions = np.empty(len(samples), dtype=float)
    for index in range(len(samples)):
        keep = np.arange(len(samples)) != index
        subset = ForegroundSamples(
            rows=samples.rows[keep],
            cols=samples.cols[keep],
            values=samples.values[keep],
            multiplicity=samples.multiplicity[keep],
            accepted_windows=int(np.sum(samples.multiplicity[keep])),
            rejected_windows=samples.rejected_windows,
            total_windows=samples.total_windows,
            pixel_scale_arcsec=samples.pixel_scale_arcsec,
            min_separation_arcsec=samples.min_separation_arcsec,
        )
        result = interpolate_foreground(
            subset,
            shape=(int(samples.rows[index]) + 1, int(samples.cols[index]) + 1),
            method=method,
            foreground_margin=foreground_margin,
            **kwargs,
        )
        predictions[index] = result.values[
            int(samples.rows[index]), int(samples.cols[index])
        ]
    residuals = predictions - (samples.values - foreground_margin)
    return {
        "mae": float(np.mean(np.abs(residuals))),
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "bias": float(np.mean(residuals)),
        "max_abs_error": float(np.max(np.abs(residuals))),
    }
