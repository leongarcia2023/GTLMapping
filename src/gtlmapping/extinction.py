"""Radiative-transfer calculation for MIREX maps."""

from __future__ import annotations

import numpy as np

from .models import UncertaintyResult

PARSEC_CM = 3.0856775814913673e18
SOLAR_MASS_G = 1.988409870698051e33
G_PER_CM2_TO_MSUN_PER_PC2 = PARSEC_CM**2 / SOLAR_MASS_G


def unresolved_transmission(observed, foreground, threshold=None) -> np.ndarray:
    """Flag nonpositive or weak transmission; this is not an opacity test.

    A threshold uses the intensity unit of the inputs. Without one, only
    nonpositive transmission is identified. Missing threshold pixels remain
    unresolved rather than becoming detections.
    """
    obs, fg = np.broadcast_arrays(np.asarray(observed, float), np.asarray(foreground, float))
    finite = np.isfinite(obs) & np.isfinite(fg)
    if threshold is None:
        return finite & (obs <= fg)
    limit = np.broadcast_to(np.asarray(threshold, float), obs.shape)
    if np.any(np.isfinite(limit) & (limit <= 0)) or np.any(np.isinf(limit)):
        raise ValueError("The detection threshold must be positive; NaN denotes missing sensitivity.")
    return finite & ((obs - fg <= limit) | ~np.isfinite(limit))


def transmission_std(observed_std, foreground_std=0.0, covariance=0.0):
    """Standard deviation of I_obs-I_fg with supplied same-pixel covariance."""
    obs, fg, cov = np.broadcast_arrays(np.asarray(observed_std, float),
                                      np.asarray(foreground_std, float), np.asarray(covariance, float))
    if np.any(obs < 0) or np.any(fg < 0):
        raise ValueError("Standard deviations must be nonnegative.")
    if np.any(np.abs(cov) > obs * fg + 1e-12):
        raise ValueError("Covariance exceeds the product of the standard deviations.")
    return np.sqrt(np.maximum(obs**2 + fg**2 - 2*cov, 0.0))


def convert_surface_density(
    surface_density: np.ndarray | np.ma.MaskedArray,
    *,
    to: str = "Msun/pc2",
) -> np.ndarray | np.ma.MaskedArray:
    """Convert surface density from g cm⁻² to solar masses pc⁻²."""

    normalized = to.lower().replace(" ", "")
    if normalized not in {"msun/pc2", "m_sun/pc^2", "solarmass/pc2"}:
        raise ValueError("Only conversion to Msun/pc2 is currently supported.")
    return surface_density * G_PER_CM2_TO_MSUN_PER_PC2


def compute_extinction(
    observed: np.ndarray,
    background: np.ndarray,
    foreground: np.ndarray,
    *,
    kappa_cm2_g: float = 7.5,
    bright_pixel_policy: str = "allow",
    saturation_policy: str = "mask",
    intensity_floor: float | np.ndarray | None = None,
    detection_threshold: float | np.ndarray | None = None,
) -> tuple[
    np.ma.MaskedArray,
    np.ma.MaskedArray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Compute optical depth and mass surface density.

    The adopted radiative-transfer equation is

    ``tau = -ln((I_obs - I_fg) / (I_bg - I_fg))``.

    Pixels with ``I_bg <= I_fg`` are always invalid. Pixels with
    ``I_obs <= I_fg`` are masked by default. With
    ``saturation_policy='lower_limit'``, an explicitly supplied positive
    ``intensity_floor`` replaces their transmitted intensity and their
    values become conditional limits. The floor also sets the detection
    threshold unless detection_threshold is supplied. The returned saturated
    mask retains the strict zero-crossing test; unresolved_transmission gives
    the full limit mask, including weak positive transmission.
    ``bright_pixel_policy='allow'`` retains negative optical depths,
    matching BT09's bias-avoidance treatment.
    """

    if kappa_cm2_g <= 0 or not np.isfinite(kappa_cm2_g):
        raise ValueError("kappa_cm2_g must be positive and finite.")
    if bright_pixel_policy not in {"allow", "zero", "mask"}:
        raise ValueError("bright_pixel_policy must be 'allow', 'zero', or 'mask'.")
    if saturation_policy not in {"mask", "lower_limit"}:
        raise ValueError("saturation_policy must be 'mask' or 'lower_limit'.")

    obs, bg, fg = np.broadcast_arrays(
        np.asarray(observed, dtype=float),
        np.asarray(background, dtype=float),
        np.asarray(foreground, dtype=float),
    )
    numerator = obs - fg
    denominator = bg - fg
    finite = np.isfinite(obs) & np.isfinite(bg) & np.isfinite(fg)
    saturated = finite & (numerator <= 0)
    selected_threshold = intensity_floor if detection_threshold is None else detection_threshold
    unresolved = unresolved_transmission(obs, fg, selected_threshold)
    invalid_background = finite & (denominator <= 0)
    valid = finite & ~unresolved & ~invalid_background
    effective_numerator = np.array(numerator, dtype=float, copy=True)

    if saturation_policy == "lower_limit":
        if intensity_floor is None:
            raise ValueError(
                "intensity_floor is required when saturation_policy='lower_limit'."
            )
        floor = np.broadcast_to(np.asarray(intensity_floor, dtype=float), obs.shape)
        if np.any(np.isfinite(floor) & (floor <= 0)) or np.any(np.isinf(floor)):
            raise ValueError("intensity_floor must be positive; NaN denotes missing sensitivity.")
        lower_limit = unresolved & finite & ~invalid_background & (floor <= denominator)
        if selected_threshold is not None:
            lower_limit &= np.isfinite(np.broadcast_to(selected_threshold, obs.shape))
        effective_numerator[lower_limit] = floor[lower_limit]
        valid |= lower_limit

    ratio = np.full(obs.shape, np.nan, dtype=float)
    np.divide(effective_numerator, denominator, out=ratio, where=valid)
    bright = valid & (ratio > 1)
    if bright_pixel_policy == "mask":
        valid &= ~bright

    tau_values = np.full(obs.shape, np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_values[valid] = -np.log(ratio[valid])
    if bright_pixel_policy == "zero":
        tau_values[bright] = 0.0

    mask = ~valid | ~np.isfinite(tau_values)
    tau = np.ma.array(tau_values, mask=mask)
    sigma = tau / float(kappa_cm2_g)
    return tau, sigma, saturated, invalid_background, bright


def propagate_uncertainty(
    observed: np.ndarray,
    background: np.ndarray,
    foreground: np.ndarray,
    *,
    observed_std: float | np.ndarray = 0.0,
    background_std: float | np.ndarray = 0.0,
    foreground_std: float | np.ndarray = 0.0,
    kappa_cm2_g: float = 7.5,
    kappa_std_cm2_g: float = 0.0,
    additional_mask: np.ndarray | None = None,
) -> UncertaintyResult:
    """Propagate independent first-order intensity and opacity uncertainties.

    For ``tau = -ln(N/D)``, where ``N = I_obs - I_fg`` and
    ``D = I_bg - I_fg``, the derivatives are ``-1/N``, ``1/D``, and
    ``1/N - 1/D`` for observed, background, and foreground intensity.
    This approximation is intentionally masked at saturated or otherwise
    non-physical pixels; Monte Carlo or censored inference is preferable
    near ``N = 0``.
    """

    if kappa_cm2_g <= 0 or not np.isfinite(kappa_cm2_g):
        raise ValueError("kappa_cm2_g must be positive and finite.")
    if kappa_std_cm2_g < 0 or not np.isfinite(kappa_std_cm2_g):
        raise ValueError("kappa_std_cm2_g must be non-negative and finite.")

    obs, bg, fg, obs_std, bg_std, fg_std = np.broadcast_arrays(
        np.asarray(observed, dtype=float),
        np.asarray(background, dtype=float),
        np.asarray(foreground, dtype=float),
        np.asarray(observed_std, dtype=float),
        np.asarray(background_std, dtype=float),
        np.asarray(foreground_std, dtype=float),
    )
    if np.any(
        (np.isfinite(obs_std) & (obs_std < 0))
        | (np.isfinite(bg_std) & (bg_std < 0))
        | (np.isfinite(fg_std) & (fg_std < 0))
    ):
        raise ValueError("Intensity standard deviations must be non-negative.")

    numerator = obs - fg
    denominator = bg - fg
    valid = (
        np.isfinite(obs)
        & np.isfinite(bg)
        & np.isfinite(fg)
        & np.isfinite(obs_std)
        & np.isfinite(bg_std)
        & np.isfinite(fg_std)
        & (numerator > 0)
        & (denominator > 0)
    )
    if additional_mask is not None:
        extra = np.asarray(additional_mask, dtype=bool)
        if extra.shape != obs.shape:
            raise ValueError("additional_mask must match the broadcast image shape.")
        valid &= ~extra

    observed_variance = np.zeros(obs.shape, dtype=float)
    background_variance = np.zeros(obs.shape, dtype=float)
    foreground_variance = np.zeros(obs.shape, dtype=float)
    tau_values = np.full(obs.shape, np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        observed_variance[valid] = (obs_std[valid] / numerator[valid]) ** 2
        background_variance[valid] = (bg_std[valid] / denominator[valid]) ** 2
        foreground_derivative = (
            1.0 / numerator[valid] - 1.0 / denominator[valid]
        )
        foreground_variance[valid] = (
            foreground_derivative * fg_std[valid]
        ) ** 2
        tau_values[valid] = -np.log(numerator[valid] / denominator[valid])

    tau_variance = (
        observed_variance + background_variance + foreground_variance
    )
    opacity_variance = np.zeros(obs.shape, dtype=float)
    opacity_variance[valid] = (
        tau_values[valid]
        * float(kappa_std_cm2_g)
        / float(kappa_cm2_g) ** 2
    ) ** 2
    sigma_components = {
        "observed": observed_variance / float(kappa_cm2_g) ** 2,
        "background": background_variance / float(kappa_cm2_g) ** 2,
        "foreground": foreground_variance / float(kappa_cm2_g) ** 2,
        "opacity": opacity_variance,
    }
    sigma_variance = (
        tau_variance / float(kappa_cm2_g) ** 2 + opacity_variance
    )
    mask = ~valid
    components = {
        name: np.ma.array(np.sqrt(variance), mask=mask)
        for name, variance in sigma_components.items()
    }
    return UncertaintyResult(
        optical_depth_std=np.ma.array(np.sqrt(tau_variance), mask=mask),
        surface_density_std=np.ma.array(np.sqrt(sigma_variance), mask=mask),
        components=components,
        diagnostics={
            "method": "independent_first_order",
            "kappa_std_cm2_g": float(kappa_std_cm2_g),
            "valid_fraction": float(np.count_nonzero(valid) / valid.size),
        },
    )
