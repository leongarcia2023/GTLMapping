from __future__ import annotations

import numpy as np
import pytest

from gtlmapping.extinction import (
    G_PER_CM2_TO_MSUN_PER_PC2,
    compute_extinction,
    convert_surface_density,
    propagate_uncertainty,
)


def test_radiative_transfer_equation_and_masks() -> None:
    foreground = np.array([[20.0, 20.0, 20.0]])
    background = np.array([[100.0, 100.0, 100.0]])
    observed = np.array([[60.0, 20.0, 120.0]])

    tau, sigma, saturated, invalid_background, bright = compute_extinction(
        observed,
        background,
        foreground,
        kappa_cm2_g=7.5,
        bright_pixel_policy="allow",
    )

    assert tau[0, 0] == pytest.approx(-np.log(0.5))
    assert sigma[0, 0] == pytest.approx(-np.log(0.5) / 7.5)
    assert tau.mask[0, 1]
    assert saturated[0, 1]
    assert not invalid_background[0, 1]
    assert bright[0, 2]
    assert tau[0, 2] < 0


def test_bright_pixels_can_be_zeroed() -> None:
    tau, _, _, _, bright = compute_extinction(
        np.array([120.0]),
        np.array([100.0]),
        np.array([20.0]),
        bright_pixel_policy="zero",
    )
    assert bright[0]
    assert tau[0] == 0


def test_saturated_pixels_can_be_reported_as_lower_limits() -> None:
    tau, _, saturated, invalid_background, _ = compute_extinction(
        np.array([20.0]),
        np.array([100.0]),
        np.array([20.0]),
        saturation_policy="lower_limit",
        intensity_floor=2.0,
    )

    assert saturated[0]
    assert not invalid_background[0]
    assert not tau.mask[0]
    assert tau[0] == pytest.approx(-np.log(2.0 / 80.0))


def test_first_order_uncertainty_matches_closed_form() -> None:
    uncertainty = propagate_uncertainty(
        np.array([60.0]),
        np.array([100.0]),
        np.array([20.0]),
        observed_std=2.0,
        background_std=4.0,
        foreground_std=1.0,
        kappa_cm2_g=10.0,
        kappa_std_cm2_g=2.0,
    )

    numerator = 40.0
    denominator = 80.0
    tau = -np.log(numerator / denominator)
    tau_variance = (
        (2.0 / numerator) ** 2
        + (4.0 / denominator) ** 2
        + ((1.0 / numerator - 1.0 / denominator) * 1.0) ** 2
    )
    expected_sigma_variance = (
        tau_variance / 10.0**2 + (tau * 2.0 / 10.0**2) ** 2
    )
    assert uncertainty.optical_depth_std[0] == pytest.approx(
        np.sqrt(tau_variance)
    )
    assert uncertainty.surface_density_std[0] == pytest.approx(
        np.sqrt(expected_sigma_variance)
    )


def test_surface_density_conversion_is_physical() -> None:
    assert G_PER_CM2_TO_MSUN_PER_PC2 == pytest.approx(4788.45, rel=1e-4)
    assert convert_surface_density(np.array([1.0]))[0] == pytest.approx(
        G_PER_CM2_TO_MSUN_PER_PC2
    )
