from __future__ import annotations

import numpy as np
import pytest

from gtlmapping.foreground import (
    detect_saturated_samples,
    estimate_bt12_foreground,
    fit_conservative_foreground,
    fit_liberal_foreground,
    fit_moderate_foreground,
    interpolate_foreground,
)
from gtlmapping.extinction import compute_extinction
from gtlmapping.mapper import GTLMapper
from gtlmapping.models import ForegroundSamples


def _sample_set() -> ForegroundSamples:
    return ForegroundSamples(
        rows=np.array([0, 0, 9, 9]),
        cols=np.array([0, 9, 0, 9]),
        values=np.array([10.0, 12.0, 14.0, 16.0]),
        multiplicity=np.ones(4, dtype=int),
        accepted_windows=4,
        rejected_windows=0,
        total_windows=4,
        pixel_scale_arcsec=1.0,
        min_separation_arcsec=4.0,
    )


def test_aliasing_edge_coverage_is_explicit_and_merges_duplicates() -> None:
    yy, xx = np.indices((42, 42))
    image = 100.0 + 2.0 * yy + 3.0 * xx
    image[41, 35] = 1.0
    image[41, 41] = 1.1

    notebook_samples = detect_saturated_samples(
        image,
        grid_n=2,
        overlap=0.5,
        noise_sigma=0.1,
        min_separation_arcsec=4.0,
        fallback_pixel_scale_arcsec=1.0,
    )
    edge_samples = detect_saturated_samples(
        image,
        grid_n=2,
        overlap=0.5,
        noise_sigma=0.1,
        min_separation_arcsec=4.0,
        fallback_pixel_scale_arcsec=1.0,
        cover_edges=True,
    )

    assert notebook_samples.total_windows == 9
    assert edge_samples.total_windows == 16
    assert (41, 35) not in set(
        zip(notebook_samples.rows, notebook_samples.cols, strict=True)
    )
    assert (41, 35) in set(zip(edge_samples.rows, edge_samples.cols, strict=True))
    assert len(edge_samples) <= edge_samples.raw_detection_count


def test_bt12_uses_mean_minus_two_sigma_inside_region() -> None:
    image = np.full((25, 25), 100.0)
    image[5, 5] = 10.0
    image[5, 6] = 10.4
    image[5, 15] = 10.8
    image[20, 20] = 11.0
    region = np.ones_like(image, dtype=bool)

    result = estimate_bt12_foreground(
        image,
        region_mask=region,
        noise_sigma=1.0,
        min_separation_arcsec=8.0,
        fallback_pixel_scale_arcsec=1.0,
    )

    expected = np.mean([10.4, 10.8, 11.0]) - 2.0
    assert result.method == "bt12"
    assert np.all(result.values == pytest.approx(expected))
    assert result.diagnostics["saturated_pixel_count"] == 3
    assert result.diagnostics["independent_saturated_pixel_count"] == 2


def test_conservative_foreground_uses_bt12_floor_and_limits_saturation() -> None:
    image = np.full((25, 25), 30.0)
    image[:, 13:] = 12.0
    image[2, 2] = 10.0
    image[2, 15] = 10.5
    samples = ForegroundSamples(
        rows=np.array([2, 2, 22, 22]),
        cols=np.array([2, 22, 2, 22]),
        values=np.array([10.0, 16.0, 11.0, 17.0]),
        multiplicity=np.ones(4, dtype=int),
        accepted_windows=4,
        rejected_windows=0,
        total_windows=4,
        pixel_scale_arcsec=1.0,
        min_separation_arcsec=8.0,
    )

    result = fit_conservative_foreground(
        samples,
        image,
        noise_sigma=1.0,
        min_separation_arcsec=8.0,
        fallback_pixel_scale_arcsec=1.0,
        maximum_local_saturation_fraction=0.01,
        maximum_strict_saturation_fraction=0.0,
    )

    assert result.method == "conservative"
    reference = result.diagnostics["reference_foreground"]
    assert np.all(result.values >= reference)
    assert np.any(result.values > reference)
    assert result.diagnostics["anchor_policy"] == "hard_bt12_floor"
    assert result.diagnostics["foreground_below_bt12_count"] == 0
    assert result.diagnostics["strict_saturation_count"] <= (
        result.diagnostics["strict_saturation_limit_count"]
    )
    assert result.diagnostics["local_saturation_count"] <= (
        result.diagnostics["local_saturation_limit_count"]
    )
    assert result.diagnostics["blend_factor"] < 1.0
    assert np.ptp(result.values) < np.ptp(samples.values)

    background = np.full_like(image, 40.0)
    _, sigma_bt12, _, _, _ = compute_extinction(
        image,
        background,
        np.full_like(image, reference),
        bright_pixel_policy="zero",
    )
    _, sigma_gtl, _, _, _ = compute_extinction(
        image,
        background,
        result.values,
        bright_pixel_policy="zero",
    )
    common = ~np.ma.getmaskarray(sigma_bt12) & ~np.ma.getmaskarray(sigma_gtl)
    assert np.all(sigma_gtl.data[common] >= sigma_bt12.data[common] - 1e-12)
    assert np.sum(sigma_gtl.data[common]) >= np.sum(sigma_bt12.data[common])


def _liberal_case() -> tuple[np.ndarray, ForegroundSamples]:
    yy, xx = np.indices((40, 40))
    image = 15.0 + 0.017 * yy + 0.023 * xx
    image[5, 5] = 5.0
    image[5, 20] = 5.5
    image[30, 30] = 5.8
    samples = ForegroundSamples(
        rows=np.array([2, 2, 37, 37, 20, 10]),
        cols=np.array([2, 37, 2, 37, 20, 28]),
        values=np.array([5.0, 9.0, 7.0, 11.0, 8.0, 8.5]),
        multiplicity=np.ones(6, dtype=int),
        accepted_windows=6,
        rejected_windows=0,
        total_windows=6,
        pixel_scale_arcsec=1.0,
        min_separation_arcsec=8.0,
    )
    return image, samples


def test_liberal_foreground_is_spatial_and_controls_saturation() -> None:
    image, samples = _liberal_case()
    result = fit_liberal_foreground(
        samples,
        image,
        noise_sigma=1.0,
        min_separation_arcsec=8.0,
        fallback_pixel_scale_arcsec=1.0,
    )

    diagnostics = result.diagnostics
    assert result.method == "liberal"
    assert np.all(np.isfinite(result.values))
    assert np.ptp(result.values) > 0
    assert diagnostics["anchor_policy"] == "ordered_one_sided_floor"
    assert diagnostics["bt12_anchor_weight"] == 0.0
    assert diagnostics["bt12_floor_enforced"]
    assert diagnostics["ordering_floor_method"] == "moderate"
    assert diagnostics["ordering_floor_enforced"]
    assert diagnostics["foreground_below_ordering_floor_count"] == 0
    assert diagnostics["more_local_saturation_than_bt12"]
    assert diagnostics["local_saturation_count"] > (
        diagnostics["reference_local_saturation_count"]
    )
    assert diagnostics["strict_saturation_count"] > 0
    assert diagnostics["strict_saturation_count"] <= (
        diagnostics["strict_saturation_limit_count"]
    )
    assert diagnostics["recommended_intensity_floor"] == 2.0


def test_mapper_compute_liberal_returns_finite_flagged_lower_limits(
    galactic_header,
) -> None:
    image, samples = _liberal_case()
    mapper = GTLMapper(image, header=galactic_header)
    mapper.fit_foreground(
        method="liberal",
        samples=samples,
        region_mask=np.ones_like(image, dtype=bool),
        noise_sigma=1.0,
        min_separation_arcsec=8.0,
        fallback_pixel_scale_arcsec=1.0,
    )
    mapper.set_background(np.full_like(image, 20.0))

    result = mapper.compute_liberal(bright_pixel_policy="zero")

    assert result.saturated_mask.any()
    assert not result.invalid_background_mask.any()
    assert not np.ma.getmaskarray(result.surface_density).any()
    assert np.all(np.isfinite(result.surface_density.data))
    assert mapper.foreground_result.diagnostics[
        "preserved_foreground_floor"
    ] == pytest.approx(
        mapper.foreground_result.diagnostics["reference_foreground"]
    )
    assert result.metadata["saturation_policy"] == "lower_limit"


def test_moderate_profile_reduces_censoring_and_computes_lower_limits(
    galactic_header,
) -> None:
    image, samples = _liberal_case()
    moderate = fit_moderate_foreground(
        samples,
        image,
        noise_sigma=1.0,
        min_separation_arcsec=8.0,
        fallback_pixel_scale_arcsec=1.0,
    )

    diagnostics = moderate.diagnostics
    assert moderate.method == "moderate"
    assert diagnostics["profile"] == "moderate"
    assert diagnostics["bt12_anchor_weight"] == 0.5
    assert diagnostics["target_local_saturation_fraction"] == 0.005
    assert diagnostics["maximum_strict_saturation_fraction"] == 0.0001
    assert diagnostics["ordering_floor_method"] == "conservative"
    assert diagnostics["ordering_floor_enforced"]
    assert diagnostics["foreground_below_ordering_floor_count"] == 0
    assert diagnostics["more_local_saturation_than_bt12"]
    assert 0 < diagnostics["strict_saturation_count"] <= (
        diagnostics["strict_saturation_limit_count"]
    )

    mapper = GTLMapper(image, header=galactic_header)
    mapper.foreground_result = moderate
    mapper.set_background(np.full_like(image, 20.0))
    result = mapper.compute_moderate(bright_pixel_policy="zero")

    assert result.saturated_mask.any()
    assert not result.invalid_background_mask.any()
    assert not np.ma.getmaskarray(result.surface_density).any()
    assert np.all(np.isfinite(result.surface_density.data))
    assert mapper.foreground_result.diagnostics[
        "preserved_foreground_floor"
    ] == pytest.approx(
        mapper.foreground_result.diagnostics["reference_foreground"]
    )


def test_named_profiles_are_pointwise_ordered() -> None:
    image, samples = _liberal_case()
    kwargs = {
        "noise_sigma": 1.0,
        "min_separation_arcsec": 8.0,
        "fallback_pixel_scale_arcsec": 1.0,
    }
    conservative = fit_conservative_foreground(samples, image, **kwargs)
    moderate = fit_moderate_foreground(samples, image, **kwargs)
    liberal = fit_liberal_foreground(samples, image, **kwargs)

    assert np.all(moderate.values >= conservative.values - 1e-12)
    assert np.all(liberal.values >= moderate.values - 1e-12)

    background = np.full_like(image, 20.0)
    surface_density = []
    for foreground in (
        conservative.values,
        moderate.values,
        liberal.values,
    ):
        _, sigma, _, _, _ = compute_extinction(
            image,
            background,
            foreground,
            bright_pixel_policy="zero",
        )
        surface_density.append(sigma)
    common = np.logical_and.reduce(
        [~np.ma.getmaskarray(values) for values in surface_density]
    )
    assert np.all(
        surface_density[1].data[common]
        >= surface_density[0].data[common] - 1e-12
    )
    assert np.all(
        surface_density[2].data[common]
        >= surface_density[1].data[common] - 1e-12
    )

    lower_limit_sums = []
    for foreground in (
        conservative.values,
        moderate.values,
        liberal.values,
    ):
        _, sigma, _, _, _ = compute_extinction(
            image,
            background,
            foreground,
            bright_pixel_policy="zero",
            saturation_policy="lower_limit",
            intensity_floor=2.0,
        )
        lower_limit_sums.append(
            float(np.sum(np.maximum(sigma.filled(0.0), 0.0)))
        )
    assert lower_limit_sums[0] <= lower_limit_sums[1]
    assert lower_limit_sums[1] <= lower_limit_sums[2]


@pytest.mark.parametrize("method", ["flat", "gaussian", "cauchy", "rbf", "kriging"])
def test_interpolators_are_finite_and_bounded(method: str) -> None:
    samples = _sample_set()
    result = interpolate_foreground(
        samples,
        (10, 10),
        method=method,
        foreground_margin=1.2,
        length_scale_pixels=5.0,
    )

    assert np.all(np.isfinite(result.values))
    assert np.nanmin(result.values) >= 8.8 - 1e-9
    assert np.nanmax(result.values) <= 14.8 + 1e-9
    if method == "kriging":
        assert result.diagnostics["kriging_duplicate_policy"] == "aggregate"
        assert result.diagnostics["kriging_pseudo_inverse"]
        assert result.diagnostics["remaining_interpolation_gap_count"] == 0
