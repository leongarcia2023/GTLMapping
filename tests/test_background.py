from __future__ import annotations

import numpy as np
from astropy.wcs import WCS

from gtlmapping.background import (
    estimate_box_background,
    estimate_smf_background,
    measure_box_background,
)
from gtlmapping.geometry import ellipse_mask
from gtlmapping.models import CloudEllipse


def test_smf_recovers_smooth_background_inside_cloud(galactic_wcs: WCS) -> None:
    yy, xx = np.indices((61, 61), dtype=float)
    truth = 100.0 + 0.05 * xx + 0.02 * yy
    cloud = CloudEllipse(
        name="G028.37+00.07",
        component="0",
        glon_deg=28.373,
        glat_deg=0.076,
        major_axis_arcmin=0.7,
        minor_axis_arcmin=0.45,
        pa_deg=30.0,
    )
    mask = ellipse_mask(truth.shape, galactic_wcs, cloud)
    observed = truth.copy()
    observed[mask] -= 30.0

    result = estimate_smf_background(
        observed,
        galactic_wcs,
        cloud,
        sampling_arcsec=5.0,
    )
    error = np.abs(result.values[mask] - truth[mask])

    assert np.all(np.isfinite(result.values))
    assert np.median(error) < 1.0
    assert result.diagnostics["sample_count"] > 3


def test_adjacent_box_background_uses_mean_of_clipped_medians() -> None:
    image = np.zeros((10, 10), dtype=float)
    image[0:3, 0:3] = 5.0
    image[7:10, 7:10] = 9.0
    image[8, 8] = 100.0

    result = estimate_box_background(
        image,
        [(0, 3, 0, 3), (7, 10, 7, 10)],
        maximum_intensity=15.0,
    )

    assert result.method == "boxes"
    assert np.all(result.values == 7.0)
    assert result.diagnostics["box_medians"] == [5.0, 9.0]

    level, diagnostics = measure_box_background(
        image,
        [(0, 3, 0, 3), (7, 10, 7, 10)],
        maximum_intensity=15.0,
    )
    assert level == 7.0
    assert diagnostics["box_median_standard_error"] == 2.0
