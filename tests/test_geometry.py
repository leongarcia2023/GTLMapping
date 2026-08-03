from __future__ import annotations

import numpy as np
import pytest
from astropy.wcs import WCS

from gtlmapping.exceptions import GridMismatchError
from gtlmapping.geometry import ellipse_mask, validate_compatible_grids
from gtlmapping.models import CloudEllipse


def test_ellipse_position_angle_is_wcs_aware(galactic_wcs: WCS) -> None:
    cloud = CloudEllipse(
        name="G028.37+00.07",
        component="0",
        glon_deg=28.373,
        glat_deg=0.076,
        major_axis_arcmin=0.8,
        minor_axis_arcmin=0.4,
        pa_deg=0.0,
    )
    mask = ellipse_mask((61, 61), galactic_wcs, cloud)
    rows, cols = np.where(mask)

    assert mask[30, 30]
    assert np.ptp(rows) > np.ptp(cols)


def test_grid_validator_detects_subtle_wcs_shift(galactic_wcs: WCS) -> None:
    shifted = galactic_wcs.deepcopy()
    shifted.wcs.crpix[0] += 1.0

    with pytest.raises(GridMismatchError, match="misaligned"):
        validate_compatible_grids(
            (61, 61),
            galactic_wcs,
            (61, 61),
            shifted,
            tolerance_pixels=0.1,
        )


def test_grid_validator_accepts_identical_grids(galactic_wcs: WCS) -> None:
    diagnostics = validate_compatible_grids(
        (61, 61),
        galactic_wcs,
        (61, 61),
        galactic_wcs.deepcopy(),
    )
    assert diagnostics["max_separation_pixels"] == pytest.approx(0.0)
