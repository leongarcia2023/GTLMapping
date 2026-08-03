from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS


@pytest.fixture
def galactic_wcs() -> WCS:
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [31.0, 31.0]
    wcs.wcs.crval = [28.373, 0.076]
    wcs.wcs.cdelt = np.array([-1.0 / 3600.0, 1.0 / 3600.0])
    wcs.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
    return wcs


@pytest.fixture
def galactic_header(galactic_wcs: WCS) -> fits.Header:
    return galactic_wcs.to_header()
