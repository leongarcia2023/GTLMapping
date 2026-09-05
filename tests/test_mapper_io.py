from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from gtlmapping.mapper import GTLMapper
from gtlmapping.models import InterpolationResult


def test_mapper_writes_auditable_multi_extension_fits(
    tmp_path: Path,
    galactic_header: fits.Header,
    galactic_wcs: WCS,
) -> None:
    observed = np.full((61, 61), 60.0)
    mapper = GTLMapper(observed, header=galactic_header, wcs=galactic_wcs)
    mapper.foreground_result = InterpolationResult(
        values=np.full_like(observed, 20.0),
        method="flat",
    )
    mapper.set_background(np.full_like(observed, 100.0))
    result = mapper.compute()
    output = result.write(tmp_path / "mapping.fits")

    with fits.open(output) as hdul:
        assert [hdu.name for hdu in hdul] == [
            "SIGMA",
            "TAU",
            "FOREGROUND",
            "BACKGROUND",
            "SATURATED",
            "INVALID_BG",
            "BRIGHT",
            "UNRESOLVED",
        ]
        assert hdul[0].header["BUNIT"] == "g cm-2"
        assert hdul[0].data[0, 0] > 0


def test_mapper_loads_error_hdu_and_writes_uncertainty_products(
    tmp_path: Path,
    galactic_header: fits.Header,
) -> None:
    source = tmp_path / "jwst_like.fits"
    observed = np.full((8, 8), 6.0)
    error = np.full((8, 8), 0.1)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(observed, header=galactic_header, name="SCI"),
            fits.ImageHDU(error, header=galactic_header, name="ERR"),
        ]
    ).writeto(source)

    mapper = GTLMapper.from_fits(source, hdu="SCI", uncertainty_hdu="ERR")
    mapper.foreground_result = InterpolationResult(
        values=np.full_like(observed, 2.0),
        method="bt12",
    )
    mapper.set_background(np.full_like(observed, 10.0), uncertainty=0.2)
    result = mapper.compute(filter_name="F480M", kappa_std_cm2_g=1.0)
    output = result.write(tmp_path / "mapping_with_uncertainty.fits")

    assert result.kappa_cm2_g == 9.76
    assert result.uncertainty is not None
    with fits.open(output) as hdul:
        assert "SIGMA_ERR" in hdul
        assert "TAU_ERR" in hdul
        assert np.all(np.isfinite(hdul["SIGMA_ERR"].data))


def test_mapper_loads_jwst_error_hdu_without_repeated_wcs(
    tmp_path: Path,
    galactic_header: fits.Header,
) -> None:
    source = tmp_path / "jwst_i2d_like.fits"
    observed = np.full((8, 8), 6.0)
    error = np.full((8, 8), 0.1)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(observed, header=galactic_header, name="SCI"),
            fits.ImageHDU(error, name="ERR"),
        ]
    ).writeto(source)

    mapper = GTLMapper.from_fits(source, hdu="SCI", uncertainty_hdu="ERR")

    assert np.array_equal(mapper.observed, observed)
    assert np.array_equal(mapper.observed_std, error)
    assert mapper.wcs.has_celestial


def test_foreground_constraint_is_explicit_and_removes_invalid_background(
    galactic_header: fits.Header,
) -> None:
    observed = np.array([[5.0, 6.0]])
    mapper = GTLMapper(observed, header=galactic_header)
    mapper.foreground_result = InterpolationResult(
        values=np.array([[7.0, 3.0]]),
        method="kriging",
    )
    mapper.set_background(np.array([[8.0, 8.0]]))

    constrained = mapper.constrain_foreground(
        minimum_transmitted_intensity=2.0,
    )
    result = mapper.compute(
        saturation_policy="lower_limit",
        intensity_floor=1.0,
    )

    assert constrained.values[0, 0] == 6.0
    assert constrained.constraint_mask[0, 0]
    assert not result.invalid_background_mask.any()
    assert not result.surface_density.mask.any()


def test_foreground_constraint_preserves_bt12_anchor(
    galactic_header: fits.Header,
) -> None:
    observed = np.array([[5.0, 6.0]])
    mapper = GTLMapper(observed, header=galactic_header)
    mapper.foreground_result = InterpolationResult(
        values=np.array([[4.5, 4.5]]),
        method="conservative",
        diagnostics={
            "bt12_floor_enforced": True,
            "reference_foreground": 4.0,
        },
    )
    mapper.set_background(np.array([[5.0, 10.0]]))

    with np.testing.assert_raises_regex(
        ValueError,
        "incompatible with the background",
    ):
        mapper.constrain_foreground(minimum_transmitted_intensity=2.0)
