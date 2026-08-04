# GTLMapping

[![CI](https://github.com/leongarcia2023/GTLMapping/actions/workflows/ci.yml/badge.svg)](https://github.com/leongarcia2023/GTLMapping/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/gtlmapping/badge/?version=latest)](https://gtlmapping.readthedocs.io/en/latest/?badge=latest)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](https://github.com/leongarcia2023/GTLMapping/blob/main/LICENSE)

GTLMapping makes mid-infrared extinction (MIREX) maps of infrared dark
clouds. It includes the constant foreground method of Butler & Tan (2012)
and three spatial models for testing how foreground structure changes the
inferred column density.

Every calculation remains inspectable. The package writes the foreground,
background, optical depth, surface density, uncertainty, and diagnostic masks
to separate extensions in one FITS file. WCS checks catch misaligned images
before they enter the radiative transfer calculation.

> GTLMapping is alpha scientific software. Its interfaces and numerical
> contracts are tested. The spatial foreground models still need validation
> on more clouds before they support general astrophysical claims.

## Installation

Install the current version from GitHub:

```bash
git clone https://github.com/leongarcia2023/GTLMapping.git
cd GTLMapping
python -m pip install .
```

Install the development tools when you want to run tests or build the
documentation:

```bash
python -m pip install -e ".[dev]"
```

The package is not yet published on PyPI.

## A first map

```python
from gtlmapping import GTLMapper

mapper = GTLMapper.from_fits("1kx1k.fits")
cloud = mapper.select_cloud("catalog.dat", "G028.37+00.07")

mapper.detect_foreground(noise_sigma=0.6)
mapper.fit_foreground(
    method="conservative",
    noise_sigma=0.6,
)
mapper.set_background_from_fits("SMFbg1.fits")

result = mapper.compute(
    kappa_cm2_g=7.5,
    bright_pixel_policy="allow",
)
result.write("cloud_c_gtl.fits")
```

The image and background must share a celestial grid. Install the `align`
extra and pass `align=True` if you want GTLMapping to reproject the
background.

## Foreground models

| Method | Foreground treatment | Best use |
|---|---|---|
| `bt12` | One value from independent saturated pixels | Published reference calculation |
| `conservative` | A broad spatial trend with BT12 as a pointwise floor | Default spatial comparison |
| `moderate` | A quadratic trend with a soft BT12 anchor and tight censoring limits | Intermediate sensitivity test |
| `liberal` | A quadratic trend driven by the samples, with no hard BT12 floor | Permissive sensitivity test |

The conservative model cannot reduce surface density relative to BT12 when
the image, background, opacity, and set of valid pixels are held fixed. Its
default guardrail also prevents new strict saturation inside the fitted
region.

Moderate and liberal fits can create pixels whose transmitted intensity is
only bounded from above. Their compute helpers assign finite lower limits and
keep those pixels marked in the `SATURATED` extension. Treat the marks as
censored measurements.

```python
mapper.detect_foreground(noise_sigma=0.6)
mapper.fit_foreground(method="moderate", noise_sigma=0.6)
mapper.set_background_from_fits("smf_background.fits")
result = mapper.compute_moderate(kappa_cm2_g=7.5)
```

Use the published BT12 calculation as the reference map:

```python
mapper.fit_foreground(
    method="bt12",
    cloud=cloud,
    noise_sigma=0.6,
)
```

Kriging, RBF, spline, Gaussian, and Cauchy interpolation remain available for
method experiments. They interpolate local minima and can imprint cloud
structure on the foreground. The main interface does not use them by default.

## JWST F480M data

Load the science and uncertainty extensions together:

```python
mapper = GTLMapper.from_fits(
    "F480M_registered.fits",
    hdu="SCI",
    uncertainty_hdu="ERR",
)

result = mapper.compute(
    filter_name="F480M",
    gas_to_dust_ratio=156,
    kappa_std_cm2_g=0.30 * 9.76,
)
```

The included OH94 opacity for moderately coagulated grains with thin ice
mantles is 9.76 cm2/g
at gas-to-dust ratio 156. The same filter convolution gives 15.2256 cm2/g
at ratio 100. Record the normalization with every result.

If a foreground estimate reaches the background outside the cloud, inspect
both surfaces. You can then impose a positive floor on transmitted intensity
and keep the adjustment in the output:

```python
mapper.constrain_foreground(
    minimum_transmitted_intensity=2 * noise_sigma,
)
result = mapper.compute(
    filter_name="F480M",
    saturation_policy="lower_limit",
    intensity_floor=2 * noise_sigma,
)
```

## Validation

The controlled Sgr C case uses one background and mask policy for both BT12
and conservative GTL. In the 123 by 123 pixel test region, the BT12 mass is
30.35 solar masses and the conservative GTL mass is 31.24 solar masses. The
filament remains visible, and the spatial model creates no new strict
saturation holes.

Run the case study with:

```bash
python examples/sgrc_f480m_compare.py F480M_registered.fits
```

The [Cloud comparison page](docs/cloud_comparisons.rst) shows the full method
on Clouds C, F, and H. The [Sgr C case study](docs/jwst_sgrc.rst) records the
fixed regions, assumptions, and sensitivity tests for F480M data.

## Data and scientific scope

Science images and `catalog.dat` stay outside the Python distribution. Supply
their paths when you run a map. This keeps large observations and local file
paths out of the package.

GTLMapping follows the MIREX work of Butler & Tan (2009, 2012), later methods
with greater dynamic range, and the dust opacities of Ossenkopf & Henning
(1994). The [method documentation](docs/scientific_method.rst) states where the
package follows those papers and where it introduces new assumptions.

See the [full documentation](https://gtlmapping.readthedocs.io/en/latest/)
and [validation page](docs/validation.rst).

## License

GTLMapping is distributed under the [BSD 3-Clause License](LICENSE).
