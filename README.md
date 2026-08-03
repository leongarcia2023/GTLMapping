# GTLMapping

[![CI](https://github.com/leongarcia2023/GTLMapping/actions/workflows/ci.yml/badge.svg)](https://github.com/leongarcia2023/GTLMapping/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](https://github.com/leongarcia2023/GTLMapping/blob/main/LICENSE)

`GTLMapping` is an astronomy package for mid-infrared extinction
(MIREX) mapping with a spatially varying foreground. It turns the
prototype notebook into an inspectable pipeline with:

- WCS-aware Simon et al. (2006) cloud ellipses;
- local saturation detection with overlapping aliasing windows;
- an explicit paper-faithful BT12 constant-foreground mode;
- a robust spatial trend with a hard BT12 foreground floor;
- experimental kriging, RBF, spline, Gaussian, Cauchy, and flat models;
- explicit notebook-compatible kriging for historical reproduction;
- BT09-style LMF/SMF and JWST-style adjacent-box backgrounds;
- filter-convolved OH94/WD01 opacity tables with explicit gas/dust scaling;
- first-order intensity, kriging, background, and opacity uncertainties;
- explicit masks, saturation lower limits, and foreground constraints;
- multi-extension FITS output containing every modeled quantity; and
- strict shape and celestial-grid checks, with optional reprojection.

The package is alpha scientific software. Its implementation and
interfaces are tested, but the new spatially varying foreground method
still requires domain review and validation across additional IRDCs
before publication claims are made.

Three spatial profiles expose the modeling trade-off explicitly:

- `conservative` stays closely anchored to the BT12 foreground while allowing a smooth spatial fluctuation.
- `moderate` permits a stronger data-driven fluctuation without erasing cloud structure or producing new non-finite pixels.
- `liberal` gives the spatial model the most freedom and is intended for sensitivity analysis, not as an automatic default.

## Install

From this repository:

```bash
git clone https://github.com/leongarcia2023/GTLMapping.git
cd GTLMapping
python -m pip install .
```

For development and documentation:

```bash
python -m pip install -e ".[dev]"
```

After a PyPI release, the intended installation command is:

```bash
python -m pip install GTLMapping
```

## Data inputs

The source distribution does not bundle science images or a cloud catalog.
Supply the FITS image and `catalog.dat` (or an equivalent catalog path) when
running a map. This keeps observational data and machine-specific paths out of
the installable package.

## License

GTLMapping is distributed under the [BSD 3-Clause License](https://github.com/leongarcia2023/GTLMapping/blob/main/LICENSE).

## Quick start

```python
from gtlmapping import GTLMapper

mapper = GTLMapper.from_fits("1kx1k.fits")
cloud = mapper.select_cloud("catalog.dat", "G028.37+00.07")

samples = mapper.detect_foreground(noise_sigma=0.6)
mapper.fit_foreground(
    method="conservative",
    noise_sigma=0.6,
)

# Use the aligned SMF product, or call estimate_background(method="smf").
mapper.set_background_from_fits("SMFbg1.fits")

result = mapper.compute(
    kappa_cm2_g=7.5,
    bright_pixel_policy="allow",
)
result.write("cloud_c_gtl.fits")
```

The conservative default fits only a robust broad gradient to the
window minima, then applies the historical GTL rule
`foreground = maximum(spatial_foreground, BT12_foreground)`.
Consequently the spatial model cannot lower the foreground or
surface density relative to BT12 when the image, background, opacity,
and valid-pixel set are held fixed. The default guardrails also prevent
new strictly saturated pixels inside the fitted region.
Exact notebook kriging compatibility is still available explicitly
with `method="kriging"`,
`kriging_duplicate_policy="repeat"` and
`clip_to_sample_range=False`.

For an intermediate result with substantially fewer censored pixels, use the
moderate profile:

```python
samples = mapper.detect_foreground(noise_sigma=0.6)
mapper.fit_foreground(method="moderate", noise_sigma=0.6)
mapper.set_background_from_fits("smf_background.fits")
result = mapper.compute_moderate(kappa_cm2_g=7.5)
```

Moderate GTL uses a 50% soft BT12 anchor, a 0.5% near-saturation budget,
and a 0.01% strict lower-limit ceiling. These defaults sit between the
hard-floor conservative profile and the sample-dominated liberal profile.

For a controlled, more permissive interpretation, select the liberal
foreground and compute censored pixels as finite lower limits:

```python
samples = mapper.detect_foreground(noise_sigma=0.6)
mapper.fit_foreground(
    method="liberal",
    noise_sigma=0.6,
    target_local_saturation_fraction=0.01,
    maximum_strict_saturation_fraction=0.001,
)
mapper.set_background_from_fits("smf_background.fits")
result = mapper.compute_liberal(
    kappa_cm2_g=7.5,
    bright_pixel_policy="allow",
)
```

Unlike the conservative method, liberal GTL has no hard BT12 floor. It fits
a robust quadratic trend to the GTL samples and raises its absolute level only
within explicit near-saturation and strict-saturation budgets. The default
BT12 soft-anchor weight is zero. `compute_liberal()` uses a `2 sigma`
transmitted-intensity floor by default, records genuinely censored pixels in
`SATURATED`, and records foreground/background feasibility adjustments in
`FG_CONSTRAINT`. Finite values at `SATURATED` pixels are lower limits, not
ordinary measurements.

For a BT12 comparison map, replace foreground detection and
interpolation with:

```python
mapper.fit_foreground(
    method="bt12",
    cloud=cloud,
    noise_sigma=0.6,
)
```

This implements the BT12 mean saturated intensity minus `2 sigma`.
The notebook's median-plus-1.0 comparison cell is not the published
BT12 prescription.

If a background map already exists, it is checked in celestial
coordinates before use:

```python
mapper.set_background_from_fits("background.fits")
```

To reproject a mismatched map, install `GTLMapping[align]` and pass
`align=True`.

## JWST F480M

The supplied Sgr C product can load its science and uncertainty
extensions together:

```python
mapper = GTLMapper.from_fits(
    "F480M_registered.fits",
    hdu="SCI",
    uncertainty_hdu="ERR",
)

# After selecting a target mask and fitting a foreground/background:
result = mapper.compute(
    filter_name="F480M",
    gas_to_dust_ratio=156,
    kappa_std_cm2_g=0.30 * 9.76,
)
```

The built-in OH94 moderately coagulated thin-ice value is
`9.76 cm2/g` at gas-to-dust ratio 156 and `15.2256 cm2/g` at ratio
100. These are two normalizations of the same filter convolution, not
competing F480M opacities.

When a spatial foreground exceeds the independently estimated
off-cloud background, the default result masks that model conflict.
To make the projection explicit and auditable:

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

For the conservative method, this projection preserves the BT12 floor
or raises an incompatibility error; it never silently lowers the
foreground below BT12. The FITS product records adjusted pixels in
`FG_CONSTRAINT` and
the censored pixels in `SATURATED`; they should not be treated as
ordinary detections.

The controlled Sgr C comparison can be reproduced with:

```bash
python examples/sgrc_f480m_compare.py F480M_registered.fits
```

Under the documented shared comparison background, the 123-by-123
Rubén box contains 30.35 solar masses with BT12 and 31.24 solar masses
with the BT12-floored spatial foreground. Every pixel retains or
increases its BT12 surface density, the comparable within-`2 sigma`
count grows from 55 to 117, and no new strict-saturation holes are
introduced. The earlier 68.26-solar-mass direct-kriging result failed
the morphology and saturation sanity checks and is withdrawn. See the
[Sgr C audit](docs/jwst_sgrc.rst).

## Command line

```bash
gtlmapping observed.fits background.fits output.fits \
  --catalog catalog.dat \
  --cloud G028.37+00.07
```

Use `--align-background` only when reprojection is intended.

## Scientific lineage

The implementation builds on the MIREX framework of Butler & Tan
(2009, 2012), subsequent high-dynamic-range and multiwavelength
extensions, the Simon et al. (2006) MSX IRDC catalog, and the
Ossenkopf & Henning (1994) dust-opacity models. See the
[scientific method](docs/scientific_method.rst) and
[references](docs/references.rst).

## Repository status

See [VALIDATION_REPORT.md](VALIDATION_REPORT.md) for the notebook and
Cloud C pressure-test findings. See
[RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) before making the
repository public or publishing to PyPI.
