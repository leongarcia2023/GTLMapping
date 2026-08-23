# GTLMapping validation report

## Scope

This report records software and numerical checks for GTLMapping. It does not
claim that a spatially varying foreground is an unbiased physical estimate
for every infrared dark cloud. The current evidence covers Clouds C, F, and H,
a fixed Sgr C F480M benchmark, and synthetic tests of the package contracts.

Every comparison holds the observed image, background, opacity, aperture, and
bright pixel policy fixed while changing the foreground model. This is
essential because optical depth increases monotonically with foreground only
on pixels that remain valid in both calculations.

## Input checks

The Cloud C image is a 918 by 918 pixel cutout with a 0.6 arcsec Galactic WCS.
It is an exact slice of a GLIMPSE IRAC4 mosaic whose header records MJy/sr as
the intensity unit. The prepared `SMFbg1.fits` background has the same shape
and an exactly matching WCS at the tested center and corners.

The Simon catalog parser reads 23,705 records and identifies:

| Cloud | Catalog name | Complete image used |
| --- | --- | --- |
| C | `G028.37+00.07` | `1kx1k.fits` |
| F | `G034.43+00.24` | `CloudF.fits` |
| H | `G035.39-00.33` | `IRDCCloudH.fits` |

Cloud F and H use the package SMF background estimator. The older Cloud H
cutout clips its catalog ellipse and is excluded from integrated comparisons.

## Foreground sample detection

The package scans the full image with overlapping windows and records both
accepted windows and unique sample coordinates. Duplicate coordinates retain
their multiplicity for methods that use sample weights.

| Cloud | Accepted windows | Unique sites |
| --- | ---: | ---: |
| C | 48 | 29 |
| F | 73 | 53 |
| H | 58 | 36 |

Held out sample tests on Cloud C give mean absolute foreground errors from
about 5.4 to 10.6 MJy/sr across the optional interpolation methods. Stable
kriging improves the worst numerical fold, but the remaining errors show why
direct interpolation of local minima is not the package default.

## BT12 reference

The BT12 estimator is implemented separately from the spatial profiles. It
finds pixels within two noise standard deviations of the cloud minimum,
requires an independent pixel at the configured angular separation, takes the
mean of the qualifying pixels, and subtracts two noise standard deviations.

For Cloud C, it finds 53 pixels near the minimum and gives a foreground of
29.3361 MJy/sr. This reference value is the conservative hard floor. Moderate
and liberal retain it indirectly through the ordered profile hierarchy.

## Spatial profile comparison

The table reports the change in summed nonnegative surface density inside each
catalog ellipse. A ratio is a controlled comparison statistic, not an
independent mass measurement.

| Profile | Cloud C | Cloud F | Cloud H |
| --- | ---: | ---: | ---: |
| Conservative | +3.91% | +1.38% | +2.26% |
| Moderate | +13.15% | +9.89% | +8.39% |
| Liberal | +22.22% | +16.67% | +25.56% |

The conservative profile creates no strict lower limits in any of these
three cases. The named profiles enforce
``BT12 <= conservative <= moderate <= liberal`` pointwise, preserving the
expected surface-density ordering on all jointly valid, uncensored pixels.

Moderate GTL produces 77, 9, and 26 strict lower limits in Clouds C, F, and H.
Liberal GTL produces 762, 81, and 259. These pixels receive finite censored
lower limits and remain marked in the `SATURATED` mask. They are not NaNs and
must not be treated as ordinary detections.

The figures and machine readable metrics are published on the
[Cloud comparisons](docs/cloud_comparisons.rst) page.

## Sgr C F480M benchmark

The fixed benchmark uses a 123 by 123 pixel cutout and five adjacent
background boxes. The median ERR is 0.1023237 MJy/sr. BT12 selects 129 pixels,
49 of which meet the 0.74 arcsec independence criterion, and gives a constant
foreground of 2.3053355 MJy/sr.

With the same background, aperture, opacity, distance, and bright pixel
policy, BT12 integrates to 30.35 solar masses and conservative GTL integrates
to 31.24 solar masses. Every cutout pixel has a nonnegative GTL minus BT12
surface density difference, and neither map has a strict lower limit.

Changing the background over the measured box scatter gives GTL to BT12
ratios from 1.026 to 1.036. Changing the scan grid from 4 through 10 gives
ratios from 1.027 to 1.038. Dilating the bright source mask by three and six
pixels gives ratios of 1.031 and 1.033.

Direct kriging fails this benchmark. It creates a foreground trough shaped
like the filament and 3,089 strict lower limits. The conservative broad trend
retains the filament structure and creates none.

## Opacity normalization

The included OH94 moderately coagulated thin ice model has an F480M opacity
of 9.76 cm²/g at gas to dust ratio 156. Scaling the same total mass opacity to
ratio 100 gives 15.2256 cm²/g. These are two normalizations of the same filter
convolution, not two different dust laws.

## Automated coverage

The test suite covers catalog parsing, WCS ellipse orientation, grid mismatch
detection, scan modes, duplicate multiplicity, BT12, all three GTL profiles,
optional interpolation methods, SMF and adjacent box backgrounds, opacity
scaling, radiative transfer, finite lower limits, first order uncertainty,
physical unit conversion, and multi extension FITS output.

Passing these tests validates the documented software contracts. Broader
scientific validation still requires independent column density data, more
clouds, a source masking policy, and a sensitivity analysis for the foreground
anchor and censoring limits.
