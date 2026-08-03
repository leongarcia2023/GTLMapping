# GTLMapping validation report

## Scope

This audit covers the prototype notebook, `1kx1k.fits`, the Simon
catalog, the prepared Cloud C background and output products, the
extinction-mapping paper lineage listed in
`extinction-mapping-papers.md`, and the package implementation. It is a
software/reproducibility audit, not a claim that the new spatially
varying foreground estimator is scientifically validated across an
IRDC population.

The local ChatGPT project mirror contained no synced files under
`sources/`. Project history was therefore reconstructed from the
notebook copies, scripts, FITS products, figures, the CASSUM project
brief, papers in `~/Downloads`, and the local IPython command history.
The December 22, 2025 history explicitly defines
`FG_floored = np.maximum(recursiveinterpfg, I_fg_BT12_map)` and then
checks the GTL/BT12 ratio from a minimum of one. It also records the
intended invariant that, with a shared background, higher foreground
must not produce lower surface density.

## Input inventory

- `1kx1k.fits` is present and is the notebook's intended observed
  image. Despite the filename, it is 918 by 918 pixels, with a
  0.6-arcsec Galactic WCS centered near
  (l, b) = (28.37317, 0.07583) degrees.
- It is an exact pixel-for-pixel slice of `TestIRDCC.fits`, rows
  831:1749 and columns 709:1627. The parent is the GLIMPSE
  `GLM_02850+0000_mosaic_I4.fits` IRAC-4 mosaic and declares
  `BUNIT=MJy/sr`; the clipping step retained the WCS but dropped that
  keyword. The intensity-unit provenance is therefore resolved.
- Its finite intensity range is 29.6751 to 5484.7144 with a median of
  74.8497 MJy/sr.
- `SMFbg1.fits` has the same shape and an exactly matching tested WCS:
  maximum center/corner separation is 0.0 arcsec.
- The local copies of `catalog.dat` and the decompressed catalog are
  byte-identical. The parser reads all 23,705 records and finds Cloud C
  at G028.37+00.07.

## Notebook reproduction

The package now reproduces the historical GTL detection contract on
`1kx1k.fits`:

| Quantity | Notebook/package result |
| --- | ---: |
| Scan windows | 81 |
| Accepted windows | 48 |
| Rejected windows | 33 |
| Raw stored minima | 48 |
| Unique coordinates | 29 |
| Raw sample min/max | 29.6751 / 66.4471 MJy/sr |

Overlapping windows repeatedly select some coordinates. The package
stores 29 unique coordinates but restores their multiplicity when
kriging, yielding the notebook's 48-row numerical fit. With the
prototype settings (`foreground_margin=1.0`,
`clip_to_sample_range=False`), the raw kriging prediction is:

- minimum: 21.8058835 MJy/sr;
- median: 48.9535671 MJy/sr; and
- maximum: 70.4255307 MJy/sr.

Those values agree with a direct execution of the notebook algorithm
to numerical precision.

Compared with the saved `GTLCloudCSigmaSMF.fits` product, the
reproduced legacy calculation over 839,833 jointly valid pixels has:

- mean absolute difference: 0.02667 g/cm²;
- bias: +0.01842 g/cm²;
- RMSE: 0.03821 g/cm²; and
- Pearson correlation: 0.97974.

The saved file is therefore closely related but not bit-for-bit
reproducible from the visible notebook state. Likely causes include
prior in-memory state, a different intermediate/background revision,
or a different execution order.

## BT12 audit

The published BT12 foreground prescription is:

1. find the global minimum inside the IRDC;
2. find every pixel strictly within `minimum + 2 sigma`;
3. require at least one qualifying pixel 8 arcsec or more from the
   minimum;
4. take the mean of all qualifying saturated pixels; and
5. subtract `2 sigma` from that mean.

On `1kx1k.fits` inside the Simon Cloud C ellipse:

- global minimum: 29.6750946 MJy/sr at row 361, column 635;
- qualifying saturated pixels: 53;
- qualifying pixels at least 8 arcsec from the minimum: 9;
- saturated-pixel mean: 30.5360785 MJy/sr; and
- BT12 foreground: 29.3360785 MJy/sr.

The notebook's comparison function instead returns the median of only
the spatially distant candidates and its map cell then adds
1.0 MJy/sr. For this image that gives 31.5882072 MJy/sr. Both the
statistic and sign differ from BT12. The package now has a distinct
`method="bt12"` path with the published mean-minus-`2 sigma`
behavior.

With `SMFbg1.fits`, the corrected BT12 surface-density map has:

- min/median/max: -0.61723 / 0.03465 / 0.70499 g/cm²; and
- valid fraction: 100%.

## Legacy kriging pressure test

The pre-0.2.2 kriging path used the paper-motivated
`2 sigma = 1.2 MJy/sr` margin, unique-coordinate pseudo-inverse
kriging, and clipping to the observed GTL sample range before applying
that margin. On the same observed/background pair:

- foreground min/median/max:
  28.47509 / 48.33781 / 65.24708 MJy/sr;
- surface-density min/median/max:
  -0.74411 / 0.05428 / 1.74017 g/cm²; and
- valid fraction: 99.6461%.

These values validate historical numerical reproduction, not the
physical foreground. Direct interpolation of local minima is no longer
the package default. Exact prototype compatibility remains available
through
`kriging_duplicate_policy="repeat"`,
`foreground_margin=1.0`, and `clip_to_sample_range=False`.

Leave-one-unique-location-out foreground errors are:

| Method | MAE | RMSE | Bias | Max absolute error |
| --- | ---: | ---: | ---: | ---: |
| Gaussian (200 px) | 9.52 | 11.68 | -4.38 | 24.07 |
| Cauchy (500 px) | 10.58 | 12.56 | -3.04 | 25.83 |
| RBF | 5.35 | 7.31 | -0.25 | 17.64 |
| Spline | 5.93 | 7.58 | +0.15 | 21.28 |
| Kriging, old direct inverse | 6.53 | 12.55 | -2.19 | 60.50 |
| Kriging, stable aggregate | 5.58 | 8.02 | -0.68 | 30.83 |
| Kriging, stable repeated rows | 5.51 | 8.01 | -1.17 | 30.83 |

The pseudo-inverse removes the worst direct-inverse instability, but
these results still do not justify declaring kriging physically
superior. Method selection and hyperparameters need multi-cloud
validation.

The corrected conservative default fits a robust plane to the local
samples and applies the historical pointwise BT12 floor. On Cloud C,
the saturation guardrail selects a blend factor of 0.17. Inside the
Simon ellipse the foreground min/median/max is
29.3361 / 30.7197 / 33.9163 MJy/sr, compared with the constant
29.3361-MJy/sr BT12 value. The common within-2-sigma count increases
from 22 to 68, while the strict-saturation count remains zero. All
761,116 ellipse pixels remain valid and satisfy
`Sigma_GTL >= Sigma_BT12`; the summed surface density increases by
3.91%. Unlike raw interpolation, the foreground contains only a broad
gradient and does not reproduce Cloud C's filament morphology.

## Package mistakes found and corrected

1. The prior `flat` interpolation was described as though it were
   BT12. It was not. BT12 is now an explicit estimator.
2. Selecting a cloud caused GTL detection to be restricted to the
   Simon ellipse. The notebook scans the full image. Full-image search
   is restored as the default; cloud restriction is opt-in.
3. The scan rounded a half-box stride and always appended edge windows,
   changing the 81-window notebook sample set. The default now uses the
   notebook's floored half-box stride and origins. Edge completion is
   opt-in.
4. Deduplication discarded the statistical effect of repeated
   overlapping-window detections. Multiplicity is now restored during
   kriging and weighted interpolation.
5. The previous report compared `CloudC550.fits` with `SMFbg1.fits`
   and attributed the resulting WCS mismatch to the intended notebook
   pair. `1kx1k.fits` and `SMFbg1.fits` are exactly aligned.
6. Documentation overclaimed equivalence between the package's
   bright-pixel policies and BT09 source cleaning. That boundary is now
   explicit.
7. The first conservative replacement centered a signed plane on BT12.
   That allowed half of the spatial model to fall below BT12 and reduced
   the Sgr C mass. The recovered project history instead used BT12 as a
   hard pointwise floor. The implementation now enforces that invariant
   and permits no new strict-saturation holes by default.

## Literature audit and method boundaries

- BT09 establishes the radiative-transfer equation, the LMF/SMF
  background framework, the fiducial 8-micron opacity, an approximately
  10% background-estimation floor, and special treatment of bright
  sources.
- BT12 establishes the empirical constant foreground from independent
  saturated cores using the mean-minus-`2 sigma` rule and an 8-arcsec
  independence criterion.
- Kainulainen & Tan (2013) combines NIR and MIR maps, adopts
  `tau_8 = 0.29 tau_K`, and estimates roughly 30% absolute MIR-opacity
  uncertainty.
- Butler, Tan & Kainulainen (2014) applies local saturation searches to
  deeper mosaics and uses a 4-arcsec spatial criterion for those data.
- Lim & Tan (2014) and Lim et al. (2016) show that the relevant noise,
  beam, foreground, background, and opacity choices are
  wavelength/dataset dependent and can materially alter derived
  distributions.
- Fedriani et al. (2025) demonstrates a JWST BT12-style implementation
  using the error extension for sigma, adjacent boxes for background,
  and a filter-specific opacity.
- André et al. (2025) demonstrates both a median-background path and
  calibration against an independent submillimeter map, while
  explicitly checking saturation.
- Fahrion & De Marchi (2023) is complementary red-clump extinction-law
  mapping, not a diffuse-background MIREX foreground algorithm.

The current package includes the supplied filter-convolved opacity
table and first-order uncertainty propagation. It does not yet
convolve arbitrary throughput/model tables, merge NIR/MIR maps, infer
spatially varying extinction laws, automate point/extended-source
masks, model correlated systematics, or perform censored Monte Carlo
inference. Those are documented limitations, not silently implied
features.

## JWST Sgr C / Rubén audit

- `F480M_registered.fits` and `F480M_registered (1).fits` are
  byte-identical. Their `SCI` and `ERR` extensions are aligned
  5660-by-2300 MJy/sr images at 0.06293 arcsec/pixel.
- Rubén's supplied filament coordinates map to the NumPy slice
  `[1980:2103, 77:200]`.
- The cutout's ERR median is 0.1023237 MJy/sr, agreeing with the
  manuscript's approximately 0.10 MJy/sr.
- The paper-faithful BT12 run finds 129 saturated pixels, 49
  independently separated pixels at 0.74 arcsec, and a constant
  foreground of 2.3053355 MJy/sr.
- An 8-by-8 edge-complete GTL scan returns 56 raw detections and 49
  unique locations. Their 2.399–6.064 MJy/sr range demonstrates that
  the local-window test does not prove every minimum is a foreground
  measurement.
- Five explicitly documented target-sized touching boxes give a shared
  comparison background of 7.2666511 MJy/sr. They are a reproducible
  pressure-test choice, not a reconstruction of Rubén's unpublished
  final box coordinates.
- Under the same background, F480M opacity (9.76 cm²/g at gas/dust
  156), 8.15-kpc distance, aperture, pixel area, and bright-pixel
  policy, the full 123-by-123 box contains 30.35 solar masses with
  BT12 and 31.24 solar masses with the hard-floor conservative GTL
  foreground.
- The GTL/BT12 mass ratio remains 1.026–1.036 when the background is
  varied over the one-standard-deviation scatter among box medians.
  Across `grid_n=4` through 10, the conservative GTL mass is
  31.16–31.49 solar masses and every ratio is 1.027–1.038; after
  three- and six-pixel bright-source-mask dilations the ratios are
  1.031 and 1.033.
- The first direct-kriging comparison failed its physical sanity
  check. It made 4,387 pixels locally saturation-consistent, made
  3,089 strict lower limits, and erased the filament morphology by
  reproducing a cloud-shaped trough in the foreground. Its
  68.26-solar-mass sum is withdrawn.
- The final replacement fits only a robust planar gradient and applies
  BT12 as a hard pointwise floor. Saturation budgets select a blend
  factor of 0.03. Under a common map-level definition, the
  within-2-sigma count increases from 55 for BT12 to 117 for GTL,
  while the strict-saturation count remains zero over every tested
  grid. Every one of the 15,129 Sgr C pixels has a nonnegative
  GTL-minus-BT12 surface-density difference, and the filament remains
  visible.
- A finite foreground does not guarantee a finite extinction solution.
  Where the foreground exceeds the separately estimated off-cloud
  background, the inputs violate the radiative-transfer ordering.
  GTLMapping now distinguishes interpolation gaps, saturated pixels,
  invalid-background pixels, and explicitly constrained pixels.
- Exact Rubén mass reproduction remains pending because the
  manuscript's adjacent background boxes are graphical, not supplied
  as machine-readable regions. No background-box coordinates were
  inferred silently. Visual inspection also finds residual stellar/PSF
  structure, so a publication result needs the final source mask.

The F480M opacity discrepancy is resolved: 9.76 cm²/g at gas/dust 156
and 15.2256 cm²/g at gas/dust 100 are the same OH94 moderately
coagulated thin-ice filter convolution under different total-mass
normalizations. The draft manuscript mixes the two normalizations in
different sections, so the package always records the requested ratio.

## Slide audit

Slides 7 and 17 of the supplied MIREX presentation invert the
radiative-transfer fraction. The notebook, package, and Sgr C
manuscript use the correct equation:

`tau = -ln[(I_obs - I_fg) / (I_bg - I_fg)]`.

The slide-level claim that the new foreground adds about 100,000 solar
masses to Cloud C is not yet pressure-tested under a fixed aperture,
distance, NIR correction, source mask, opacity normalization, and
uncertainty model. It should remain a hypothesis rather than a package
validation result.

## Cloud F and H profile portability check

- `CloudF.fits` and `IRDCCloudH.fits` contain the complete WCS-aware Simon
  ellipses for `G034.43+00.24` and `G035.39-00.33`. The older `CloudH.fits`
  footprint clips the H ellipse and is excluded from integrated results.
- The default liberal foreground uses a robust quadratic GTL-sample trend,
  zero BT12 anchor weight, a 1% local-saturation budget, and a 0.1% strict
  lower-limit budget.
- Cloud F changes from 36 BT12-selected pixels to 833 locally
  saturation-consistent pixels, of which 81 are strict lower limits. Its
  nonnegative summed surface-density proxy increases by 16.61% relative to
  BT12.
- Cloud H changes from 53 BT12-selected pixels to 616 locally
  saturation-consistent pixels, of which 259 are strict lower limits. Its
  nonnegative summed surface-density proxy increases by 25.56% relative to
  BT12.
- Moderate GTL reduces the result to 268 local/9 strict pixels for Cloud F and
  112 local/26 strict pixels for Cloud H. Its nonnegative summed
  surface-density proxies are 9.89% and 8.38% above BT12, compared with
  1.38%/2.26% for conservative and 16.61%/25.56% for liberal.
- Hollow magenta validation markers identify finite censored lower limits, not
  NaNs. The displayed surface-density color scale is capped at 0.5 g/cm², so
  different values above the ceiling can share the same red/pink color even
  though their FITS values differ.
- All finite input pixels inside both ellipses have finite liberal surface
  densities and neither ellipse contains an invalid-background pixel. The
  filamentary structures remain visible. Strict pixels are explicitly marked
  as censored lower limits; they are not reclassified as ordinary detections.
- These checks establish portability and numerical behavior only. The liberal
  saturation budgets remain explicit modeling assumptions requiring a
  sensitivity study and independent physical validation.

## Automated coverage

The test suite covers catalog parsing, WCS ellipse orientation, grid
mismatch detection, notebook-compatible and edge-complete scan modes,
duplicate multiplicity, the paper-faithful BT12 estimator, conservative and
liberal foreground contracts, all primary GTL interpolators, stable
pseudo-inverse kriging, synthetic SMF and
adjacent-box background recovery, filter/gas-dust opacity scaling,
radiative transfer, saturation lower limits, first-order uncertainty,
physical-unit conversion, and auditable FITS output.

Passing these tests validates software contracts. A public scientific
release still needs domain review, more clouds, independent comparison
data, uncertainty propagation, a source-masking policy, and
confirmation of the input intensity units.
