# Changelog

## Unreleased

- Add BSD 3-Clause licensing and standards-compliant SPDX package metadata.
- Add canonical GitHub project links, artifact validation, and clean-wheel
  smoke testing to the release workflow.
- Keep the packaging-consistency tests compatible with supported Python 3.10.

## 0.4.0

- Add a documented ``moderate`` foreground profile between conservative and
  liberal GTL.
- Default moderate GTL to a 50% soft BT12 anchor, 0.5% near-saturation budget,
  and 0.01% strict lower-limit ceiling.
- Add ``GTLMapper.compute_moderate`` and matching command-line behavior so its
  small censored set stays finite and explicitly flagged.
- Clarify that magenta validation overlays are lower-limit markers and that the
  plotting color ceiling can make distinct high-column pixels share a color.

## 0.3.0

- Add an explicit ``liberal`` foreground option alongside the unchanged
  ``conservative`` default.
- Fit the liberal foreground from a robust quadratic GTL sample trend with no
  BT12 floor and an optional, visible soft BT12 anchor weight.
- Control near-saturated and strictly saturated pixels with separate budgets
  so a more permissive map cannot silently become raw, unconstrained kriging.
- Add ``GTLMapper.compute_liberal`` to project the foreground onto the
  foreground/background feasible range and report censored pixels as finite
  lower limits rather than NaNs.
- Add matching Python and command-line options, diagnostics, tests, and
  scientific-interpretation documentation.

## 0.2.3

- Restore the hard BT12 foreground floor explicitly present in the
  December 2025 prototype history.
- Replace the signed, median-centered spatial correction with
  `maximum(spatial foreground, BT12 foreground)`, so GTL cannot lower
  foreground or surface density on a jointly valid BT12 pixel.
- Permit no new strictly saturated pixels in the fitted region by
  default, preventing detected-mass decreases caused by newly masked
  holes.
- Add pointwise foreground and surface-density monotonicity tests.
- Revalidate Cloud C and Sgr C: summed surface density increases by
  3.91% and 2.95%, respectively, while both morphologies remain intact.
- Distinguish the BT12 global-minimum selection count from the common
  map-level within-two-sigma count.
- Preserve the BT12 anchor when explicitly constraining a conservative
  foreground against a background model.

## 0.2.2

- Withdraw the aphysical direct-kriging Sgr C mass comparison from
  0.2.1 after it erased the filament and produced thousands of
  artificial saturation classifications.
- Make a BT12-anchored robust planar trend the high-level foreground
  default. Local minima now constrain a broad gradient rather than the
  absolute foreground surface.
- Add explicit local- and strict-saturation budgets that automatically
  regularize the spatial trend.
- Retain raw kriging and the other interpolators only as explicit
  experimental and notebook-compatibility methods.
- Replace the Sgr C validation products and numerical conclusions with
  the morphology-preserving conservative result.

The signed-centered anchoring in this release was superseded in 0.2.3
because it allowed the spatial foreground and mass to fall below BT12.

## 0.2.1

- Add a reproducible, controlled Sgr C F480M comparison of the BT12
  and spatially varying GTL foreground methods.
- Add low-memory adjacent-box background measurement with per-box
  medians, scatter, and standard error.
- Record background, grid, and source-mask sensitivity tests and
  clarify independent-sample versus map-level saturation counts.
- Confirm the OH94 F480M and IRAC-4 filter-convolved opacity values
  against Rubén Fedriani's supplied figure.

The direct-kriging scientific comparison in this release was
superseded and withdrawn in 0.2.2.

## 0.2.0

- Add an explicit BT12 foreground estimator using the cloud-wide
  saturated-pixel mean minus two noise standard deviations.
- Restore exact notebook-compatible GTL scan origins and full-image
  detection defaults.
- Preserve overlapping-window multiplicity in kriging and weighted
  interpolation.
- Correct the Cloud C validation against the intended `1kx1k.fits`
  input and its exactly aligned `SMFbg1.fits` background.
- Resolve `1kx1k.fits` as an exact slice of a BUNIT-declared GLIMPSE
  IRAC-4 MJy/sr mosaic.
- Stabilize kriging with unique sample coordinates, a pseudo-inverse,
  finite-gap recovery, and explicit notebook repeated-row mode.
- Add JWST adjacent-box backgrounds and direct SCI/ERR FITS loading.
- Add OH94/WD01 filter-convolved opacity lookup with gas-to-dust
  normalization.
- Add saturation lower limits, explicit physical foreground
  constraints, and first-order uncertainty products.
- Audit the supplied Sgr C data, Rubén correspondence, manuscript, and
  MIREX slides.
- Expand the scientific-lineage and method-boundary documentation.

## 0.1.0

- Extract the prototype notebook into a modular package.
- Add WCS-aware Simon catalog parsing and ellipse masks.
- Add spatially varying foreground detection and six interpolation
  modes.
- Add BT09-style LMF and SMF background estimators.
- Add strict grid validation and optional reprojection.
- Add masked radiative transfer and multi-extension FITS output.
- Add tests, examples, Sphinx documentation, and release tooling.
