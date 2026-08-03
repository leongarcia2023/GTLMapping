Validation
==========

Automated tests
---------------

The test suite covers:

- the fixed-width Simon catalog format and Cloud C lookup;
- WCS-aware ellipses and sub-pixel grid mismatch detection;
- notebook-compatible and edge-complete aliasing scan modes;
- duplicate merging, stable pseudo-inverse kriging, and explicit
  repeated-coordinate compatibility;
- the paper-faithful BT12 constant foreground;
- the hard-BT12-floor spatial trend, monotonicity, and saturation
  guardrails;
- the sample-driven liberal trend, separate saturation budgets, and finite
  flagged lower limits;
- all primary foreground interpolators;
- synthetic SMF and adjacent-box background recovery;
- filter-convolved opacity and gas/dust scaling;
- radiative transfer, invalid masks, lower limits, and bright-pixel
  policies;
- first-order uncertainty propagation;
- the physical surface-density unit conversion; and
- multi-extension FITS output.

Run it with:

.. code-block:: console

   python -m pytest

Cloud C pressure test
---------------------

The intended ``1kx1k.fits`` image is available. The package reproduces
the notebook's 81 windows, 48 accepted minima, 33 rejected windows, and
29 unique foreground coordinates. Its raw notebook-compatible kriging
range agrees with a direct execution to numerical precision.

Inside the Simon Cloud C ellipse, the corrected BT12 implementation
finds 53 saturated pixels and a foreground of 29.3361 MJy
sr\ :sup:`-1`. The notebook's BT12 comparison cell used the wrong
statistic and offset sign.

Leave-one-location-out errors remain large—roughly 5.4 to 10.6 MJy
sr\ :sup:`-1` in mean absolute error depending on interpolation
method. Stable pseudo-inverse kriging reduces the historical worst
fold from 60.5 to 30.8 MJy sr\ :sup:`-1`, but multi-cloud validation
remains a release requirement.

The intended ``1kx1k.fits`` and ``SMFbg1.fits`` grids match exactly.
``1kx1k.fits`` is an exact slice of the parent GLIMPSE IRAC-4 mosaic,
whose header declares MJy/sr. Complete numerical results and method
boundaries are in ``VALIDATION_REPORT.md``.

With the corrected hard-floor spatial model, the Cloud C foreground
inside the Simon ellipse is 29.3361/30.7197/33.9163 MJy
sr\ :sup:`-1` (min/median/max), compared with the constant
29.3361-MJy sr\ :sup:`-1` BT12 value. The within-:math:`2\sigma`
count increases from 22 to 68, no new strictly saturated pixels are
introduced, and the summed surface density increases by 3.91 percent.
All 761,116 ellipse pixels satisfy
:math:`\Sigma_\mathrm{GTL}\geq\Sigma_\mathrm{BT12}`.

JWST Sgr C pressure test
------------------------

The package reproduces the supplied F480M SCI/ERR interface and Rubén
BT12 foreground in the 123-by-123 filament cutout. The ERR median is
0.1023237 MJy sr\ :sup:`-1`; 129 saturated pixels produce a foreground
of 2.3053355 MJy sr\ :sup:`-1`; and 49 remain independent at 0.74
arcsec.

The first direct-kriging comparison failed the physical sanity check:
it made 4,387 pixels locally saturation-consistent, made 3,089 strict
lower limits, and erased the filament morphology. That result is
withdrawn. The hard-floor conservative default retains the filament,
increases the comparable within-:math:`2\sigma` count from 55 to 117,
and introduces no strict lower-limit holes. Its mass is
31.24 M\ :sub:`sun` versus 30.35 M\ :sub:`sun` for BT12. See
:doc:`jwst_sgrc`.

The exact Rubén mass remains pending the machine-readable adjacent
background boxes and aperture. See :doc:`jwst_sgrc`.

Cloud F and H profile portability test
---------------------------------------

The Simon catalog entries ``G034.43+00.24`` (Cloud F) and
``G035.39-00.33`` (Cloud H) load without hand-entering ellipse geometry.
``CloudF.fits`` contains the complete F ellipse; ``IRDCCloudH.fits`` contains
the complete H ellipse. The older ``CloudH.fits`` footprint clips the H
ellipse and is not used for the integrated comparison.

With the default liberal budgets, Cloud F increases from 36 BT12-selected
near-minimum pixels to 833 locally saturation-consistent pixels, including 81
strictly censored lower limits. Cloud H increases from 53 BT12-selected pixels
to 616 locally saturation-consistent pixels, including 259 lower limits. The
shared-background, nonnegative summed surface-density proxies increase by
16.61% and 25.56%, respectively, relative to BT12. Both filamentary
morphologies remain visible.

The moderate profile retains 268 locally saturation-consistent pixels and 9
strict lower limits for Cloud F; for Cloud H it retains 112 and 26. Its summed
surface-density proxies are 9.89% and 8.38% above BT12, between the
conservative increases of 1.38%/2.26% and liberal increases of
16.61%/25.56%.

In the validation figures, hollow magenta symbols mark finite censored lower
limits; they are not NaNs. Independently, the surface-density color map is
displayed only through 0.5 g cm\ :sup:`-2`, so distinct values above that
ceiling share the same red/pink display color. The FITS arrays retain their
actual values.

Every finite input pixel inside both ellipses receives a finite surface
density; no invalid-background pixels are introduced. The lower-limit pixels
remain flagged rather than being treated as detections. These are portability
and numerical-contract results, not a claim that the default 1%/0.1%
saturation budgets are scientifically optimal. Those budgets require
sensitivity analysis and independent validation before publication.

Interpretation boundary
-----------------------

Passing software tests establishes that interfaces, units, masks,
file handling, and numerical contracts behave as specified. It does
not establish that a spatially varying interpolated foreground is an
unbiased physical estimator. That requires comparison with independent
data and review by the scientific team.
