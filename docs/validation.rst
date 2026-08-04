Validation
==========

Software tests
--------------

The automated suite checks four parts of the package:

.. list-table::
   :header-rows: 1

   * - Area
     - Tests
   * - Data and geometry
     - Simon catalog parsing, cloud lookup, WCS ellipses, grid mismatch, and FITS output
   * - Foreground and background
     - BT12, the three GTL profiles, legacy interpolation, SMF, LMF, and adjacent boxes
   * - Radiative transfer
     - Opacity scaling, masks, lower limits, bright-pixel policies, and physical units
   * - Uncertainty
     - ERR input, foreground variance, background error, and opacity error

Run the suite from the repository root:

.. code-block:: console

   python -m pytest

The tests establish that the implementation follows its documented numerical
rules. The cases below ask whether those rules behave sensibly on real cloud
images.

Cloud C regression case
-----------------------

The ``1kx1k.fits`` image reproduces the prototype notebook's 81 windows, 48
accepted minima, 33 rejected windows, and 29 unique foreground coordinates.
The notebook-compatible kriging surface agrees with a direct execution to
numerical precision.

Inside the Simon Cloud C ellipse, the corrected BT12 routine finds 53
near-minimum pixels and a foreground of 29.3361 MJy sr\ :sup:`-1`. The
prototype's comparison cell used a median and the wrong offset sign, which
explains the earlier disagreement.

Cross-validation remains demanding. Leave-one-location-out mean absolute
errors range from about 5.4 to 10.6 MJy sr\ :sup:`-1` across the interpolation
methods. Stable kriging reduces the worst historical fold from 60.5 to 30.8
MJy sr\ :sup:`-1`, but the remaining errors rule out treating a numerically
stable surface as a validated physical foreground.

The conservative foreground inside the ellipse has a min / median / max of
29.3361 / 30.7197 / 33.9163 MJy sr\ :sup:`-1`, compared with the constant
BT12 value of 29.3361 MJy sr\ :sup:`-1`. It increases the
within-:math:`2\sigma` count from 22 to 68 and the summed surface density by
3.91%. It creates no new strict saturation, and all 761,116 ellipse pixels
satisfy :math:`\Sigma_\mathrm{GTL}\geq\Sigma_\mathrm{BT12}`.

The input image and ``SMFbg1.fits`` share the same grid. ``1kx1k.fits`` is an
exact slice of its parent GLIMPSE IRAC-4 mosaic, whose FITS header gives
MJy/sr as the intensity unit. ``VALIDATION_REPORT.md`` contains the complete
numerical record.

Sgr C F480M case
----------------

The F480M test loads the SCI and ERR extensions together and reproduces the
BT12 foreground in a 123 by 123 pixel filament cutout. The median ERR is
0.1023237 MJy sr\ :sup:`-1`. The BT12 selection has 129 pixels, 49 of which
meet the 0.74 arcsec independence criterion, and gives a foreground of
2.3053355 MJy sr\ :sup:`-1`.

Raw kriging fails the morphology check: 3,089 pixels become strict lower
limits and the fitted foreground copies the filament. The conservative model
retains the filament, raises the shared within-:math:`2\sigma` count from 55
to 117, and creates no strict lower-limit holes. With one background and mask
policy, its mass is 31.24 M\ :sub:`sun` versus 30.35 M\ :sub:`sun` for BT12.

The :doc:`jwst_sgrc` page gives the comparison boxes and sensitivity tests.
The final publication mass remains pending the author-defined aperture and
background regions.

Cloud F and H portability cases
-------------------------------

The Simon catalog entries ``G034.43+00.24`` and ``G035.39-00.33`` load Cloud
F and Cloud H without hand-entered ellipses. ``CloudF.fits`` covers the full F
ellipse. ``IRDCCloudH.fits`` covers the full H ellipse; the older
``CloudH.fits`` clips it and is excluded from integrated comparisons.

.. list-table::
   :header-rows: 1

   * - Profile
     - Cloud F increase over BT12
     - Cloud H increase over BT12
   * - Conservative
     - 1.38%
     - 2.26%
   * - Moderate
     - 9.89%
     - 8.38%
   * - Liberal
     - 16.61%
     - 25.56%

Moderate GTL produces 9 strict lower limits in Cloud F and 26 in Cloud H.
Liberal GTL produces 81 and 259. Every finite input pixel inside both ellipses
still receives a finite surface-density value; the lower-limit mask preserves
which values are censored.

The validation plots use hollow magenta symbols for finite lower limits. Their
surface-density color scale ends at 0.5 g cm\ :sup:`-2`, so pixels above that
value share one display color even though the FITS file keeps their distinct
values.

These cases show that the same catalog and mapping interface works across
three clouds. They do not establish that the default saturation budgets are
optimal. Publication work should vary the budgets and compare the results
with an independent column-density tracer.

Scientific boundary
-------------------

Passing the test suite establishes file handling, units, masks, and numerical
behavior. It cannot show that an interpolated foreground is an unbiased
physical estimate. That claim requires more clouds, independent data, and
review by the scientific team.
