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
     - BT12, the three GTL profiles, optional interpolation methods, SMF, LMF, and adjacent boxes
   * - Radiative transfer
     - Opacity scaling, masks, lower limits, bright pixel policies, and physical units
   * - Uncertainty
     - ERR input, foreground variance, background error, and opacity error

Run the suite from the repository root:

.. code-block:: console

   python -m pytest

These tests establish that the implementation follows its documented
numerical rules. The image checks below test whether those rules remain stable
on real clouds.

Cloud C
-------

The Cloud C image produces 81 scan windows, 48 accepted minima, and 29 unique
foreground coordinates. Its prepared SMF background shares the image grid.

Inside the Simon catalog ellipse, BT12 finds 53 pixels near the minimum and a
foreground of 29.3361 MJy sr\ :sup:`-1`. The conservative foreground has a
minimum, median, and maximum of 29.3361, 30.7197, and 33.9163 MJy
sr\ :sup:`-1`. It increases the count within :math:`2\sigma` from 22 to 68,
creates no strict lower limits, and raises the summed nonnegative surface
density by 3.91%. All 761,116 finite ellipse pixels satisfy
:math:`\Sigma_\mathrm{GTL}\geq\Sigma_\mathrm{BT12}`.

Held out sample tests remain difficult. Depending on the interpolation
method, the mean absolute foreground error ranges from about 5.4 to 10.6 MJy
sr\ :sup:`-1`. A numerically stable interpolation is therefore not sufficient
evidence that the inferred foreground is physically correct.

Sgr C F480M
-----------

The F480M benchmark loads its SCI and ERR extensions together. In the 123 by
123 pixel test region, the median ERR is 0.1023237 MJy sr\ :sup:`-1`. The BT12
selection contains 129 pixels, including 49 that meet the 0.74 arcsec
independence criterion, and gives a foreground of 2.3053355 MJy
sr\ :sup:`-1`.

Direct kriging fails the morphology check because it makes 3,089 pixels
strict lower limits and copies the filament into the foreground. The
conservative model retains the filament, raises the count within
:math:`2\sigma` from 55 to 117, and creates no strict lower limits. With the
same background, aperture, opacity, and mask policy, its integrated value is
31.24 M\ :sub:`sun` compared with 30.35 M\ :sub:`sun` for BT12.

The :doc:`jwst_sgrc` page defines the fixed benchmark regions and reports its
sensitivity tests.

Clouds F and H
--------------

The catalog entries ``G034.43+00.24`` and ``G035.39-00.33`` load Clouds F and
H without manually entering ellipse parameters. ``CloudF.fits`` and
``IRDCCloudH.fits`` cover their complete catalog ellipses.

.. list-table:: Change in summed nonnegative surface density relative to BT12
   :header-rows: 1

   * - Profile
     - Cloud F
     - Cloud H
   * - Conservative
     - +1.38%
     - +2.26%
   * - Moderate
     - +9.89%
     - +8.38%
   * - Liberal
     - +16.61%
     - +25.56%

Moderate GTL produces 9 strict lower limits in Cloud F and 26 in Cloud H.
Liberal GTL produces 81 and 259. Every finite input pixel inside both ellipses
still receives a finite surface density value. The ``SATURATED`` mask records
which values are censored.

See :doc:`cloud_comparisons` for the input images, sample sites, foreground
changes, and surface density maps for all three clouds.

Scientific boundary
-------------------

The current evidence checks file handling, units, masks, reproducibility, and
the expected ordering relative to BT12. It does not show that a spatial
foreground is an unbiased physical estimate. That assessment requires more
clouds, an independent tracer of column density, a source masking policy, and
a sensitivity analysis for the foreground anchor and censoring limits.
