Sgr C F480M case study
======================

Dataset and region
------------------

The validation file contains a 5660 by 2300 pixel ``SCI`` image and a
matching ``ERR`` image. Both extensions use MJy sr\ :sup:`-1`. The eastern
filament test region uses rows 1980:2103 and columns 77:200 with NumPy's
stop-exclusive slicing. The resulting cutout is 123 by 123 pixels at 0.06293
arcsec per pixel.

The final manuscript background regions were available only as a figure.
GTLMapping therefore uses five adjacent, target-sized boxes for this case
study. The calculation records their pixel coordinates so another analyst can
replace them with the publication regions.

BT12 reference
--------------

The cutout gives:

* a median ERR of 0.1023237 MJy sr\ :sup:`-1`;
* a minimum intensity of 2.3683054 MJy sr\ :sup:`-1`;
* 129 pixels between the minimum and :math:`I_\mathrm{min}+2\sigma`;
* 49 qualifying pixels at least 0.74 arcsec from the minimum; and
* a BT12 mean-minus-:math:`2\sigma` foreground of 2.3053355 MJy
  sr\ :sup:`-1`.

These values reproduce the error-based BT12 foreground calculation. The
mass comparison below uses the same background, opacity, distance, and
bright-pixel policy for BT12 and GTL.

Shared background
-----------------

The five comparison boxes touch the north, northeast, east, southeast, and
south sides of the target. Pixels above 15 MJy sr\ :sup:`-1` are excluded.
The box medians are 6.858, 7.823, 7.285, 8.684, and 5.683 MJy
sr\ :sup:`-1`; their mean is 7.267 MJy sr\ :sup:`-1`.

The spread among the boxes is part of the background uncertainty. These boxes
make the package test reproducible, but they do not recover unpublished region
coordinates.

BT12 and conservative GTL
-------------------------

An 8 by 8 overlapping GTL scan with edge completion finds 56 accepted windows
and 49 unique coordinates. The accepted minima range from 2.399 to 6.064 MJy
sr\ :sup:`-1`. That range shows why the mapper cannot treat every local
minimum as a direct foreground measurement.

The conservative fit uses a robust plane, keeps BT12 as a pointwise floor,
and reduces its spatial amplitude until the saturation limits are met. At
:math:`\kappa_\mathrm{F480M}=9.76` cm\ :sup:`2` g\ :sup:`-1`, gas-to-dust
ratio 156, and a distance of 8.15 kpc, the results are:

.. list-table::
   :header-rows: 1

   * - Quantity in the 123 by 123 pixel region
     - BT12
     - Conservative GTL
   * - Foreground min / median / max (MJy sr\ :sup:`-1`)
     - 2.305 / 2.305 / 2.305
     - 2.330 / 2.371 / 2.411
   * - Detected-pixel mass (M\ :sub:`sun`)
     - 30.35
     - 31.24
   * - Mass including finite lower limits (M\ :sub:`sun`)
     - 30.35
     - 31.24
   * - Mass ratio to BT12
     - 1.000
     - 1.030

All 15,129 pixels satisfy
:math:`\Sigma_\mathrm{GTL}\geq\Sigma_\mathrm{BT12}`. The integrated mass
therefore increases under the shared mask. The foreground fluctuation remains
broad enough that it does not divide the filament out of the extinction map.

.. image:: _static/sgrc_f480m_comparison.png
   :alt: Controlled Sgr C F480M comparison of BT12 and GTL maps
   :width: 100%

Sensitivity checks
------------------

Changing the background by one standard deviation of the box medians gives
GTL/BT12 mass ratios from 1.026 to 1.036. Values of ``grid_n`` from 4 through
10 give GTL masses from 31.16 to 31.49 M\ :sub:`sun` and ratios from 1.027 to
1.038. Dilating the bright-source mask by three and six pixels gives ratios of
1.031 and 1.033.

These tests support the sign and scale of the conservative correction in this
cutout. A publication mass still requires the final aperture, source mask,
and background regions.

Saturation counts
-----------------

The word *saturated* refers to two related masks in this analysis. BT12 first
selects pixels near the global minimum to estimate the foreground. The final
map then identifies pixels whose observed intensity falls below the fitted
foreground.

The BT12 selection contains 129 near-minimum pixels, including 49 that meet
the independence criterion. The GTL scan also has 49 unique minima, so it does
not add independent foreground measurements. Relative to the fitted
foreground, BT12 has 55 pixels within :math:`2\sigma`; conservative GTL has
117. Neither map has a strict lower-limit pixel in this cutout.

A direct kriging fit failed this check. It produced 4,387 near-saturation
pixels, 3,089 strict lower limits, and a filament-shaped foreground trough.
The current conservative model avoids that failure with a broad trend and a
blend factor of 0.03.

Uncertainty
-----------

The FITS products include first-order uncertainty from the ``ERR`` image,
background-box standard error, foreground-plane covariance, and a 30% opacity
term. Several terms are correlated across pixels, so an integrated mass
cannot use the quadrature sum of every pixel error.

Report the detected-pixel mass separately from a censored lower-limit sum.
Give the background and source-mask sensitivity ranges, then quote the opacity
scale as a correlated systematic. In this example the 30% opacity terms are
9.10 M\ :sub:`sun` for BT12 and 9.37 M\ :sub:`sun` for GTL. Saturated lower
limits do not have symmetric Gaussian errors.

A publication posterior should vary the shared background level, foreground
surface, source mask, opacity, and distance. A censored Bayesian model or a
Monte Carlo calculation can preserve the one-sided information at saturated
pixels.

Opacity normalization
---------------------

The OH94 moderately coagulated thin-ice F480M opacity is 9.76
cm\ :sup:`2` g\ :sup:`-1` at gas-to-dust ratio 156. At ratio 100, the same
filter convolution gives:

.. math::

   9.76 \times \frac{156}{100} = 15.2256\
   \;\mathrm{cm^2\,g^{-1}}.

Changing the gas-to-dust ratio changes the normalization of one opacity model.
It does not create a second F480M dust law.

Run the case study
------------------

Install the plotting extra, then run from the repository root:

.. code-block:: console

   python examples/sgrc_f480m_compare.py F480M_registered.fits

The script writes BT12 and GTL FITS maps, a comparison figure, and a JSON file
under ``validation/sgrc_f480m``. Replace ``touching_background_boxes`` when
the final machine-readable background regions are available.
