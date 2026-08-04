Sgr C F480M case study
======================

Dataset and region
------------------

This benchmark uses a FITS file with a 5660 by 2300 pixel ``SCI`` image and a
matching ``ERR`` image. Both extensions use MJy sr\ :sup:`-1`. The eastern
filament region uses rows 1980:2103 and columns 77:200 with NumPy's exclusive
upper bounds. The resulting cutout is 123 by 123 pixels at 0.06293 arcsec per
pixel.

Five adjacent boxes with the same dimensions as the target define a fixed
background measurement for this example. Their coordinates are recorded by
the script so the calculation can be repeated exactly. Analyses of other
regions should choose background boxes from the relevant image and document
that choice.

BT12 reference
--------------

The cutout gives:

* a median ERR of 0.1023237 MJy sr\ :sup:`-1`;
* a minimum intensity of 2.3683054 MJy sr\ :sup:`-1`;
* 129 pixels between the minimum and :math:`I_\mathrm{min}+2\sigma`;
* 49 qualifying pixels at least 0.74 arcsec from the minimum; and
* a BT12 foreground of 2.3053355 MJy sr\ :sup:`-1`, calculated as the mean
  minus :math:`2\sigma`.

The comparison uses the same background, opacity, distance, aperture, and
bright pixel policy for BT12 and GTL.

Shared background
-----------------

The five comparison boxes touch the north, northeast, east, southeast, and
south sides of the target. Pixels above 15 MJy sr\ :sup:`-1` are excluded.
The box medians are 6.858, 7.823, 7.285, 8.684, and 5.683 MJy
sr\ :sup:`-1`; their mean is 7.267 MJy sr\ :sup:`-1`.

The spread among these medians contributes to the background uncertainty.
The boxes make this a reproducible package benchmark, not a universal
background prescription for Sgr C.

BT12 and conservative GTL
-------------------------

An 8 by 8 overlapping GTL scan with edge completion finds 56 accepted windows
and 49 unique coordinates. The accepted minima range from 2.399 to 6.064 MJy
sr\ :sup:`-1`. This range is one reason not to treat every local minimum as a
direct measurement of foreground emission.

The conservative fit uses a robust plane, keeps BT12 as a pointwise floor,
and reduces its spatial amplitude until the saturation limits are met. At
:math:`\kappa_\mathrm{F480M}=9.76` cm\ :sup:`2` g\ :sup:`-1`, a gas to dust
ratio of 156, and a distance of 8.15 kpc, the results are:

.. list-table::
   :header-rows: 1

   * - Quantity in the 123 by 123 pixel region
     - BT12
     - Conservative GTL
   * - Foreground min / median / max (MJy sr\ :sup:`-1`)
     - 2.305 / 2.305 / 2.305
     - 2.330 / 2.371 / 2.411
   * - Integrated surface density (M\ :sub:`sun`)
     - 30.35
     - 31.24
   * - Ratio to BT12
     - 1.000
     - 1.030

All 15,129 pixels satisfy
:math:`\Sigma_\mathrm{GTL}\geq\Sigma_\mathrm{BT12}`. The integrated value
therefore increases under the shared mask. The foreground fluctuation remains
broad enough that it does not divide the filament out of the extinction map.

.. image:: _static/sgrc_f480m_comparison.png
   :alt: Sgr C F480M comparison showing the observed image, background boxes, BT12 foreground, conservative foreground, and both surface density maps
   :width: 100%

Sensitivity checks
------------------

Changing the background by one standard deviation of the box medians gives
GTL to BT12 ratios from 1.026 to 1.036. Values of ``grid_n`` from 4 through 10
give GTL values from 31.16 to 31.49 M\ :sub:`sun` and ratios from 1.027 to
1.038. Dilating the mask around bright sources by three and six pixels gives
ratios of 1.031 and 1.033.

These checks support the sign and scale of the conservative correction in
this cutout. Results for a larger scientific region still depend on the
chosen aperture, source mask, and background regions.

Saturation counts
-----------------

The word *saturated* refers to two related masks. BT12 first selects pixels
near the global minimum to estimate the foreground. The final map then
identifies pixels whose observed intensity falls below the fitted foreground.

The BT12 selection contains 129 pixels near the minimum, including 49 that
meet the independence criterion. The GTL scan also has 49 unique minima, so
it does not add independent measurements in this cutout. Relative to the
fitted foreground, BT12 has 55 pixels within :math:`2\sigma`; conservative
GTL has 117. Neither map has a strict lower limit.

A direct kriging fit fails this check. It produces 4,387 pixels near
saturation, 3,089 strict lower limits, and a foreground trough shaped like the
filament. The conservative model avoids that behavior with a broad trend and
a blend factor of 0.03.

Uncertainty
-----------

The FITS products include first order uncertainty from the ``ERR`` image, the
standard error among background boxes, the covariance of the foreground
plane, and a 30% opacity term. Several terms are correlated across pixels, so
an integrated uncertainty cannot use the quadrature sum of every pixel error.

Report detections separately from censored lower limits. Give sensitivity
ranges for the background and source mask, then quote the opacity scale as a
correlated systematic. In this example the 30% opacity terms are 9.10
M\ :sub:`sun` for BT12 and 9.37 M\ :sub:`sun` for GTL. Censored lower limits
do not have symmetric Gaussian errors.

A fuller inference can vary the shared background level, foreground surface,
source mask, opacity, and distance. A censored Bayesian model or Monte Carlo
calculation can preserve the information carried by lower limits.

Opacity normalization
---------------------

The F480M opacity for moderately coagulated OH94 grains with thin ice mantles
is 9.76 cm\ :sup:`2` g\ :sup:`-1` at gas to dust ratio 156. At ratio 100, the
same filter convolution gives:

.. math::

   9.76 \times \frac{156}{100} = 15.2256\
   \;\mathrm{cm^2\,g^{-1}}.

Changing the gas to dust ratio changes the normalization of one opacity
model. It does not create a second F480M dust law.

Run the case study
------------------

Install the plotting extra, then run from the repository root:

.. code-block:: console

   python examples/sgrc_f480m_compare.py F480M_registered.fits

The script writes BT12 and GTL FITS maps, a comparison figure, and a JSON file
under ``validation/sgrc_f480m``. Edit the target and background regions in the
script when applying the calculation to another dataset.
