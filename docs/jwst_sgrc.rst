Sgr C F480M case study
========================

Data and geometry
-------------------

The example reads SCI and ERR extensions of a registered 5660 by 2300 pixel
F480M image in MJy/sr. The target uses columns 77:200 and rows 1980:2103
with exclusive upper bounds: 123 by 123 pixels at 0.06293 arcsec per pixel.

Five touching boxes define a fixed background. These are array coordinates;
they do not assume that increasing row means north.

.. list-table:: Background boxes
   :header-rows: 1

   * - Rows
     - Columns
     - Median (MJy/sr)
   * - 1857:1980
     - 77:200
     - 6.858
   * - 1857:1980
     - 200:323
     - 7.823
   * - 1980:2103
     - 200:323
     - 7.285
   * - 2103:2226
     - 200:323
     - 8.684
   * - 2103:2226
     - 77:200
     - 5.683

Pixels above 15 MJy/sr are excluded from box measurements. The mean of
medians is 7.26665 MJy/sr and their sample standard deviation is 1.11735.
These boxes define a benchmark; they have not been established as an unbiased
background for a larger Sgr C region.

Foregrounds and limits
------------------------

The adopted noise is the median target ERR, 0.1023237 MJy/sr. BT12 selects
129 pixels near the minimum, including 49 beyond the 0.74 arcsec partner
separation, and gives a foreground of 2.3053355 MJy/sr.

An 8 by 8 GTL scan with edge completion yields 49 unique sites.
Conservative GTL selects a planar blend of 0.03, with foreground ranging
from 2.33023 to 2.41089 MJy/sr.

Both profiles use a threshold of 0.2046474 MJy/sr. Positive residuals below
that threshold are unresolved. Neither map has a strict zero crossing.

.. list-table:: Gas mass at 8.15 kpc and adopted F480M opacity 9.76 cm²/g
   :header-rows: 1

   * - Quantity
     - BT12
     - Conservative
   * - Detected-only mass (solar masses)
     - 29.768
     - 29.973
   * - Detected mass plus sensitivity limits (solar masses)
     - 30.300
     - 31.101
   * - Unresolved pixels
     - 55
     - 117
   * - Strict zero-crossing pixels
     - 0
     - 0

The finite-map ratio is 1.02644. Detected-only totals use different resolved
subsets and are not a fixed-pixel recovery comparison. The true column
beyond each limit is not measured.

.. image:: _static/sgrc_f480m_comparison.png
   :alt: Sgr C intensity, BT12 and conservative surface densities, their difference, grid and background sensitivity
   :width: 100%

Sensitivity and masks
-----------------------

Scan grids from 4 through 10 give ratios of 1.0242 to 1.0333. Varying the
background over the box scatter gives BT12 masses of 21.29 to 38.40 solar
masses, GTL masses of 21.96 to 39.29, and ratios of 1.0234 to 1.0314.

The common bright-pixel rule identifies intensity above the background.
It excludes 2507 pixels (16.57%) from positive-mass inference, but is not a
complete point-source classifier. Three- and six-pixel mask dilations give
ratios of 1.0274 and 1.0287.

The image-noise threshold omits foreground uncertainty and resampling
covariance. Analytic FITS uncertainties are masked on unresolved pixels.
Their quadrature sum is not a cloud-mass error. The nominal 30% opacity
term is a shared scale uncertainty; background and mask changes are
sensitivity tests.

Opacity and reproduction
--------------------------

The benchmark adopts a gas-mass convention for the supplied F480M value:
9.76 cm²/g at gas/dust=156. At gas/dust=100, the same adopted dust opacity
gives 15.2256 per gram of gas. Exact total-mass normalization uses R+1 in
the denominator. See :doc:`scientific_method`.

The rounded value is a benchmark input, not an independent filter
convolution. A calibrated application should archive the dust and
throughput curves.

.. code-block:: console

   python examples/sgrc_f480m_compare.py F480M_registered.fits

The script writes FITS maps and a JSON record of regions, masks, thresholds
and sensitivity tests. Choose and document new regions for a new dataset.
