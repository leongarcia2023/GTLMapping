Scientific method
=================

Radiative transfer
------------------

GTLMapping combines the observed intensity :math:`I_\mathrm{obs}`, foreground
intensity :math:`I_\mathrm{fg}`, and off-cloud background model
:math:`I_\mathrm{bg}` as follows:

.. math::

   \tau = -\ln\left(
       \frac{I_\mathrm{obs} - I_\mathrm{fg}}
            {I_\mathrm{bg} - I_\mathrm{fg}}
   \right),

and:

.. math::

   \Sigma = \frac{\tau}{\kappa_{\mathrm{filter}}}.

If you supply neither a filter nor an opacity, GTLMapping uses the
compatibility value :math:`\kappa_{8\mu\mathrm{m}}=7.5`
cm\ :sup:`2` g\ :sup:`-1`, following the BT09/BT12 fiducial model. It
is a configurable reference value.

Foreground samples
------------------

The GTL method scans the full image with overlapping windows. The
minimum in a window is accepted when another pixel is:

1. within two instrumental-noise standard deviations of the minimum;
2. inside an optional user-supplied region mask; and
3. at least the configured independent separation away.

The GTL prototype uses four beam radii, or 4 arcsec for a 2-arcsec FWHM
beam. BT12 used an 8-arcsec independent-core criterion. Set
``min_separation_arcsec`` to match the dataset and include the choice in
sensitivity tests.

Overlapping windows may find the same coordinate repeatedly.
GTLMapping merges those duplicates for storage and records their
multiplicity. Weighted-kernel methods can use that multiplicity.
Stable kriging fits the unique coordinates with a pseudo-inverse;
``kriging_duplicate_policy="repeat"`` restores the repeated rows and
is numerically equivalent to the prototype's raw list.

The default scan origins reproduce the prototype, leaving a few edge pixels
outside a complete window.
``cover_edges=True`` appends a final valid origin. The search is not
restricted to the Simon ellipse unless ``restrict_to_cloud=True`` is
requested.

BT12 constant foreground
------------------------

GTLMapping implements BT12 independently of the ``flat`` summary of GTL
samples. Inside the selected IRDC, the BT12 routine:

1. finds the global observed-intensity minimum;
2. labels all pixels strictly between the minimum and
   :math:`I_\mathrm{min}+2\sigma` as saturated if at least one is
   8 arcsec or more from the minimum;
3. evaluates the mean of all labeled saturated pixels; and
4. sets :math:`I_\mathrm{fg}` to that mean minus :math:`2\sigma`.

Use ``fit_foreground(method="bt12")`` for the published comparison. The
prototype notebook used a median and added 1.0 MJy sr\ :sup:`-1`, so it does
not reproduce the BT12 prescription.

Conservative spatial foreground
-------------------------------

A structured cloud can place a local minimum above the true foreground.
Interpolating through every window minimum then copies cloud morphology into
:math:`I_\mathrm{fg}` and removes that structure from the extinction map.

The ``conservative`` model assigns BT12 and the local minima different roles:

1. BT12 fixes the absolute foreground level.
2. The local minima constrain a robust planar spatial candidate.
3. The same foreground margin is subtracted from that candidate.
4. The final model is the pointwise maximum of the spatial candidate
   and the BT12 foreground.
5. A blend factor is reduced until no more than the configured
   fractions of pixels are locally saturation-consistent or strictly
   below the foreground.

The default limits are 1% for
:math:`I_\mathrm{obs}\leq I_\mathrm{fg}+2\sigma` and zero new pixels
for :math:`I_\mathrm{obs}\leq I_\mathrm{fg}`. The reference BT12
counts are always allowed, so the regularizer cannot make an already
saturated dataset impossible to fit. The resulting plane covariance
is retained as a foreground-uncertainty estimate.

For a fixed observed image and background with
:math:`I_\mathrm{bg}>I_\mathrm{obs}>I_\mathrm{fg}`, optical depth is
monotonic in foreground:

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{fg}} =
   \frac{1}{I_\mathrm{obs}-I_\mathrm{fg}} -
   \frac{1}{I_\mathrm{bg}-I_\mathrm{fg}} > 0.

The pointwise BT12 floor makes the conservative surface density greater than
or equal to BT12 on every jointly valid pixel. A stronger foreground can still
reduce a detected-only sum by creating masked saturated pixels. The default
strict-saturation limit prevents that loss inside the fitted region.

Moderate spatial foreground
---------------------------

The ``moderate`` profile fits a quadratic GTL trend with a 50% soft BT12
anchor. It allows 0.5% of pixels near saturation and caps strict censoring at
0.01%. These settings admit more foreground variation than the conservative
fit while producing fewer lower-limit pixels than liberal GTL.

``compute_moderate`` applies the transmitted-intensity floor and
foreground/background feasibility projection used for liberal products.
Strict pixels remain finite lower limits marked in ``SATURATED``; feasibility
adjustments remain marked in ``FG_CONSTRAINT``. The moderate profile is a
parameter preset. Include its anchor and censoring limits in sensitivity
tests.

Liberal spatial foreground
--------------------------

The ``liberal`` option lets the GTL minima set the broad foreground level. It
fits a robust quadratic surface to the accepted samples, clips extrapolation
to their intensity range, subtracts the foreground margin, and applies two
censoring budgets:

* ``target_local_saturation_fraction`` limits pixels satisfying
  :math:`I_\mathrm{obs}\leq I_\mathrm{fg}+2\sigma`; and
* ``maximum_strict_saturation_fraction`` separately limits censored pixels
  satisfying :math:`I_\mathrm{obs}\leq I_\mathrm{fg}`.

The defaults are 1% and 0.1%, respectively. BT12 is evaluated for comparison
diagnostics but has zero weight in the liberal foreground by default. A nonzero
``bt12_anchor_weight`` creates a documented soft anchor; it never becomes the
conservative method's hard floor.

The saturation budgets are user-set assumptions. Publication analyses should
report sensitivity to both fractions and
inspect the foreground, residual, ``SATURATED``, and ``FG_CONSTRAINT`` maps.
The quadratic trend avoids the cloud-shaped foreground produced by raw
interpolation through every minimum in the Sgr C test.

Strictly saturated pixels do not have measured optical depths. For them,
``compute_liberal`` substitutes a positive transmitted-intensity floor
(default :math:`2\sigma`) and returns a finite lower limit while preserving the
``SATURATED`` flag. Before the radiative-transfer calculation it also enforces

.. math::

   I_\mathrm{fg} \leq I_\mathrm{bg} - I_\mathrm{trans,min},

recording every adjustment in ``FG_CONSTRAINT``. The result contains finite
lower limits instead of model-generated NaNs. Non-finite input data remain
non-finite, and a
background below the requested transmitted-intensity floor raises an error
because no physical projection exists.

Legacy foreground interpolation
-------------------------------

Kriging, RBF, spline, Gaussian, Cauchy, and the constant sample summary
called ``flat`` remain available for experiments. Kriging uses unique
coordinates and a pseudo-inverse; numerical gaps are nearest-sample
filled and recorded. These methods interpolate the minima directly and
are no longer the high-level default.

The prototype used a 1.0-MJy sr\ :sup:`-1` margin. Set
``method="kriging"``, ``foreground_margin=1.0``,
``clip_to_sample_range=False``, and
``kriging_duplicate_policy="repeat"`` only when reproducing that
historical calculation. A finite or numerically stable interpolation
does not establish that the foreground is physically valid.

Background models
-----------------

The LMF implementation approximates the published procedure by
sampling a large trimmed median on a
24-arcsecond grid. The SMF implementation:

1. uses a local filter one third of the Simon major axis outside the
   cloud ellipse;
2. excludes the cloud from those local samples; and
3. fills the ellipse using inverse-square weighting within one
   semi-major-axis radius.

The notebook used a sparse raw-pixel stencil that did not implement the
published BT09 procedure. GTLMapping follows the paper's algorithmic
description; the original unpublished code was unavailable for a bitwise
comparison.

``estimate_box_background`` implements the JWST Sgr C prescription:
measure the median in each adjacent region after rejecting intensities
above 15 MJy sr\ :sup:`-1`, then take the mean of those medians. Users supply
the threshold and pixel boxes.

Filter-aware opacity
--------------------

``get_filter_opacity`` contains the filter-convolved table from the
supplied Sgr C analysis for IRAC2, F480M, IRAC4, F770W, and F2100W
under five OH94/WD01 dust models. The reference normalization is a
gas-to-dust ratio of 156. Since the value is per total mass:

.. math::

   \kappa(R) = \kappa(156)\frac{156}{R}.

Thus the moderately coagulated OH94 thin-ice F480M opacity is 9.76
cm\ :sup:`2` g\ :sup:`-1` at ratio 156 and 15.2256
cm\ :sup:`2` g\ :sup:`-1` at ratio 100. The manuscript's two values
are consistent once this normalization is stated.

Masks and negative values
-------------------------

Pixels with :math:`I_\mathrm{obs} \le I_\mathrm{fg}` or
:math:`I_\mathrm{bg} \le I_\mathrm{fg}` are masked by default. They
are not silently converted to zero. ``saturation_policy="lower_limit"``
uses an explicit positive transmitted-intensity floor for the first
case and preserves the ``SATURATED`` mask.

Kriging can be numerically finite while disagreeing with a separately
estimated background. ``constrain_foreground`` enforces

.. math::

   I_\mathrm{fg} \le I_\mathrm{bg} - I_\mathrm{trans,min}

and writes an ``FG_CONSTRAINT`` mask. The projection enforces the feasible
radiative-transfer range. The mask identifies pixels where the original
foreground and background disagreed.

When :math:`I_\mathrm{obs} > I_\mathrm{bg}`, the optical depth is
negative. ``bright_pixel_policy="allow"`` keeps those values and preserves
BT09's bias-avoidance treatment. ``"zero"`` and ``"mask"`` provide other
policies. None reproduces BT09's
separate bright-source cleaning rule, which set only pixels above a
half-FWHM intensity threshold to zero. Significant stellar and
extended emission should therefore be masked upstream.

One g cm\ :sup:`-2` is approximately 4,788.45 solar masses
pc\ :sup:`-2`.

Uncertainty propagation
-----------------------

For :math:`N=I_\mathrm{obs}-I_\mathrm{fg}` and
:math:`D=I_\mathrm{bg}-I_\mathrm{fg}`, the first-order derivatives are:

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{obs}}=-\frac{1}{N},\quad
   \frac{\partial\tau}{\partial I_\mathrm{bg}}=\frac{1}{D},\quad
   \frac{\partial\tau}{\partial I_\mathrm{fg}}=
   \frac{1}{N}-\frac{1}{D}.

``propagate_uncertainty`` combines independent observed, background,
foreground, and opacity terms. ``GTLMapper`` uses an
attached ERR image and an available foreground-model variance. The
total arrays are written as ``TAU_ERR`` and ``SIGMA_ERR``.

The approximation is masked at saturated pixels and becomes unreliable
near :math:`N=0`. For publication inference, treat saturation as
censoring and use Monte Carlo or a hierarchical model that can include
foreground/background covariance. Separate opacity, gas/dust ratio,
and background systematics from random pixel noise.

Scientific scope and remaining limits
-------------------------------------

The package implements a parameterized single-band MIREX calculation and the
supplied filter-convolved opacity table. It does not yet convolve arbitrary
throughput curves with raw dust tables or infer a spatially varying extinction
law. NIR/MIR merging, covariance models, and censored Monte Carlo inference
also remain outside the current release. The compatibility default
:math:`\kappa_{8\mu\mathrm{m}}=7.5` cm\ :sup:`2` g\ :sup:`-1` has an
estimated absolute systematic uncertainty of about 30 percent in the
BT09/BT12/KT13 lineage.

For filters outside the registry, users must supply an appropriate
opacity and noise estimate and verify that cold-cloud emission is
negligible. Independent submillimeter calibration and automated
point-source masks remain future work.
