Scientific method
=================

Radiative transfer
------------------

With observed intensity :math:`I_\mathrm{obs}`, foreground intensity
:math:`I_\mathrm{fg}`, and the observed off-cloud background model
:math:`I_\mathrm{bg}`, GTLMapping evaluates:

.. math::

   \tau = -\ln\left(
       \frac{I_\mathrm{obs} - I_\mathrm{fg}}
            {I_\mathrm{bg} - I_\mathrm{fg}}
   \right),

and:

.. math::

   \Sigma = \frac{\tau}{\kappa_{\mathrm{filter}}}.

When neither a filter nor an explicit opacity is supplied, the
compatibility default :math:`\kappa_{8\mu\mathrm{m}}` is 7.5
cm\ :sup:`2` g\ :sup:`-1`, following the BT09/BT12 fiducial model. It
is a parameter, not a universal constant.

Foreground samples
------------------

The GTL method scans the full image with overlapping windows. The
minimum in a window is accepted when another pixel is:

1. within two instrumental-noise standard deviations of the minimum;
2. inside an optional user-supplied region mask; and
3. at least the configured independent separation away.

The GTL prototype defines the default as four beam radii: 4 arcsec for
a 2-arcsec FWHM beam. BT12 used an 8-arcsec independent-core
criterion. ``min_separation_arcsec`` is explicit so this scientific
choice can be tested.

Overlapping windows may find the same coordinate repeatedly.
GTLMapping merges those duplicates for storage and records their
multiplicity. Weighted-kernel methods can use that multiplicity.
Stable kriging fits the unique coordinates with a pseudo-inverse;
``kriging_duplicate_policy="repeat"`` restores the repeated rows and
is numerically equivalent to the prototype's raw list.

The default scan origins reproduce the prototype exactly. This means a
few edge pixels may not fall inside a complete window.
``cover_edges=True`` appends a final valid origin. The search is not
restricted to the Simon ellipse unless ``restrict_to_cloud=True`` is
requested.

BT12 constant foreground
------------------------

BT12 is a separate scientific method, not a flat interpolation of GTL
samples. Inside the selected IRDC, GTLMapping:

1. finds the global observed-intensity minimum;
2. labels all pixels strictly between the minimum and
   :math:`I_\mathrm{min}+2\sigma` as saturated if at least one is
   8 arcsec or more from the minimum;
3. evaluates the mean of all labeled saturated pixels; and
4. sets :math:`I_\mathrm{fg}` to that mean minus :math:`2\sigma`.

Use ``fit_foreground(method="bt12")`` for this comparison. The
prototype notebook's BT12 cell instead returned a median and then
added 1.0 MJy sr\ :sup:`-1`; that is not the BT12 prescription and is
not reproduced as a package method.

Conservative spatial foreground
-------------------------------

A local minimum is not automatically a saturated foreground
measurement. In a structured cloud, exact interpolation through all
window minima can reproduce the cloud morphology in
:math:`I_\mathrm{fg}` and divide it out of the extinction map.

The default ``conservative`` model therefore separates the identifiable
parts of the problem and follows the BT12-floor rule recovered from the
later prototype history:

1. BT12 fixes the absolute foreground level.
2. The local minima constrain only a robust planar spatial candidate.
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

This construction makes the intended comparison mathematically
testable. For a fixed observed image and background with
:math:`I_\mathrm{bg}>I_\mathrm{obs}>I_\mathrm{fg}`, optical depth is
monotonic in foreground:

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{fg}} =
   \frac{1}{I_\mathrm{obs}-I_\mathrm{fg}} -
   \frac{1}{I_\mathrm{bg}-I_\mathrm{fg}} > 0.

Because the conservative GTL foreground is never below BT12, its
surface density cannot be below BT12 on a jointly valid pixel. A
detected-only sum can otherwise appear to decrease if a more
aggressive foreground creates newly masked saturated pixels; the
zero-new-strict-saturation default prevents that bookkeeping failure
inside the fitted region.

Moderate spatial foreground
---------------------------

The ``moderate`` profile uses the controlled quadratic GTL trend with a 50%
soft BT12 anchor. Its default near-saturation budget is 0.5%, and its strict
censoring ceiling is 0.01%. It therefore permits broad spatial foreground
structure and more locally saturation-consistent pixels than the BT12
comparison while producing far fewer lower-limit markers than liberal GTL.

``compute_moderate`` applies the same explicit transmitted-intensity floor and
foreground/background feasibility projection used for liberal products.
Strict pixels remain finite lower limits marked in ``SATURATED``; feasibility
adjustments remain marked in ``FG_CONSTRAINT``. The moderate profile is a
documented parameter preset, not an independent physical law, and its three
defaults should be included in sensitivity tests.

Liberal spatial foreground
--------------------------

The ``liberal`` option asks a different scientific question: how much more
extinction is supported if the GTL local minima set the broad foreground level
without a hard BT12 floor? It fits a robust quadratic surface to the accepted
GTL samples, clips extrapolation to their observed intensity range, subtracts
the foreground margin, and shifts the resulting surface subject to two
explicit budgets:

* ``target_local_saturation_fraction`` limits pixels satisfying
  :math:`I_\mathrm{obs}\leq I_\mathrm{fg}+2\sigma`; and
* ``maximum_strict_saturation_fraction`` separately limits censored pixels
  satisfying :math:`I_\mathrm{obs}\leq I_\mathrm{fg}`.

The defaults are 1% and 0.1%, respectively. BT12 is evaluated for comparison
diagnostics but has zero weight in the liberal foreground by default. A nonzero
``bt12_anchor_weight`` creates a documented soft anchor; it never becomes the
conservative method's hard floor.

The saturation budgets are modeling assumptions, not new empirical laws.
Publication analyses should report sensitivity to both fractions and should
inspect the foreground, residual, ``SATURATED``, and ``FG_CONSTRAINT`` maps.
The liberal mode intentionally avoids raw interpolation through every minimum,
which previously reproduced cloud structure in the foreground.

Strictly saturated pixels do not have measured optical depths. For them,
``compute_liberal`` substitutes a positive transmitted-intensity floor
(default :math:`2\sigma`) and returns a finite lower limit while preserving the
``SATURATED`` flag. Before the radiative-transfer calculation it also enforces

.. math::

   I_\mathrm{fg} \leq I_\mathrm{bg} - I_\mathrm{trans,min},

recording every adjustment in ``FG_CONSTRAINT``. Thus liberal products avoid
model-generated NaNs without pretending that a censored lower limit is an
ordinary measurement. Non-finite input data remain non-finite, and a
background below the requested transmitted-intensity floor raises an error
because no physical projection exists.

Legacy foreground interpolation
-------------------------------

Kriging, RBF, spline, Gaussian, Cauchy, and a constant sample summary
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

This is intentionally different from the notebook's sparse raw-pixel
stencil, which did not implement the published BT09 procedure. The
implementation follows the published algorithmic description, but it
is not claimed to be bit-for-bit identical to the authors' original
unpublished code.

``estimate_box_background`` implements the JWST Sgr C prescription:
measure the median in each adjacent region after rejecting intensities
above 15 MJy sr\ :sup:`-1`, then take the mean of those medians. The
threshold and pixel boxes are explicit parameters.

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

and writes an ``FG_CONSTRAINT`` mask. This is a transparent projection
onto the feasible radiative-transfer range, not evidence that the
original foreground/background models were correct.

When :math:`I_\mathrm{obs} > I_\mathrm{bg}`, the optical depth is
negative. ``bright_pixel_policy="allow"`` keeps all such values,
preserving BT09's bias-avoidance motivation. ``"zero"`` and ``"mask"``
are explicit alternatives. These policies do not reproduce BT09's
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
foreground, and opacity terms. ``GTLMapper`` automatically uses an
attached ERR image and an available foreground-model variance. The
total arrays are written as ``TAU_ERR`` and ``SIGMA_ERR``.

The approximation is masked at saturated pixels and becomes unreliable
near :math:`N=0`. For publication inference, treat saturation as
censoring and use Monte Carlo or a hierarchical model that can include
foreground/background covariance. Separate opacity, gas/dust ratio,
and background systematics from random pixel noise.

Scientific scope and remaining limits
-------------------------------------

The package currently implements a parameterized single-band MIREX
radiative-transfer calculation. It provides the supplied
filter-convolved opacity table but does not yet convolve arbitrary
throughput curves with raw dust-model tables, infer a spatially varying
extinction law, merge NIR and MIR maps, model covariance, or perform
censored Monte Carlo inference. The compatibility default
:math:`\kappa_{8\mu\mathrm{m}}=7.5` cm\ :sup:`2` g\ :sup:`-1` has an
estimated absolute systematic uncertainty of about 30 percent in the
BT09/BT12/KT13 lineage.

For filters outside the registry, users must supply an appropriate
opacity and noise estimate and verify that cold-cloud emission is
negligible. Independent submillimeter calibration and automated
point-source masks remain future work.
