Scientific method
=================

Radiative transfer
------------------

GTLMapping combines the observed intensity :math:`I_\mathrm{obs}`, foreground
intensity :math:`I_\mathrm{fg}`, and background model outside the cloud
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

1. within two standard deviations of the instrumental noise from the minimum;
2. inside an optional region mask supplied by the user; and
3. at least the configured independent separation away.

The legacy calculation uses four beam radii, or 4 arcsec for a beam with a 2
arcsec FWHM. BT12 required independent cores to be separated by 8 arcsec. Set
``min_separation_arcsec`` to match the dataset and include the choice in
sensitivity tests.

Overlapping windows may find the same coordinate repeatedly.
GTLMapping merges those duplicates for storage and records their
multiplicity. Kernel methods can use the multiplicity as a weight.
Stable kriging fits the unique coordinates with a pseudoinverse.
``kriging_duplicate_policy="repeat"`` restores the repeated rows and is
numerically equivalent to the original calculation.

The default scan origins preserve the legacy sample grid, leaving a few edge
pixels outside a complete window.
``cover_edges=True`` appends a final valid origin. The search is not
restricted to the Simon ellipse unless ``restrict_to_cloud=True`` is
requested.

BT12 constant foreground
------------------------

GTLMapping implements BT12 independently of the ``flat`` summary of GTL
samples. Inside the selected IRDC, the BT12 routine:

1. finds the minimum observed intensity in the full image;
2. labels all pixels strictly between the minimum and
   :math:`I_\mathrm{min}+2\sigma` as saturated if at least one is
   8 arcsec or more from the minimum;
3. evaluates the mean of all labeled saturated pixels; and
4. sets :math:`I_\mathrm{fg}` to that mean minus :math:`2\sigma`.

Use ``fit_foreground(method="bt12")`` for the BT12 comparison. This estimator
is separate from the legacy ``flat`` sample summary.

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
   fractions of pixels are consistent with local saturation or fall strictly
   below the foreground.

The default limits are 1% for
:math:`I_\mathrm{obs}\leq I_\mathrm{fg}+2\sigma` and zero new pixels
for :math:`I_\mathrm{obs}\leq I_\mathrm{fg}`. The reference BT12
counts are always allowed, so the regularizer cannot make an already
saturated dataset impossible to fit. The covariance of the fitted plane
provides an estimate of foreground uncertainty.

For a fixed observed image and background with
:math:`I_\mathrm{bg}>I_\mathrm{obs}>I_\mathrm{fg}`, optical depth is
monotonic in foreground:

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{fg}} =
   \frac{1}{I_\mathrm{obs}-I_\mathrm{fg}} -
   \frac{1}{I_\mathrm{bg}-I_\mathrm{fg}} > 0.

The named profiles enforce the pointwise ordering

.. math::

   I_{\mathrm{fg,BT12}} \leq I_{\mathrm{fg,cons}}
   \leq I_{\mathrm{fg,mod}} \leq I_{\mathrm{fg,lib}}.

Their surface densities therefore have the same ordering on every jointly
valid, uncensored pixel. Moderate and liberal use progressively larger,
explicit budgets for censored lower limits; comparisons must retain those
masks rather than treating the limits as ordinary detections.

Moderate spatial foreground
---------------------------

The ``moderate`` profile fits a quadratic GTL trend with a 50% two-sided BT12
pull, then applies conservative GTL as a one-sided pointwise floor. It allows
0.5% of pixels near saturation and caps strict censoring at 0.01%. These
settings admit more foreground variation than the conservative fit while
producing fewer pixels with lower limits than liberal GTL.

``compute_moderate`` applies the floor on transmitted intensity and
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

The defaults are 1% and 0.1%, respectively. The raw liberal quadratic has zero
two-sided BT12 weight by default, but the final named profile uses moderate GTL
as a one-sided pointwise floor. Thus the samples determine every positive
spatial enhancement while the method cannot undo mass already inferred by a
less permissive profile. A nonzero ``bt12_anchor_weight`` adds a documented
two-sided pull to the raw quadratic.

The user chooses the saturation budgets. Scientific analyses should report
sensitivity to both fractions and
inspect the foreground, residual, ``SATURATED``, and ``FG_CONSTRAINT`` maps.
The quadratic trend avoids the foreground shaped like the cloud that raw
interpolation produced in the Sgr C test.

Strictly saturated pixels do not have measured optical depths. For them,
``compute_liberal`` substitutes a positive floor on transmitted intensity
(default :math:`2\sigma`) and returns a finite lower limit while preserving the
``SATURATED`` flag. Before calculating radiative transfer, it also enforces

.. math::

   I_\mathrm{fg} \leq I_\mathrm{bg} - I_\mathrm{trans,min},

recording every adjustment in ``FG_CONSTRAINT``. The result contains finite
lower limits instead of NaNs created by the model. Input values that are not
finite remain unchanged. A background below the requested floor on transmitted
intensity raises an error because no physical projection exists.

Legacy foreground interpolation
-------------------------------

Kriging, RBF, spline, Gaussian, Cauchy, and the constant sample summary
called ``flat`` remain available for experiments. Kriging uses unique
coordinates and a pseudoinverse; the nearest sample fills any numerical gaps,
which the result records. These methods interpolate the minima directly and
are no longer the default in the main interface.

The legacy calculation used a margin of 1.0 MJy sr\ :sup:`-1`. Set
``method="kriging"``, ``foreground_margin=1.0``,
``clip_to_sample_range=False``, and
``kriging_duplicate_policy="repeat"`` only when reproducing that calculation.
A finite or numerically stable interpolation
does not establish that the foreground is physically valid.

Background models
-----------------

The LMF implementation approximates the published procedure by
sampling a large trimmed median on a grid spaced by 24 arcseconds. The SMF
implementation:

1. uses a local filter one third of the Simon major axis outside the
   cloud ellipse;
2. excludes the cloud from those local samples; and
3. fills the ellipse using weights proportional to inverse squared distance
   within one semimajor axis.

Earlier experiments used a sparse stencil of raw pixels that did not implement
the published BT09 procedure. GTLMapping follows the algorithmic description
in the paper.

``estimate_box_background`` measures a background from adjacent regions:
measure the median in each adjacent region after rejecting intensities
above 15 MJy sr\ :sup:`-1`, then take the mean of those medians. Users supply
the threshold and pixel boxes.

Opacity by filter
-----------------

``get_filter_opacity`` contains filter convolved values for IRAC2, F480M,
IRAC4, F770W, and F2100W. The table covers five OH94 and WD01 dust models. The
reference normalization is a gas to dust ratio of 156. Since the value is per
total mass:

.. math::

   \kappa(R) = \kappa(156)\frac{156}{R}.

Thus the F480M opacity for moderately coagulated OH94 grains with thin ice
mantles is 9.76
cm\ :sup:`2` g\ :sup:`-1` at ratio 156 and 15.2256
cm\ :sup:`2` g\ :sup:`-1` at ratio 100. Both values describe the same dust
model under different total mass normalizations.

Masks and negative values
-------------------------

Pixels with :math:`I_\mathrm{obs} \le I_\mathrm{fg}` or
:math:`I_\mathrm{bg} \le I_\mathrm{fg}` are masked by default. They
are not silently converted to zero. ``saturation_policy="lower_limit"``
uses an explicit positive floor on transmitted intensity for the first
case and preserves the ``SATURATED`` mask.

Kriging can be numerically finite while disagreeing with a separately
estimated background. ``constrain_foreground`` enforces

.. math::

   I_\mathrm{fg} \le I_\mathrm{bg} - I_\mathrm{trans,min}

and writes an ``FG_CONSTRAINT`` mask. The projection enforces the feasible
range for radiative transfer. The mask identifies pixels where the original
foreground and background disagreed.

When :math:`I_\mathrm{obs} > I_\mathrm{bg}`, the optical depth is
negative. ``bright_pixel_policy="allow"`` keeps those values and preserves
the treatment BT09 used to avoid bias. ``"zero"`` and ``"mask"`` provide
other policies. None reproduces the separate BT09 cleaning rule for bright
sources, which set only pixels above a threshold measured within half the FWHM
to zero. Significant stellar and
extended emission should therefore be masked upstream.

One g cm\ :sup:`-2` is approximately 4,788.45 solar masses
pc\ :sup:`-2`.

Uncertainty propagation
-----------------------

For :math:`N=I_\mathrm{obs}-I_\mathrm{fg}` and
:math:`D=I_\mathrm{bg}-I_\mathrm{fg}`, the derivatives to first order are:

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{obs}}=-\frac{1}{N},\quad
   \frac{\partial\tau}{\partial I_\mathrm{bg}}=\frac{1}{D},\quad
   \frac{\partial\tau}{\partial I_\mathrm{fg}}=
   \frac{1}{N}-\frac{1}{D}.

``propagate_uncertainty`` combines independent observed, background,
foreground, and opacity terms. ``GTLMapper`` uses an
attached ERR image and an available variance from the foreground model. The
total arrays are written as ``TAU_ERR`` and ``SIGMA_ERR``.

The approximation is masked at saturated pixels and becomes unreliable
near :math:`N=0`. For scientific inference, treat saturation as
censoring and use Monte Carlo or a hierarchical model that can include
foreground/background covariance. Separate opacity, gas/dust ratio,
and background systematics from random pixel noise.

Scientific scope and limits
---------------------------

The package implements a parameterized MIREX calculation for one band and a
registry of filter convolved opacity values. It does not yet
convolve arbitrary throughput curves with raw dust tables or infer a spatially
varying extinction law. NIR/MIR merging, covariance models, and censored Monte
Carlo inference also remain outside the current release. The compatibility default
:math:`\kappa_{8\mu\mathrm{m}}=7.5` cm\ :sup:`2` g\ :sup:`-1` has an
estimated absolute systematic uncertainty of about 30 percent in the
BT09/BT12/KT13 lineage.

For filters outside the registry, users must supply an appropriate
opacity and noise estimate and verify that emission from cold clouds is
negligible. Independent submillimeter calibration and automated
masks for point sources remain future work.
