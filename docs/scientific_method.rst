Scientific method
===================

Radiative transfer
--------------------

The unobscured intensity includes the foreground. With transmitted intensity
N = I_obs - I_fg and D = I_bg - I_fg, GTLMapping computes

.. math::

   \tau=-\ln(N/D), \qquad \Sigma=\tau/\kappa.

A detected absorption signal requires D > N > 0 and enough sensitivity to
measure N. The compatibility opacity is 7.5 cm²/g. Its mass convention must
match the reported surface density. On jointly detected pixels,

.. math::

   \frac{\partial\tau}{\partial I_\mathrm{fg}}=\frac{1}{N}-\frac{1}{D}>0.

A higher foreground therefore raises optical depth. This does not establish
that it is more accurate, and does not order stored values once a pixel
becomes a limit.

Candidate foreground sites
----------------------------

The scan searches overlapping windows across the input image. It accepts a
window minimum when a second eligible pixel lies within twice the noise of
that minimum and beyond the configured angular separation. The default is
4 arcsec for a 2 arcsec beam; the IRAC comparisons use 8 arcsec for BT12.

A separated partner does not make all sites mutually independent. A minimum
can still contain transmitted background. Tests without opaque cores
show the resulting bias. Repeated coordinates are merged and their
multiplicities retained as fit weights; overlapping windows are not
independent measurements of the same site.

The default scan leaves incomplete edge windows unused. ``cover_edges=True``
includes a final valid origin. Use ``restrict_to_cloud=True`` or an explicit
``region_mask`` to limit the search. Otherwise, changing the crop can change
the samples even when the catalog ellipse stays fixed.

Planes require three independent design columns and quadratics require six.
Rank-deficient fits raise ``InsufficientSamplesError``. Diagnostics record
the condition number and the numbers of sites inside and outside the fit
region. Full rank alone does not establish adequate spatial coverage.

BT12 and the spatial profiles
-------------------------------

BT12 finds the minimum inside the selected aperture. It selects pixels
strictly between the minimum and the minimum plus twice the noise, requires
a qualifying pixel beyond the angular separation, and subtracts twice the
noise from their mean. This estimator is separate from a constant fit to
GTL window minima.

Conservative GTL fits a plane with a soft-L1 loss, subtracts the foreground
margin, and retains positive enhancement above BT12. It chooses the largest
allowed blend on a 101-point grid.

Moderate GTL fits a quadratic trend, clips extrapolation to the sample
intensity range, and applies a 50% two-sided BT12 pull. Its final foreground
retains positive enhancement above conservative GTL. Liberal GTL uses the
same quadratic family without that two-sided pull, retaining positive
enhancement above moderate GTL. Both reduce the enhancement to meet their
budgets.

.. list-table:: Default fitting budgets
   :header-rows: 1

   * - Profile
     - Near saturation
     - Strict saturation
   * - Conservative
     - 1%
     - 0%
   * - Moderate
     - 0.5%
     - 0.01%
   * - Liberal
     - 1%
     - 0.1%

Near saturation means N at or below twice the adopted noise; strict
saturation means N at or below zero. Counts already present in the preceding
reference are allowed. These budgets can become active fitting ceilings.
Their achieved counts are not independent discoveries of opaque gas.

The foreground hierarchy is imposed:

.. math::

   I_\mathrm{fg}^{BT12}\leq I_\mathrm{fg}^{cons}
   \leq I_\mathrm{fg}^{mod}\leq I_\mathrm{fg}^{lib}.

All named GTL profiles require a usable BT12 reference, including liberal
GTL. They are sensitivity presets, not confidence levels. Kriging, RBF,
splines and kernel interpolation remain available for experiments.
Numerical stability does not establish a physical foreground.

Detections and unresolved transmission
----------------------------------------

A residual near zero cannot be treated as a measured column simply because
a logarithm can be evaluated. GTLMapping writes separate masks:

* ``SATURATED`` retains the strict condition N at or below zero.
* ``UNRESOLVED`` identifies transmission at or below the adopted threshold,
  including positive residuals.
* ``TRANS_LIM`` records the threshold in MJy/sr.

With ``saturation_policy="lower_limit"``, unresolved pixels receive

.. math::

   \Sigma_\mathrm{lim}=-\frac{1}{\kappa}
   \ln\left(\frac{I_\mathrm{trans,min}}{D}\right).

This is a conditional sensitivity limit, not a calibrated confidence bound.
A wrong foreground or background can invalidate its interpretation.
``saturation_policy="mask"`` masks these pixels instead.
Use ``result.detection_mask`` for positive detections only.

``GTLMapper.compute`` defaults to finite limits when sensitivity information
is available. The threshold comes from an explicit ``detection_threshold``,
an ``intensity_floor``, or twice the supplied ``residual_std``, observed
uncertainty, or fit noise, in that order. ``detection_sigma`` changes the
multiplier for standard deviations. The observed-noise default omits
foreground uncertainty; metadata state the basis. Without noise information,
only strict masking is possible.

The low-level ``compute_extinction`` does not infer a noise model. Pass a
floor or threshold for measured-column analyses. A floor also sets the
threshold unless a different threshold is supplied. Document both if they
differ. ``compute_moderate`` and ``compute_liberal`` project the foreground
below the background by the floor and retain finite limits.
``FG_CONSTRAINT`` records adjustments. Missing data and backgrounds that
cannot support a physical inversion remain invalid.

Background and source masks
-----------------------------

LMF samples a large trimmed median on a grid spaced by 24 arcsec. SMF
measures local medians outside the cloud, with a filter one third of the
catalog major axis, and fills the cloud using inverse-square-distance
weights within one semimajor axis. An aligned external background can be
supplied instead.

The box estimator averages adjacent-box medians after user-specified bright
rejection. The Sgr C example uses 15 MJy/sr. Choose boxes from the relevant
field and retain their coordinates and masks.

For observed intensity above the background, optical depth is negative.
``bright_pixel_policy="allow"`` retains it; ``"zero"`` and ``"mask"`` are
alternatives. This rule does not identify all point sources. Apply a
separate source mask when the analysis requires one.

Opacity and mass convention
-----------------------------

The registry contains adopted, rounded inputs for five dust prescriptions
and five filters. It does not independently reproduce their convolutions.
Archive the dust table, extinction-versus-absorption choice, weighting
spectrum, throughput convention and normalization for a calibrated opacity.

The registry adopts a gas-mass convention at gas/dust ratio R = 156:

.. math::

   \kappa_\mathrm{gas}(R)=\kappa_\mathrm{gas}(156)\frac{156}{R},
   \qquad
   \kappa_\mathrm{total}(R)=\kappa_\mathrm{gas}(156)\frac{156}{R+1}.

The adopted F480M value is 9.76 at ratio 156 and 15.2256 at ratio 100 per
gram of gas. This is not exact total-mass scaling. Use ``mass_basis="total"``
for gas plus dust. A reference already per total mass instead rescales by
(R0+1)/(R+1); ``FilterOpacity`` supports that explicit reference basis.
Direct opacities must use the mass basis declared by the caller.

Uncertainty
-------------

For a supplied foreground error and same-pixel covariance:

.. code-block:: python

   from gtlmapping.extinction import transmission_std

   residual_error = transmission_std(
       observed_std, foreground_std, covariance=observed_foreground_covariance
   )
   result = mapper.compute(residual_std=residual_error)

The helper evaluates

.. math::

   \sigma_N^2=\sigma_\mathrm{obs}^2+\sigma_\mathrm{fg}^2
       -2\,\mathrm{Cov}(I_\mathrm{obs},I_\mathrm{fg}).

It does not estimate covariance from the data. Analytic ``TAU_ERR`` and
``SIGMA_ERR`` combine independent first-order intensity and opacity terms
and are masked on unresolved pixels. The trend variance is conditional on
the selected sites; it omits selection uncertainty and BT12 reference error.

Integrated inference must preserve shared systematics and spatial noise
correlation. Draw background and opacity parameters jointly per realization,
repeat sample selection and foreground fitting, and separate detected mass
from unresolved contributions. The package does not claim calibrated
censored posterior intervals. Validation reports sensitivity and recovery
errors, not a pixel-noise error bar on a cloud mass.
