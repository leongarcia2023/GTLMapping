Quick start
===========

The high-level interface keeps every stage available for inspection:

.. code-block:: python

   from gtlmapping import GTLMapper

   mapper = GTLMapper.from_fits("1kx1k.fits")
   cloud = mapper.select_cloud("catalog.dat", "G028.37+00.07")

   samples = mapper.detect_foreground(
       noise_sigma=0.6,
   )
   foreground = mapper.fit_foreground(
       method="conservative",
       noise_sigma=0.6,
   )
   background = mapper.set_background_from_fits("SMFbg1.fits")

   result = mapper.compute(
       kappa_cm2_g=7.5,
       bright_pixel_policy="allow",
   )
   result.write("cloud_c_gtl.fits")

The output contains:

``SIGMA``
   Mass surface density in g cm\ :sup:`-2`.

``TAU``
   Optical depth.

``FOREGROUND`` and ``BACKGROUND``
   The two modeled intensity surfaces.

``SATURATED``, ``INVALID_BG``, and ``BRIGHT``
   Diagnostic masks.

``SIGMA_ERR`` and ``TAU_ERR``
   Present when any uncertainty source is supplied.

``FG_CONSTRAINT``
   Present when the interpolated foreground was explicitly constrained
   against the observed background.

Using a prepared background
---------------------------

.. code-block:: python

   mapper.set_background_from_fits("smf_background.fits")

Matching array dimensions are not enough: the method compares the
celestial grid at the center and corners. A mismatch raises
``GridMismatchError``.

If reprojection is intended:

.. code-block:: python

   mapper.set_background_from_fits(
       "smf_background.fits",
       align=True,
   )

Install the ``align`` extra first.

Foreground alternatives
-----------------------

``conservative`` is the spatially varying default. It fits a robust
planar foreground candidate to the GTL window samples and applies BT12
as a hard pointwise floor. Saturation-fraction guardrails shrink the
spatial enhancement when necessary and, by default, permit no new
strictly saturated pixels in the fitted region.

``moderate`` is the middle profile. It uses a robust quadratic GTL trend with
a 50% soft BT12 anchor, a 0.5% near-saturation budget, and a 0.01% strict
lower-limit ceiling.

.. code-block:: python

   samples = mapper.detect_foreground(noise_sigma=0.6)
   foreground = mapper.fit_foreground(
       method="moderate",
       noise_sigma=0.6,
   )
   mapper.set_background_from_fits("smf_background.fits")
   result = mapper.compute_moderate(kappa_cm2_g=7.5)

``liberal`` is the controlled, more permissive alternative. It fits a robust
quadratic trend directly to the GTL samples, uses no hard BT12 floor, and has
separate budgets for pixels within :math:`2\sigma` of the foreground and
pixels strictly below it. The default soft BT12-anchor weight is zero.

.. code-block:: python

   samples = mapper.detect_foreground(noise_sigma=0.6)
   foreground = mapper.fit_foreground(
       method="liberal",
       noise_sigma=0.6,
       target_local_saturation_fraction=0.01,
       maximum_strict_saturation_fraction=0.001,
   )
   mapper.set_background_from_fits("smf_background.fits")
   result = mapper.compute_liberal(
       kappa_cm2_g=7.5,
       bright_pixel_policy="allow",
   )

``compute_liberal`` applies the foreground/background feasibility projection
and a :math:`2\sigma` transmitted-intensity floor. Censored pixels therefore
receive finite surface-density lower limits while remaining marked in
``SATURATED``; adjusted foreground pixels remain marked in ``FG_CONSTRAINT``.
They must not be analyzed as ordinary detections.

The original ``kriging``, ``rbf``, ``spline``, ``gaussian``, and
``cauchy`` interpolators remain available for experiments and notebook
reproduction. They fit local minima directly and can imprint cloud
structure into the foreground. ``flat`` is a constant summary of the
GTL samples; it is not the published BT12 method.

For the paper-faithful BT12 comparison:

.. code-block:: python

   foreground = mapper.fit_foreground(
       method="bt12",
       cloud=cloud,
       noise_sigma=0.6,
   )

BT12 searches inside the cloud ellipse for pixels within
:math:`2\sigma` of the global minimum, requires an independent pixel at
least 8 arcsec away, and sets the constant foreground to the mean of
all qualifying pixels minus :math:`2\sigma`.

Legacy GTL interpolation is clipped to the observed sample range by
default to prevent unphysical RBF/spline extrapolation. Set
``clip_to_sample_range=False`` only for a deliberate compatibility run.

Reproducing the prototype notebook
----------------------------------

The default GTL scan reproduces the notebook's 5-by-5, half-box-stride,
full-image search. Stable kriging aggregates duplicate coordinates and
uses a pseudo-inverse. To reproduce the notebook's repeated-coordinate
foreground surface exactly, use:

.. code-block:: python

   samples = mapper.detect_foreground(noise_sigma=0.6)
   foreground = mapper.fit_foreground(
       method="kriging",
       foreground_margin=1.0,
       clip_to_sample_range=False,
       kriging_duplicate_policy="repeat",
   )

The compatibility run uses a :math:`2\sigma=1.2` margin unless
overridden. ``cover_edges=True`` adds final scan windows that the
notebook omitted; ``restrict_to_cloud=True`` limits GTL detection to
the selected cloud. Both are explicit opt-ins because they change the
historical sample set.

JWST F480M with uncertainty
---------------------------

The Sgr C file stores intensity and calibrated uncertainty in separate
extensions:

.. code-block:: python

   import numpy as np
   mapper = GTLMapper.from_fits(
       "F480M_registered.fits",
       hdu="SCI",
       uncertainty_hdu="ERR",
   )

   noise_sigma = float(np.nanmedian(mapper.observed_std[target_mask]))
   foreground = mapper.fit_foreground(
       method="bt12",
       region_mask=target_mask,
       noise_sigma=noise_sigma,
       min_separation_arcsec=0.74,
   )
   background = mapper.estimate_background(
       method="boxes",
       boxes=adjacent_pixel_boxes,
       maximum_intensity=15.0,
   )
   result = mapper.compute(
       filter_name="F480M",
       gas_to_dust_ratio=156,
       background_std=background_uncertainty,
       kappa_std_cm2_g=0.30 * 9.76,
   )

The opacity lookup returns 9.76 cm\ :sup:`2` g\ :sup:`-1` for gas/dust
156 and 15.2256 cm\ :sup:`2` g\ :sup:`-1` for gas/dust 100. State the
normalization in every result rather than substituting one number for
the other.

Handling model-conflict pixels
------------------------------

A finite kriging prediction does not guarantee a physical radiative
transfer solution. If :math:`I_\mathrm{fg} \ge I_\mathrm{bg}`, inspect
the foreground and background models first. An explicit constrained
run is available:

.. code-block:: python

   mapper.constrain_foreground(
       minimum_transmitted_intensity=2 * noise_sigma,
   )
   result = mapper.compute(
       filter_name="F480M",
       saturation_policy="lower_limit",
       intensity_floor=2 * noise_sigma,
   )

Constraint and saturation masks remain in memory and in the FITS
output. Lower-limit pixels are censored measurements; do not include
them in ordinary Gaussian fitting.
