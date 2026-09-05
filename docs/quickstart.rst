Quick start
===========

A complete map
--------------

This example selects Cloud C from the Simon catalog, fits the conservative
foreground, loads a prepared background, and writes the result:

.. code-block:: python

   from gtlmapping import GTLMapper

   mapper = GTLMapper.from_fits("1kx1k.fits")
   cloud = mapper.select_cloud("catalog.dat", "G028.37+00.07")

   mapper.detect_foreground(noise_sigma=0.6)
   mapper.fit_foreground(
       method="conservative",
       noise_sigma=0.6,
   )
   mapper.set_background_from_fits("SMFbg1.fits")

   result = mapper.compute(
       kappa_cm2_g=7.5,
       bright_pixel_policy="allow",
   )
   result.write("cloud_c_gtl.fits")

``result`` keeps the arrays in memory. The FITS file stores the same arrays in
named extensions:

.. list-table::
   :header-rows: 1

   * - Extension
     - Contents
   * - ``SIGMA``
     - Mass surface density in g cm\ :sup:`-2`
   * - ``TAU``
     - Optical depth
   * - ``FOREGROUND`` / ``BACKGROUND``
     - The intensity surfaces used in the transfer equation
   * - ``SATURATED`` / ``INVALID_BG`` / ``BRIGHT``
     - Masks for nonpositive transmission, invalid backgrounds, and negative
       optical depth
   * - ``UNRESOLVED`` / ``TRANS_LIM``
     - Full sensitivity mask and the adopted transmission threshold
   * - ``SIGMA_ERR`` / ``TAU_ERR``
     - Uncertainty propagated to first order, when the calculation has
       uncertainty inputs
   * - ``FG_CONSTRAINT``
     - Pixels adjusted by a requested foreground feasibility constraint

Check the background grid
-------------------------

Array shape alone cannot establish that two images cover the same sky.
``set_background_from_fits`` compares their celestial coordinates at the
center and corners. A mismatch raises ``GridMismatchError``.

.. code-block:: python

   mapper.set_background_from_fits("smf_background.fits")

If reprojection belongs in the analysis, install ``GTLMapping[align]`` and
request it:

.. code-block:: python

   mapper.set_background_from_fits(
       "smf_background.fits",
       align=True,
   )

Choose a foreground profile
---------------------------

GTLMapping separates the published BT12 reference from three spatial
profiles:

.. list-table::
   :header-rows: 1

   * - Method
     - Behavior
     - Use
   * - ``bt12``
     - One foreground value from minima passing a partner-separation test
     - Reference calculation
   * - ``conservative``
     - Broad spatial trend with a pointwise BT12 floor
     - Default spatial comparison
   * - ``moderate``
     - Quadratic trend ordered above conservative GTL
     - Intermediate sensitivity test
   * - ``liberal``
     - Sample-driven quadratic ordered above moderate GTL
     - Permissive sensitivity test

The profiles enforce ``BT12 <= conservative <= moderate <= liberal``
pointwise. At a fixed background and opacity, their surface densities follow
the same order on jointly valid, uncensored pixels. Moderate and liberal use
progressively larger explicit censoring budgets.

Moderate GTL allows 0.5% of pixels near saturation and caps strict censoring
at 0.01%:

.. code-block:: python

   mapper.detect_foreground(noise_sigma=0.6)
   mapper.fit_foreground(
       method="moderate",
       noise_sigma=0.6,
   )
   mapper.set_background_from_fits("smf_background.fits")
   result = mapper.compute_moderate(kappa_cm2_g=7.5)

Liberal GTL gives the accepted minima more influence. Its default budgets are
1% near saturation and 0.1% strict censoring:

.. code-block:: python

   mapper.detect_foreground(noise_sigma=0.6)
   mapper.fit_foreground(
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

``compute_moderate`` and ``compute_liberal`` use a positive floor on
transmitted intensity for censored pixels. The same default applies to ``compute`` when a threshold can be inferred.
Use ``UNRESOLVED`` to identify all limits, including weak positive residuals;
``SATURATED`` records only nonpositive residuals. ``result.detection_mask``
selects positive, resolved columns. These finite values are conditional
sensitivity scales, not calibrated confidence bounds.

Use BT12 for the published comparison:

.. code-block:: python

   mapper.fit_foreground(
       method="bt12",
       cloud=cloud,
       noise_sigma=0.6,
   )

BT12 searches the cloud ellipse for pixels within :math:`2\sigma` of the
global minimum. It requires a second pixel at least 8 arcsec away, then sets
the foreground to the mean qualifying intensity minus :math:`2\sigma`.

Legacy interpolation
--------------------

Kriging, RBF, spline, Gaussian, Cauchy, and ``flat`` remain available for
method studies. They fit the accepted minima directly and can transfer cloud
structure into the foreground. The package clips legacy interpolation to the
sample range unless you disable that safeguard.

The following settings reproduce the legacy kriging calculation, which
retained repeated coordinates:

.. code-block:: python

   mapper.detect_foreground(noise_sigma=0.6)
   mapper.fit_foreground(
       method="kriging",
       foreground_margin=1.0,
       clip_to_sample_range=False,
       kriging_duplicate_policy="repeat",
   )

The default scan uses 5 by 5 windows and a stride of half a box.
``cover_edges=True`` adds final windows at uncovered edges.
``restrict_to_cloud=True`` limits detection to the selected ellipse. Both
settings change the legacy sample set.

JWST F480M with an ERR image
----------------------------

JWST products often store the science image and calibrated uncertainty in
separate extensions. In this fragment, supply your field's ``target_mask``,
``adjacent_pixel_boxes`` and ``background_uncertainty``. A complete runnable
case is provided in :doc:`jwst_sgrc`.

.. code-block:: python

   import numpy as np

   mapper = GTLMapper.from_fits(
       "F480M_registered.fits",
       hdu="SCI",
       uncertainty_hdu="ERR",
   )

   noise_sigma = float(np.nanmedian(mapper.observed_std[target_mask]))
   mapper.fit_foreground(
       method="bt12",
       region_mask=target_mask,
       noise_sigma=noise_sigma,
       min_separation_arcsec=0.74,
   )
   mapper.estimate_background(
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

The opacity registry returns 9.76 cm\ :sup:`2` g\ :sup:`-1` at gas-to-dust
ratio 156 and 15.2256 cm\ :sup:`2` g\ :sup:`-1` at ratio 100. Include the
gas-mass normalization in reported results. These rounded opacities are
adopted inputs, not independently reproduced filter convolutions.

Handle conflicts between model surfaces
---------------------------------------

The transfer equation has no physical solution where
:math:`I_\mathrm{fg} \ge I_\mathrm{bg}`. Inspect the foreground and
background before applying a constraint. If a positive floor on transmitted
intensity is justified, record the projection in the output:

.. code-block:: python

   mapper.constrain_foreground(
       minimum_transmitted_intensity=2 * noise_sigma,
   )
   result = mapper.compute(
       filter_name="F480M",
       saturation_policy="lower_limit",
       intensity_floor=2 * noise_sigma,
   )

``FG_CONSTRAINT`` identifies adjusted foreground pixels. ``UNRESOLVED``
marks all sensitivity limits, including weak positive transmission;
``SATURATED`` marks only strict zero crossings. Exclude unresolved pixels
from ordinary Gaussian fits.
