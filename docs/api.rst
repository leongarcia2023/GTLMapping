API reference
=============

Main workflow
-------------

.. autoclass:: gtlmapping.GTLMapper
   :members:

Catalog and geometry
--------------------

.. autofunction:: gtlmapping.read_simon_catalog

.. autofunction:: gtlmapping.find_cloud

Foreground
----------

.. autofunction:: gtlmapping.detect_saturated_samples

.. autofunction:: gtlmapping.estimate_bt12_foreground

.. autofunction:: gtlmapping.fit_conservative_foreground

.. autofunction:: gtlmapping.fit_liberal_foreground

.. autofunction:: gtlmapping.fit_moderate_foreground

.. autofunction:: gtlmapping.interpolate_foreground

.. autofunction:: gtlmapping.cross_validate_foreground

Background
----------

.. autofunction:: gtlmapping.measure_box_background

.. autofunction:: gtlmapping.estimate_box_background

.. autofunction:: gtlmapping.estimate_lmf_background

.. autofunction:: gtlmapping.estimate_smf_background

Extinction and units
--------------------

.. autofunction:: gtlmapping.compute_extinction

.. autofunction:: gtlmapping.propagate_uncertainty

.. autofunction:: gtlmapping.transmission_std

.. autofunction:: gtlmapping.unresolved_transmission

.. autofunction:: gtlmapping.convert_surface_density

Opacity
-------

.. autofunction:: gtlmapping.get_filter_opacity

.. autofunction:: gtlmapping.list_filter_opacities

Result types
------------

.. autoclass:: gtlmapping.CloudEllipse
   :members:

.. autoclass:: gtlmapping.ForegroundSamples
   :members:

.. autoclass:: gtlmapping.InterpolationResult
   :members:

.. autoclass:: gtlmapping.BackgroundResult
   :members:

.. autoclass:: gtlmapping.MappingResult
   :members:

.. autoclass:: gtlmapping.UncertaintyResult
   :members:
