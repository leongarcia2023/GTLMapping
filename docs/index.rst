GTLMapping
==========

GTLMapping makes mid-infrared extinction maps with either the constant
foreground of Butler & Tan (2012) or a foreground that varies across the
image. It writes every modeled surface and diagnostic mask to FITS, so you can
trace a mass estimate back to the image, foreground, background, and opacity
that produced it.

The package offers three spatial profiles. ``conservative`` keeps BT12 as a
pointwise foreground floor. ``moderate`` allows a broader fluctuation with
tight censoring limits. ``liberal`` gives the local minima more influence and
belongs in sensitivity analyses.

.. warning::

   GTLMapping is alpha scientific software. The code has passed its software
   and numerical tests. The spatial foreground models need validation on more
   infrared dark clouds before they support general astrophysical claims.

Where to start
--------------

Install the package, then work through one complete map in :doc:`quickstart`.
The :doc:`scientific_method` page defines the radiative-transfer calculation
and the assumptions behind each foreground profile. :doc:`jwst_sgrc` records
the controlled F480M comparison against BT12.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Method and evidence

   scientific_method
   jwst_sgrc
   validation
   references

.. toctree::
   :maxdepth: 2
   :caption: Package reference

   api
   release
