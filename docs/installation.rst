Installation
============

Current repository
------------------

Clone or download the repository, enter its top-level directory, and
install it with:

.. code-block:: console

   git clone https://github.com/leongarcia2023/GTLMapping.git
   cd GTLMapping
   python -m pip install .

For tests and documentation:

.. code-block:: console

   python -m pip install -e ".[dev]"

Science data are intentionally not bundled in the Python distribution.
Provide the FITS image and ``catalog.dat`` (or another catalog path)
separately when running GTLMapping.

Optional alignment support
--------------------------

Strict WCS checking is part of the base package. Automatic reprojection
is optional:

.. code-block:: console

   python -m pip install "GTLMapping[align]"

Future PyPI release
-------------------

After the package has been published, the intended command is:

.. code-block:: console

   python -m pip install GTLMapping

Do not advertise the PyPI command as active until the first release is
actually published.

Requirements
------------

The package requires Python 3.10 or later, NumPy, SciPy, Astropy, and
PyKrige.
