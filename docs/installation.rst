Installation
============

Install from GitHub
-------------------

GTLMapping has not reached PyPI yet. Install the current release from its
GitHub repository:

.. code-block:: console

   git clone https://github.com/leongarcia2023/GTLMapping.git
   cd GTLMapping
   python -m pip install .

Use an editable install for development, tests, and documentation:

.. code-block:: console

   python -m pip install -e ".[dev]"

GTLMapping requires Python 3.10 or later. NumPy, SciPy, Astropy, and PyKrige
are installed as dependencies.

Bring your own data
-------------------

The Python distribution does not contain science images or ``catalog.dat``.
Pass their paths to GTLMapping when you run an analysis. Keeping the images
outside the package avoids a large download and prevents local paths from
entering a release.

Optional reprojection
---------------------

The base package checks whether the observed image and background share a
celestial grid. Install the ``align`` extra if you want GTLMapping to reproject
a mismatched background:

.. code-block:: console

   python -m pip install "GTLMapping[align]"

Pass ``align=True`` only when the reprojection is part of your analysis.

PyPI release
------------

After the first PyPI release, installation will use:

.. code-block:: console

   python -m pip install GTLMapping

Until then, use the GitHub command above.
