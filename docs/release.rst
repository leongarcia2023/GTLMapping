Release and hosting
===================

GitHub
------

Before the first public push, confirm the public author metadata. The
software is distributed under the BSD 3-Clause License, and the canonical repository is
https://github.com/leongarcia2023/GTLMapping.

Build artifacts locally:

.. code-block:: console

   python -m build
   python -m twine check dist/*
   check-wheel-contents dist/*.whl

Install the wheel into a clean environment and exercise both the Python
import and ``gtlmapping --help`` before uploading it.

Read the Docs
-------------

The repository includes ``.readthedocs.yaml``. After pushing to GitHub:

1. import the repository into Read the Docs;
2. build the default branch;
3. confirm that warnings fail the build; and
4. set the canonical documentation URL in the repository metadata.

PyPI
----

Publish to TestPyPI first, install the resulting wheel into a clean
environment, run the quick start, and only then publish to PyPI.

The repository-level ``RELEASE_CHECKLIST.md`` tracks the remaining
scientific and governance decisions.
