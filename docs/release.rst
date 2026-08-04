Release and hosting
===================

Repository
----------

The canonical repository is
https://github.com/leongarcia2023/GTLMapping. GTLMapping uses the BSD
3-Clause License.

Build and check a release with:

.. code-block:: console

   python -m build
   python -m twine check dist/*
   check-wheel-contents dist/*.whl

Install the wheel in a clean environment. Check the Python import and run
``gtlmapping --help`` before uploading the package.

Read the Docs
-------------

Read the Docs builds ``docs/conf.py`` through ``.readthedocs.yaml`` and
publishes the result at https://gtlmapping.readthedocs.io/en/latest/.
Warnings fail the build, which keeps broken references out of the public site.

After a documentation change:

1. build the site locally with ``python -m sphinx -W -b html docs docs/_build/html``;
2. push the checked source to the default branch;
3. confirm that Read the Docs completed the build; and
4. open the public site and test the edited pages.

PyPI
----

Publish to TestPyPI first. Install that wheel in a clean environment and run
the quick-start example before publishing the same build to PyPI.

``RELEASE_CHECKLIST.md`` tracks the remaining scientific and authorship
decisions.
