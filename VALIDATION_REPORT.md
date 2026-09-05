# Validation

The [validation page](docs/validation.rst) reports software checks and a
synthetic recovery suite with 220 images and 880 attempted fits. It includes
failed fits and cases where a spatial foreground worsens the inferred column.

The [Cloud C/F/H comparisons](docs/cloud_comparisons.rst) and
[Sagittarius C benchmark](docs/jwst_sgrc.rst) use a shared transmission
threshold within each comparison. Their tables distinguish detections,
finite sensitivity limits, and nonpositive transmission. Increasing a
finite-map sum does not establish better physical recovery.

Machine-readable inputs to the tables and figures are in `docs/_static`.
The examples regenerate comparisons from user-supplied FITS files. The
[method documentation](docs/scientific_method.rst) describes the opacity
convention, sample support, active constraints and uncertainty limitations.
