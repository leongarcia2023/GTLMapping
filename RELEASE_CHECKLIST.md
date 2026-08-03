# Public release checklist

The code is structured for GitHub, PyPI, and Read the Docs, but these
decisions still require the project owner's approval:

- [x] Choose and add the BSD 3-Clause open-source license.
- [x] Confirm the canonical GitHub organization/repository URL.
- [ ] Confirm the public author list, affiliations, and contact route.
- [ ] Decide whether `catalog.dat` may be redistributed or should be
      downloaded separately from CDS/VizieR.
- [ ] Obtain scientific review of the four-beam-radius independence
      criterion versus the 8-arcsecond BT12 criterion.
- [x] Trace `1kx1k.fits` to the parent GLIMPSE IRAC-4 MJy/sr mosaic.
- [ ] Add the recovered `BUNIT=MJy/sr` keyword to a canonical input
      product outside the read-only project mirror.
- [ ] Define and validate point-source, extended-emission, and
      instrumental-artifact masks.
- [ ] Validate the foreground model and uncertainty on multiple IRDCs.
- [x] Implement transparent first-order kriging/background/opacity
      uncertainty propagation.
- [ ] Add censored Monte Carlo validation near saturation and document
      covariance/systematic priors.
- [ ] Obtain machine-readable Sgr C background boxes/apertures for an
      exact Rubén mass reproduction.
- [ ] Add a software DOI/citation after creating a Zenodo release.
- [x] Configure test, documentation, build, metadata, and wheel-smoke checks in CI.
- [x] Confirm the first GitHub CI run passes.
- [x] Build and inspect both wheel and source distribution locally.
- [x] Run `twine check` and `check-wheel-contents` on clean distributions.
- [ ] Publish to TestPyPI before PyPI.
- [x] Connect the GitHub repository to Read the Docs.
