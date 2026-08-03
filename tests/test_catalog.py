from __future__ import annotations

from pathlib import Path

import pytest

from gtlmapping.catalog import find_cloud, read_simon_catalog
from gtlmapping.exceptions import CatalogError


CATALOG_LINE = (
    "MSXDC G028.37+00.07 0  28.373  0.076 12.0  9.3  78 41.80 0.61      10.91\n"
)


def test_fixed_width_catalog_and_flexible_lookup(tmp_path: Path) -> None:
    path = tmp_path / "catalog.dat"
    path.write_text(CATALOG_LINE, encoding="ascii")

    entries = read_simon_catalog(path)
    cloud = find_cloud(entries, "028.37+00.07")

    assert len(entries) == 1
    assert cloud.name == "G028.37+00.07"
    assert cloud.is_cloud
    assert cloud.major_axis_arcmin == pytest.approx(12.0)
    assert cloud.minor_axis_arcmin == pytest.approx(9.3)
    assert cloud.contrast_snr is None


def test_missing_catalog_match_is_explicit(tmp_path: Path) -> None:
    path = tmp_path / "catalog.dat"
    path.write_text(CATALOG_LINE, encoding="ascii")

    with pytest.raises(CatalogError, match="found 0"):
        find_cloud(path, "G999.99+99.99")
