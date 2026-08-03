from __future__ import annotations

import pytest

from gtlmapping.opacity import get_filter_opacity, list_filter_opacities


def test_f480m_opacity_tracks_gas_to_dust_normalization() -> None:
    assert get_filter_opacity("F480M") == pytest.approx(9.76)
    assert get_filter_opacity(
        "JWST/NIRCam F480M",
        gas_to_dust_ratio=100.0,
    ) == pytest.approx(15.2256)


def test_oh94_filter_table_contains_spitzer_and_jwst_values() -> None:
    table = list_filter_opacities()

    assert table["IRAC4"] == pytest.approx(7.80)
    assert table["F770W"] == pytest.approx(7.25)
    assert table["F2100W"] == pytest.approx(7.86)
