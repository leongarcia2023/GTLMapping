"""Filter-convolved dust opacities used by MIREX mapping."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class FilterOpacity:
    """One filter-convolved opacity normalized to a reference gas/dust ratio."""

    filter_name: str
    wavelength_um: float
    dust_model: str
    kappa_reference_cm2_g: float
    reference_gas_to_dust_ratio: float = 156.0

    def at_gas_to_dust_ratio(self, gas_to_dust_ratio: float) -> float:
        """Return opacity per total mass for another gas-to-dust ratio."""

        ratio = float(gas_to_dust_ratio)
        if ratio <= 0:
            raise ValueError("gas_to_dust_ratio must be positive.")
        return (
            self.kappa_reference_cm2_g
            * self.reference_gas_to_dust_ratio
            / ratio
        )


_FILTER_ALIASES = {
    "IRAC2": "IRAC2",
    "IRACCH2": "IRAC2",
    "SPITZERIRAC2": "IRAC2",
    "F480M": "F480M",
    "JWSTF480M": "F480M",
    "NIRCAMF480M": "F480M",
    "JWSTNIRCAMF480M": "F480M",
    "IRAC4": "IRAC4",
    "IRACCH4": "IRAC4",
    "SPITZERIRAC4": "IRAC4",
    "F770W": "F770W",
    "JWSTF770W": "F770W",
    "MIRIF770W": "F770W",
    "F2100W": "F2100W",
    "JWSTF2100W": "F2100W",
    "MIRIF2100W": "F2100W",
}

_WAVELENGTHS_UM = {
    "IRAC2": 4.5,
    "F480M": 4.8,
    "IRAC4": 8.0,
    "F770W": 7.7,
    "F2100W": 21.0,
}

# Filter convolved OH94 and WD01 values used by the included benchmarks.
# Every row is normalized to gas to dust ratio 156.
_REFERENCE_OPACITIES = {
    "wd01_rv31": {
        "IRAC2": 6.67,
        "F480M": 5.85,
        "IRAC4": 8.02,
        "F770W": 5.93,
        "F2100W": 6.30,
    },
    "wd01_rv55": {
        "IRAC2": 10.80,
        "F480M": 9.50,
        "IRAC4": 10.47,
        "F770W": 7.98,
        "F2100W": 6.91,
    },
    "wd01_rv55_case_b": {
        "IRAC2": 12.15,
        "F480M": 11.43,
        "IRAC4": 13.20,
        "F770W": 10.96,
        "F2100W": 7.58,
    },
    "oh94_thin_ice_0yr": {
        "IRAC2": 8.62,
        "F480M": 7.90,
        "IRAC4": 6.60,
        "F770W": 6.04,
        "F2100W": 5.83,
    },
    "oh94_thin_ice_coagulated": {
        "IRAC2": 10.71,
        "F480M": 9.76,
        "IRAC4": 7.80,
        "F770W": 7.25,
        "F2100W": 7.86,
    },
}

_MODEL_ALIASES = {
    "OH94": "oh94_thin_ice_coagulated",
    "OH94THINICE": "oh94_thin_ice_coagulated",
    "OH94THINICECOAGULATED": "oh94_thin_ice_coagulated",
    "OH94THINICE1E5YR": "oh94_thin_ice_coagulated",
    "OH94THINICE0YR": "oh94_thin_ice_0yr",
    "WD01RV31": "wd01_rv31",
    "WD01RV55": "wd01_rv55",
    "WD01RV55CASEB": "wd01_rv55_case_b",
}


def _compact(value: str) -> str:
    return "".join(character for character in value.upper() if character.isalnum())


def _canonical_filter(filter_name: str) -> str:
    try:
        return _FILTER_ALIASES[_compact(filter_name)]
    except KeyError as exc:
        choices = ", ".join(_WAVELENGTHS_UM)
        raise KeyError(
            f"Unknown filter {filter_name!r}; available filters are {choices}."
        ) from exc


def _canonical_model(dust_model: str) -> str:
    if dust_model in _REFERENCE_OPACITIES:
        return dust_model
    try:
        return _MODEL_ALIASES[_compact(dust_model)]
    except KeyError as exc:
        choices = ", ".join(_REFERENCE_OPACITIES)
        raise KeyError(
            f"Unknown dust model {dust_model!r}; available models are {choices}."
        ) from exc


def get_filter_opacity(
    filter_name: str,
    *,
    dust_model: str = "oh94_thin_ice_coagulated",
    gas_to_dust_ratio: float = 156.0,
) -> float:
    """Return a filter-convolved opacity in cm² g⁻¹ of total mass.

    The tabulated values use a gas-to-dust ratio of 156. Opacity per total
    mass scales inversely with the requested ratio, so the supplied F480M
    value is 9.76 cm² g⁻¹ at 156 and 15.2256 cm² g⁻¹ at 100.
    """

    filter_key = _canonical_filter(filter_name)
    model_key = _canonical_model(dust_model)
    record = FilterOpacity(
        filter_name=filter_key,
        wavelength_um=_WAVELENGTHS_UM[filter_key],
        dust_model=model_key,
        kappa_reference_cm2_g=_REFERENCE_OPACITIES[model_key][filter_key],
    )
    return record.at_gas_to_dust_ratio(gas_to_dust_ratio)


def list_filter_opacities(
    *,
    dust_model: str = "oh94_thin_ice_coagulated",
    gas_to_dust_ratio: float = 156.0,
) -> dict[str, float]:
    """Return all supported filter opacities for one model and normalization."""

    model_key = _canonical_model(dust_model)
    return {
        filter_name: get_filter_opacity(
            filter_name,
            dust_model=model_key,
            gas_to_dust_ratio=gas_to_dust_ratio,
        )
        for filter_name in _WAVELENGTHS_UM
    }
