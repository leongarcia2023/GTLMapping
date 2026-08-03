"""Reader for the Simon et al. (2006) MSX IRDC catalog."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from .exceptions import CatalogError
from .models import CloudEllipse


def _optional_float(text: str) -> float | None:
    stripped = text.strip()
    return float(stripped) if stripped else None


def _parse_line(line: str, line_number: int) -> CloudEllipse:
    if len(line.rstrip("\n")) < 61 or line[0:5] != "MSXDC":
        raise CatalogError(f"Malformed Simon catalog record at line {line_number}.")
    padded = line.rstrip("\n").ljust(72)
    try:
        return CloudEllipse(
            name=padded[6:19].strip(),
            component=padded[20:21].strip() or "0",
            glon_deg=float(padded[22:29]),
            glat_deg=float(padded[30:36]),
            major_axis_arcmin=float(padded[37:41]),
            minor_axis_arcmin=float(padded[42:46]),
            pa_deg=float(padded[47:50]),
            area_arcmin2=_optional_float(padded[51:56]),
            peak_contrast=_optional_float(padded[57:61]),
            contrast_snr=_optional_float(padded[62:66]),
            integrated_contrast_arcmin2=_optional_float(padded[67:72]),
        )
    except ValueError as exc:
        raise CatalogError(
            f"Invalid numeric field in Simon catalog at line {line_number}."
        ) from exc


def read_simon_catalog(path: str | Path) -> list[CloudEllipse]:
    """Read ``catalog.dat`` using its published fixed-width definition."""

    catalog_path = Path(path)
    if not catalog_path.is_file():
        raise CatalogError(f"Catalog file does not exist: {catalog_path}")

    entries: list[CloudEllipse] = []
    with catalog_path.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            entries.append(_parse_line(line, line_number))
    if not entries:
        raise CatalogError(f"No records were found in {catalog_path}.")
    return entries


def _canonical_name(name: str) -> str:
    value = name.strip()
    if value.startswith("MSXDC "):
        value = value[6:]
    if not value.startswith("G"):
        value = f"G{value}"
    return value


def find_cloud(
    catalog: str | Path | Iterable[CloudEllipse],
    name: str,
    *,
    component: str = "0",
) -> CloudEllipse:
    """Find one cloud/core by catalog name.

    ``name`` may be written as ``G028.37+00.07``,
    ``028.37+00.07``, or ``MSXDC G028.37+00.07``.
    """

    entries = (
        read_simon_catalog(catalog)
        if isinstance(catalog, (str, Path))
        else list(catalog)
    )
    target = _canonical_name(name)
    matches = [
        entry
        for entry in entries
        if _canonical_name(entry.name) == target and entry.component == component
    ]
    if len(matches) != 1:
        raise CatalogError(
            f"Expected one match for {target!r} component {component!r}; "
            f"found {len(matches)}."
        )
    return matches[0]
