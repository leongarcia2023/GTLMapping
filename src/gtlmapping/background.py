"""BT09-style large- and small-scale background estimators."""

from __future__ import annotations

import numpy as np
from astropy.wcs import WCS
from scipy.interpolate import griddata
from scipy.signal import fftconvolve

from .geometry import ellipse_mask, pixel_scale_arcsec
from .models import BackgroundResult, CloudEllipse


def _grid_positions(length: int, step: int) -> np.ndarray:
    positions = list(range(0, length, step))
    if positions[-1] != length - 1:
        positions.append(length - 1)
    return np.asarray(positions, dtype=int)


def _trimmed_median(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    if finite.size < 16 or np.nanmin(finite) == np.nanmax(finite):
        return float(np.median(finite))
    counts, edges = np.histogram(finite, bins="fd")
    modal_index = int(np.argmax(counts))
    mode = 0.5 * (edges[modal_index] + edges[modal_index + 1])
    trimmed = finite[finite <= 2.0 * mode] if mode > 0 else finite
    return float(np.median(trimmed if trimmed.size else finite))


def _sample_local_background(
    data: np.ndarray,
    *,
    filter_size_pixels: int,
    sampling_step_pixels: int,
    excluded_mask: np.ndarray | None = None,
    centers_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = data.shape
    half = max(1, int(round(filter_size_pixels / 2)))
    ys = _grid_positions(ny, sampling_step_pixels)
    xs = _grid_positions(nx, sampling_step_pixels)
    points: list[tuple[float, float]] = []
    values: list[float] = []
    for row in ys:
        for col in xs:
            if centers_mask is not None and not centers_mask[row, col]:
                continue
            y0, y1 = max(0, row - half), min(ny, row + half + 1)
            x0, x1 = max(0, col - half), min(nx, col + half + 1)
            window = data[y0:y1, x0:x1]
            if excluded_mask is not None:
                keep = ~excluded_mask[y0:y1, x0:x1]
                window = window[keep]
            value = _trimmed_median(window)
            if np.isfinite(value):
                points.append((float(col), float(row)))
                values.append(value)
    if len(points) < 3:
        raise ValueError("Background estimation produced fewer than three samples.")
    return np.asarray(points, dtype=float), np.asarray(values, dtype=float)


def _interpolate_sample_grid(
    points: np.ndarray,
    values: np.ndarray,
    shape: tuple[int, int],
) -> np.ndarray:
    yy, xx = np.indices(shape, dtype=float)
    linear = griddata(points, values, (xx, yy), method="linear")
    if np.any(~np.isfinite(linear)):
        nearest = griddata(points, values, (xx, yy), method="nearest")
        linear = np.where(np.isfinite(linear), linear, nearest)
    return np.asarray(linear, dtype=float)


def estimate_lmf_background(
    image: np.ndarray,
    wcs: WCS,
    *,
    filter_size_arcmin: float = 13.0,
    sampling_arcsec: float = 24.0,
) -> BackgroundResult:
    """Estimate a BT09-style large-scale median-filter background."""

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be two-dimensional.")
    if filter_size_arcmin <= 0 or sampling_arcsec <= 0:
        raise ValueError("Filter size and sampling must be positive.")
    scale = pixel_scale_arcsec(wcs)
    filter_pixels = max(3, int(round(filter_size_arcmin * 60.0 / scale)))
    step_pixels = max(1, int(round(sampling_arcsec / scale)))
    points, values = _sample_local_background(
        data,
        filter_size_pixels=filter_pixels,
        sampling_step_pixels=step_pixels,
    )
    background = _interpolate_sample_grid(points, values, data.shape)
    return BackgroundResult(
        values=background,
        method="lmf",
        diagnostics={
            "filter_size_arcmin": float(filter_size_arcmin),
            "filter_size_pixels": filter_pixels,
            "sampling_arcsec": float(sampling_arcsec),
            "sample_count": len(values),
        },
    )


def measure_box_background(
    image: np.ndarray,
    boxes: list[tuple[int, int, int, int]]
    | tuple[tuple[int, int, int, int], ...],
    *,
    maximum_intensity: float | None = 15.0,
    minimum_intensity: float | None = None,
) -> tuple[float, dict[str, object]]:
    """Measure an observed-background level from adjacent pixel boxes.

    Each box is ``(row_start, row_stop, col_start, col_stop)`` with a
    stop-exclusive convention matching NumPy slices. The median of each
    valid box is measured, then those medians are averaged. The default
    15 MJy/sr ceiling follows the supplied Sgr C F480M analysis and removes
    bright stellar emission. Unlike :func:`estimate_box_background`, this
    function does not allocate a full-size background image, making it useful
    when the boxes are measured on a large mosaic before mapping a cutout.
    """

    data = np.asarray(image)
    if data.ndim != 2:
        raise ValueError("image must be two-dimensional.")
    if not boxes:
        raise ValueError("At least one adjacent background box is required.")
    if (
        maximum_intensity is not None
        and minimum_intensity is not None
        and maximum_intensity <= minimum_intensity
    ):
        raise ValueError("maximum_intensity must exceed minimum_intensity.")

    ny, nx = data.shape
    medians: list[float] = []
    counts: list[int] = []
    normalized_boxes: list[tuple[int, int, int, int]] = []
    for box in boxes:
        if len(box) != 4:
            raise ValueError(
                "Each box must contain row_start, row_stop, col_start, col_stop."
            )
        row_start, row_stop, col_start, col_stop = map(int, box)
        if not (
            0 <= row_start < row_stop <= ny
            and 0 <= col_start < col_stop <= nx
        ):
            raise ValueError(f"Background box {box!r} falls outside image bounds.")
        values = data[row_start:row_stop, col_start:col_stop]
        valid = np.isfinite(values)
        if maximum_intensity is not None:
            valid &= values < float(maximum_intensity)
        if minimum_intensity is not None:
            valid &= values > float(minimum_intensity)
        selected = values[valid]
        if selected.size == 0:
            raise ValueError(
                f"Background box {box!r} has no pixels after intensity cuts."
            )
        medians.append(float(np.median(selected)))
        counts.append(int(selected.size))
        normalized_boxes.append((row_start, row_stop, col_start, col_stop))

    level = float(np.mean(medians))
    scatter = float(np.std(medians, ddof=1)) if len(medians) > 1 else 0.0
    diagnostics: dict[str, object] = {
        "boxes": normalized_boxes,
        "box_medians": medians,
        "box_valid_pixel_counts": counts,
        "background_level": level,
        "box_median_scatter": scatter,
        "box_median_standard_error": (
            scatter / np.sqrt(len(medians)) if medians else np.nan
        ),
        "maximum_intensity": maximum_intensity,
        "minimum_intensity": minimum_intensity,
    }
    return level, diagnostics


def estimate_box_background(
    image: np.ndarray,
    boxes: list[tuple[int, int, int, int]]
    | tuple[tuple[int, int, int, int], ...],
    *,
    maximum_intensity: float | None = 15.0,
    minimum_intensity: float | None = None,
) -> BackgroundResult:
    """Estimate a constant background image from adjacent pixel boxes."""

    data = np.asarray(image, dtype=float)
    level, diagnostics = measure_box_background(
        data,
        boxes,
        maximum_intensity=maximum_intensity,
        minimum_intensity=minimum_intensity,
    )
    return BackgroundResult(
        values=np.full(data.shape, level, dtype=float),
        method="boxes",
        diagnostics=diagnostics,
    )


def _inverse_square_fill(
    baseline: np.ndarray,
    outside_mask: np.ndarray,
    *,
    radius_pixels: int,
) -> np.ndarray:
    coordinates = np.arange(-radius_pixels, radius_pixels + 1, dtype=float)
    dy, dx = np.meshgrid(coordinates, coordinates, indexing="ij")
    radius2 = dx**2 + dy**2
    kernel = np.zeros_like(radius2)
    eligible = (radius2 > 0) & (radius2 <= radius_pixels**2)
    kernel[eligible] = 1.0 / radius2[eligible]

    valid = outside_mask & np.isfinite(baseline)
    numerator = fftconvolve(
        np.where(valid, baseline, 0.0),
        kernel,
        mode="same",
    )
    denominator = fftconvolve(valid.astype(float), kernel, mode="same")
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(baseline, np.nan, dtype=float),
        where=denominator > 0,
    )


def estimate_smf_background(
    image: np.ndarray,
    wcs: WCS,
    cloud: CloudEllipse,
    *,
    filter_fraction: float = 1.0 / 3.0,
    sampling_arcsec: float = 24.0,
) -> BackgroundResult:
    """Estimate a BT09-style small-scale median-filter background.

    Outside the Simon ellipse, local trimmed medians are sampled on the
    published 24-arcsecond grid. Inside the ellipse, those exterior values
    are interpolated with inverse-square separation weighting out to one
    semi-major-axis radius.
    """

    data = np.asarray(image, dtype=float)
    if data.ndim != 2:
        raise ValueError("image must be two-dimensional.")
    if filter_fraction <= 0 or sampling_arcsec <= 0:
        raise ValueError("Filter fraction and sampling must be positive.")
    scale = pixel_scale_arcsec(wcs)
    cloud_mask = ellipse_mask(data.shape, wcs, cloud)
    outside = ~cloud_mask
    filter_arcmin = cloud.major_axis_arcmin * filter_fraction
    filter_pixels = max(3, int(round(filter_arcmin * 60.0 / scale)))
    step_pixels = max(1, int(round(sampling_arcsec / scale)))

    points, values = _sample_local_background(
        data,
        filter_size_pixels=filter_pixels,
        sampling_step_pixels=step_pixels,
        excluded_mask=cloud_mask,
        centers_mask=outside,
    )
    baseline = _interpolate_sample_grid(points, values, data.shape)
    radius_pixels = max(
        1,
        int(round((cloud.major_axis_arcmin / 2.0) * 60.0 / scale)),
    )
    inside_interpolation = _inverse_square_fill(
        baseline,
        outside,
        radius_pixels=radius_pixels,
    )
    background = np.where(cloud_mask, inside_interpolation, baseline)
    return BackgroundResult(
        values=background,
        method="smf",
        diagnostics={
            "filter_fraction": float(filter_fraction),
            "filter_size_arcmin": float(filter_arcmin),
            "filter_size_pixels": filter_pixels,
            "sampling_arcsec": float(sampling_arcsec),
            "interpolation_radius_pixels": radius_pixels,
            "sample_count": len(values),
        },
    )
