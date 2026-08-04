"""Create the Cloud C, F, and H comparison figures used in the documentation."""

from __future__ import annotations

import argparse
import gc
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from gtlmapping import GTLMapper
from gtlmapping.geometry import ellipse_mask


@dataclass(frozen=True)
class CloudCase:
    label: str
    catalog_name: str
    image_path: Path
    background_path: Path | None = None


def _crop_bounds(mask: np.ndarray) -> tuple[slice, slice]:
    rows, cols = np.where(mask)
    if rows.size == 0:
        raise ValueError("The catalog ellipse does not overlap the image.")
    height = int(rows.max() - rows.min() + 1)
    width = int(cols.max() - cols.min() + 1)
    pad_y = max(10, int(round(0.04 * height)))
    pad_x = max(10, int(round(0.04 * width)))
    y0 = max(0, int(rows.min()) - pad_y)
    y1 = min(mask.shape[0], int(rows.max()) + pad_y + 1)
    x0 = max(0, int(cols.min()) - pad_x)
    x1 = min(mask.shape[1], int(cols.max()) + pad_x + 1)
    return slice(y0, y1), slice(x0, x1)


def _profile_result(
    base: GTLMapper,
    *,
    method: str,
    cloud,
    samples,
    background: np.ndarray,
    noise_sigma: float,
):
    mapper = GTLMapper(
        base.observed,
        header=base.header,
        wcs=base.wcs,
        observed_std=base.observed_std,
    )
    mapper.cloud = cloud
    mapper.samples = samples
    mapper.set_background(background, method="shared_smf")
    mapper.fit_foreground(
        method=method,
        samples=samples,
        cloud=cloud,
        noise_sigma=noise_sigma,
        min_separation_arcsec=8.0,
    )
    if method == "moderate":
        return mapper.compute_moderate(
            kappa_cm2_g=7.5,
            bright_pixel_policy="zero",
        )
    if method == "liberal":
        return mapper.compute_liberal(
            kappa_cm2_g=7.5,
            bright_pixel_policy="zero",
        )
    return mapper.compute(
        kappa_cm2_g=7.5,
        bright_pixel_policy="zero",
    )


def _positive_sum(values: np.ma.MaskedArray, mask: np.ndarray) -> float:
    data = values.filled(np.nan)
    selected = data[mask & np.isfinite(data)]
    return float(np.sum(np.maximum(selected, 0.0)))


def _plot_case(
    case: CloudCase,
    *,
    catalog_path: Path,
    output_dir: Path,
    noise_sigma: float,
) -> dict[str, object]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import AsinhNorm, Normalize, TwoSlopeNorm
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "Install GTLMapping with the plot extra to create these figures."
        ) from exc

    print(f"Cloud {case.label}: loading {case.image_path.name}", flush=True)
    base = GTLMapper.from_fits(case.image_path)
    cloud = base.select_cloud(catalog_path, case.catalog_name)
    mask = ellipse_mask(base.observed.shape, base.wcs, cloud)
    samples = base.detect_foreground(noise_sigma=noise_sigma)
    if case.background_path is None:
        print(f"Cloud {case.label}: estimating the SMF background", flush=True)
        background_result = base.estimate_background(method="smf", cloud=cloud)
    else:
        print(f"Cloud {case.label}: loading the prepared background", flush=True)
        background_result = base.set_background_from_fits(case.background_path)
    background = background_result.values
    crop = _crop_bounds(mask)

    profile_names = ("bt12", "conservative", "moderate", "liberal")
    metrics: dict[str, dict[str, float | int]] = {}
    bt12_sigma = None
    bt12_foreground = None
    moderate_sigma = None
    moderate_foreground = None
    moderate_saturated = None

    for method in profile_names:
        print(f"Cloud {case.label}: calculating {method}", flush=True)
        result = _profile_result(
            base,
            method=method,
            cloud=cloud,
            samples=samples,
            background=background,
            noise_sigma=noise_sigma,
        )
        foreground = np.asarray(result.foreground, dtype=float)
        strict = np.asarray(result.saturated_mask, dtype=bool) & mask
        local = mask & np.isfinite(base.observed) & (
            base.observed <= foreground + 2.0 * noise_sigma
        )
        metrics[method] = {
            "surface_density_sum": _positive_sum(result.surface_density, mask),
            "local_saturation_pixels": int(np.count_nonzero(local)),
            "strict_lower_limit_pixels": int(np.count_nonzero(strict)),
            "foreground_min": float(np.nanmin(foreground[mask])),
            "foreground_median": float(np.nanmedian(foreground[mask])),
            "foreground_max": float(np.nanmax(foreground[mask])),
        }
        if method == "bt12":
            bt12_sigma = result.surface_density[crop].filled(np.nan).copy()
            bt12_foreground = foreground[crop].copy()
        elif method == "moderate":
            moderate_sigma = result.surface_density[crop].filled(np.nan).copy()
            moderate_foreground = foreground[crop].copy()
            moderate_saturated = strict[crop].copy()
        del result
        gc.collect()

    if any(
        item is None
        for item in (
            bt12_sigma,
            bt12_foreground,
            moderate_sigma,
            moderate_foreground,
            moderate_saturated,
        )
    ):
        raise RuntimeError("The BT12 or moderate profile was not calculated.")

    bt12_sum = metrics["bt12"]["surface_density_sum"]
    for method in profile_names:
        profile_sum = metrics[method]["surface_density_sum"]
        metrics[method]["ratio_to_bt12"] = float(profile_sum / bt12_sum)

    observed_crop = base.observed[crop]
    mask_crop = mask[crop]
    y0 = 0 if crop[0].start is None else crop[0].start
    x0 = 0 if crop[1].start is None else crop[1].start
    sample_rows = samples.rows - y0
    sample_cols = samples.cols - x0
    inside_crop = (
        (sample_rows >= 0)
        & (sample_rows < observed_crop.shape[0])
        & (sample_cols >= 0)
        & (sample_cols < observed_crop.shape[1])
    )

    finite_observed = observed_crop[np.isfinite(observed_crop)]
    obs_low, obs_high = np.nanpercentile(finite_observed, (1.0, 99.5))
    observed_norm = Normalize(vmin=float(obs_low), vmax=float(obs_high))

    foreground_change = np.where(
        mask_crop,
        moderate_foreground - bt12_foreground,
        np.nan,
    )
    bt12_sigma = np.where(mask_crop, bt12_sigma, np.nan)
    moderate_sigma = np.where(mask_crop, moderate_sigma, np.nan)
    finite_change = np.abs(foreground_change[np.isfinite(foreground_change)])
    change_limit = max(float(np.nanpercentile(finite_change, 99.0)), 1.0e-6)
    change_norm = TwoSlopeNorm(
        vmin=-change_limit,
        vcenter=0.0,
        vmax=change_limit,
    )

    density_values = np.concatenate(
        (
            bt12_sigma[mask_crop & np.isfinite(bt12_sigma)],
            moderate_sigma[mask_crop & np.isfinite(moderate_sigma)],
        )
    )
    density_high = max(float(np.nanpercentile(density_values, 99.5)), 0.05)
    density_norm = AsinhNorm(
        linear_width=max(0.01 * density_high, 1.0e-4),
        vmin=0.0,
        vmax=density_high,
    )

    panel_aspect = observed_crop.shape[1] / observed_crop.shape[0]
    figure_height = max(5.6, min(9.0, 1.0 + 10.0 / panel_aspect))
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(12.0, figure_height),
        constrained_layout=True,
    )
    observed_image = axes[0, 0].imshow(
        observed_crop,
        origin="lower",
        cmap="gray",
        norm=observed_norm,
        interpolation="nearest",
    )
    axes[0, 0].contour(mask_crop, levels=[0.5], colors="white", linewidths=1.0)
    axes[0, 0].scatter(
        sample_cols[inside_crop],
        sample_rows[inside_crop],
        s=24,
        facecolors="none",
        edgecolors="#22d3ee",
        linewidths=1.2,
        label=(
            f"{np.count_nonzero(inside_crop)} shown of "
            f"{len(samples)} accepted sites"
        ),
    )
    axes[0, 0].legend(loc="upper right", framealpha=0.85, fontsize=9)
    axes[0, 0].set_title(f"Cloud {case.label}: observed 8 micron intensity")
    figure.colorbar(observed_image, ax=axes[0, 0], label="MJy sr$^{-1}$")

    change_image = axes[0, 1].imshow(
        foreground_change,
        origin="lower",
        cmap="coolwarm",
        norm=change_norm,
        interpolation="nearest",
    )
    axes[0, 1].contour(mask_crop, levels=[0.5], colors="black", linewidths=0.8)
    axes[0, 1].set_title("Moderate foreground minus BT12")
    figure.colorbar(change_image, ax=axes[0, 1], label="MJy sr$^{-1}$")

    bt12_image = axes[1, 0].imshow(
        bt12_sigma,
        origin="lower",
        cmap="magma",
        norm=density_norm,
        interpolation="nearest",
    )
    axes[1, 0].contour(mask_crop, levels=[0.5], colors="white", linewidths=0.8)
    axes[1, 0].set_title("BT12 surface density")
    figure.colorbar(bt12_image, ax=axes[1, 0], label="$\\Sigma$ (g cm$^{-2}$)")

    moderate_image = axes[1, 1].imshow(
        moderate_sigma,
        origin="lower",
        cmap="magma",
        norm=density_norm,
        interpolation="nearest",
    )
    axes[1, 1].contour(mask_crop, levels=[0.5], colors="white", linewidths=0.8)
    saturated_rows, saturated_cols = np.where(moderate_saturated)
    if saturated_rows.size:
        axes[1, 1].scatter(
            saturated_cols,
            saturated_rows,
            s=26,
            facecolors="none",
            edgecolors="#f0abfc",
            linewidths=1.2,
            label="Censored lower limit",
        )
        axes[1, 1].legend(loc="upper right", framealpha=0.85, fontsize=9)
    axes[1, 1].set_title("Moderate GTL surface density")
    figure.colorbar(moderate_image, ax=axes[1, 1], label="$\\Sigma$ (g cm$^{-2}$)")

    for axis in axes.flat:
        axis.set_xlabel("Image column")
        axis.set_ylabel("Image row")

    output_path = output_dir / f"cloud_{case.label.lower()}_method.png"
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    print(f"Cloud {case.label}: wrote {output_path}", flush=True)

    return {
        "cloud": case.label,
        "catalog_name": case.catalog_name,
        "image": case.image_path.name,
        "background": (
            case.background_path.name
            if case.background_path is not None
            else "GTLMapping SMF estimate"
        ),
        "noise_sigma_mjy_sr": noise_sigma,
        "accepted_windows": samples.accepted_windows,
        "accepted_unique_sites": len(samples),
        "profiles": metrics,
        "figure": output_path.name,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--cloud-c", required=True, type=Path)
    parser.add_argument("--cloud-c-background", required=True, type=Path)
    parser.add_argument("--cloud-f", required=True, type=Path)
    parser.add_argument("--cloud-h", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--noise-sigma", type=float, default=0.6)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cases = (
        CloudCase(
            "C",
            "G028.37+00.07",
            args.cloud_c,
            args.cloud_c_background,
        ),
        CloudCase("F", "G034.43+00.24", args.cloud_f),
        CloudCase("H", "G035.39-00.33", args.cloud_h),
    )
    summary = [
        _plot_case(
            case,
            catalog_path=args.catalog,
            output_dir=args.output_dir,
            noise_sigma=args.noise_sigma,
        )
        for case in cases
    ]
    (args.output_dir / "cloud_method_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
