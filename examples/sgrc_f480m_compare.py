"""Compare BT12 and GTL foregrounds in Rubén Fedriani's Sgr C box.

The five touching background boxes used here are a reproducible comparison
choice, not a claim that they are Rubén's unpublished final background boxes.
Replace ``touching_background_boxes`` when those machine-readable regions are
available.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation

from gtlmapping import (
    G_PER_CM2_TO_MSUN_PER_PC2,
    GTLMapper,
    compute_extinction,
    measure_box_background,
)
from gtlmapping.geometry import pixel_scale_arcsec


def touching_background_boxes(
    x_range: tuple[int, int],
    y_range: tuple[int, int],
) -> list[tuple[int, int, int, int]]:
    """Return the five in-bounds neighbors used for the edge-adjacent target."""

    x0, x1 = x_range
    y0, y1 = y_range
    width = x1 - x0
    height = y1 - y0
    return [
        (y0 - height, y0, x0, x1),  # north
        (y0 - height, y0, x1, x1 + width),  # northeast
        (y0, y1, x1, x1 + width),  # east
        (y1, y1 + height, x1, x1 + width),  # southeast
        (y1, y1 + height, x0, x1),  # south
    ]


def integrated_mass_msun(
    surface_density: np.ma.MaskedArray,
    *,
    pixel_scale: float,
    distance_kpc: float,
) -> float:
    """Integrate surface density over all unmasked pixels in a box."""

    pixel_pc = distance_kpc * 1000.0 * np.deg2rad(pixel_scale / 3600.0)
    return float(
        np.sum(surface_density.filled(0.0))
        * G_PER_CM2_TO_MSUN_PER_PC2
        * pixel_pc**2
    )


def run_comparison(
    input_path: Path,
    output_dir: Path,
    *,
    x_range: tuple[int, int] = (77, 200),
    y_range: tuple[int, int] = (1980, 2103),
    distance_kpc: float = 8.15,
) -> dict[str, object]:
    """Run the controlled BT12/GTL comparison and write FITS/PNG products."""

    output_dir.mkdir(parents=True, exist_ok=True)
    section = np.s_[y_range[0] : y_range[1], x_range[0] : x_range[1]]
    with fits.open(input_path, memmap=True) as hdul:
        full_science = hdul["SCI"].data
        observed = np.asarray(full_science[section], dtype=float)
        observed_std = np.asarray(hdul["ERR"].data[section], dtype=float)
        cutout_wcs = WCS(hdul["SCI"].header).slice(section)
        header = cutout_wcs.to_header()
        boxes = touching_background_boxes(x_range, y_range)
        background_level, background_diagnostics = measure_box_background(
            full_science,
            boxes,
            maximum_intensity=15.0,
        )

    noise_sigma = float(np.nanmedian(observed_std))
    background_std = float(
        background_diagnostics["box_median_standard_error"]
    )
    background = np.full(observed.shape, background_level, dtype=float)
    region_mask = np.ones(observed.shape, dtype=bool)

    bt12 = GTLMapper(
        observed,
        header=header,
        wcs=cutout_wcs,
        observed_std=observed_std,
    )
    bt12.set_background(background, uncertainty=background_std)
    bt12_foreground = bt12.fit_foreground(
        method="bt12",
        region_mask=region_mask,
        noise_sigma=noise_sigma,
        min_separation_arcsec=0.74,
    )
    bt12_result = bt12.compute(
        filter_name="F480M",
        gas_to_dust_ratio=156.0,
        kappa_std_cm2_g=0.30 * 9.76,
        bright_pixel_policy="zero",
    )

    gtl = GTLMapper(
        observed,
        header=header,
        wcs=cutout_wcs,
        observed_std=observed_std,
    )
    gtl.set_background(background, uncertainty=background_std)
    samples = gtl.detect_foreground(
        grid_n=8,
        cover_edges=True,
        noise_sigma=noise_sigma,
        beam_fwhm_arcsec=0.18,
        min_separation_arcsec=0.74,
    )
    gtl_foreground = gtl.fit_foreground(
        method="conservative",
        region_mask=region_mask,
        noise_sigma=noise_sigma,
        min_separation_arcsec=0.74,
    )
    gtl_detected = gtl.compute(
        filter_name="F480M",
        gas_to_dust_ratio=156.0,
        kappa_std_cm2_g=0.30 * 9.76,
        bright_pixel_policy="zero",
    )
    gtl_lower_limits = gtl.compute(
        filter_name="F480M",
        gas_to_dust_ratio=156.0,
        kappa_std_cm2_g=0.30 * 9.76,
        bright_pixel_policy="zero",
        saturation_policy="lower_limit",
        intensity_floor=2.0 * noise_sigma,
    )

    minimum = float(bt12_foreground.diagnostics["global_minimum"])
    bt12_saturated = (
        (observed > minimum) & (observed < minimum + 2.0 * noise_sigma)
    )
    bt12_local_saturated = observed <= (
        bt12_foreground.values + 2.0 * noise_sigma
    )
    gtl_local_saturated = observed <= (
        gtl_foreground.values + 2.0 * noise_sigma
    )
    scale = pixel_scale_arcsec(cutout_wcs)
    bt12_mass = integrated_mass_msun(
        bt12_result.surface_density,
        pixel_scale=scale,
        distance_kpc=distance_kpc,
    )
    gtl_detected_mass = integrated_mass_msun(
        gtl_detected.surface_density,
        pixel_scale=scale,
        distance_kpc=distance_kpc,
    )
    gtl_lower_limit_mass = integrated_mass_msun(
        gtl_lower_limits.surface_density,
        pixel_scale=scale,
        distance_kpc=distance_kpc,
    )
    grid_sensitivity: list[dict[str, float | int]] = []
    for grid_n in range(4, 11):
        trial = GTLMapper(
            observed,
            header=header,
            wcs=cutout_wcs,
            observed_std=observed_std,
        )
        trial.set_background(background, uncertainty=background_std)
        trial_samples = trial.detect_foreground(
            grid_n=grid_n,
            cover_edges=True,
            noise_sigma=noise_sigma,
            beam_fwhm_arcsec=0.18,
            min_separation_arcsec=0.74,
        )
        trial_foreground = trial.fit_foreground(
            method="conservative",
            region_mask=region_mask,
            noise_sigma=noise_sigma,
            min_separation_arcsec=0.74,
        )
        trial_result = trial.compute(
            filter_name="F480M",
            gas_to_dust_ratio=156.0,
            bright_pixel_policy="zero",
            saturation_policy="lower_limit",
            intensity_floor=2.0 * noise_sigma,
        )
        trial_mass = integrated_mass_msun(
            trial_result.surface_density,
            pixel_scale=scale,
            distance_kpc=distance_kpc,
        )
        grid_sensitivity.append(
            {
                "grid_n": grid_n,
                "raw_sample_count": trial_samples.raw_detection_count,
                "unique_sample_count": len(trial_samples),
                "local_saturated_pixels": int(
                    np.count_nonzero(
                        observed <= trial_foreground.values + 2.0 * noise_sigma
                    )
                ),
                "strict_lower_limit_pixels": int(
                    np.count_nonzero(trial_result.saturated_mask)
                ),
                "mass_msun": trial_mass,
                "gtl_to_bt12_mass_ratio": trial_mass / bt12_mass,
            }
        )

    source_mask_sensitivity: dict[str, dict[str, float | int]] = {}
    for radius in (0, 3, 6):
        if radius == 0:
            source_mask = bt12_result.bright_mask
        else:
            yy, xx = np.indices((2 * radius + 1, 2 * radius + 1)) - radius
            disk = xx**2 + yy**2 <= radius**2
            source_mask = binary_dilation(
                bt12_result.bright_mask,
                structure=disk,
            )
        keep = ~source_mask
        bt12_masked_mass = integrated_mass_msun(
            np.ma.array(
                bt12_result.surface_density.data,
                mask=np.ma.getmaskarray(bt12_result.surface_density) | source_mask,
            ),
            pixel_scale=scale,
            distance_kpc=distance_kpc,
        )
        gtl_masked_mass = integrated_mass_msun(
            np.ma.array(
                gtl_lower_limits.surface_density.data,
                mask=(
                    np.ma.getmaskarray(gtl_lower_limits.surface_density)
                    | source_mask
                ),
            ),
            pixel_scale=scale,
            distance_kpc=distance_kpc,
        )
        source_mask_sensitivity[str(radius)] = {
            "masked_pixel_count": int(np.count_nonzero(source_mask)),
            "kept_pixel_fraction": float(np.mean(keep)),
            "bt12_mass_msun": bt12_masked_mass,
            "gtl_mass_msun": gtl_masked_mass,
            "gtl_to_bt12_mass_ratio": gtl_masked_mass / bt12_masked_mass,
        }

    background_sensitivity: list[dict[str, float | int]] = []
    background_scatter = float(
        background_diagnostics["box_median_scatter"]
    )
    for level in (
        background_level - background_scatter,
        background_level,
        background_level + background_scatter,
    ):
        level_map = np.full(observed.shape, level)
        _, bt12_sigma, _, _, _ = compute_extinction(
            observed,
            level_map,
            bt12_foreground.values,
            kappa_cm2_g=9.76,
            bright_pixel_policy="zero",
        )
        _, gtl_sigma, _, invalid_background, _ = compute_extinction(
            observed,
            level_map,
            gtl_foreground.values,
            kappa_cm2_g=9.76,
            bright_pixel_policy="zero",
            saturation_policy="lower_limit",
            intensity_floor=2.0 * noise_sigma,
        )
        bt12_sensitivity_mass = integrated_mass_msun(
            bt12_sigma,
            pixel_scale=scale,
            distance_kpc=distance_kpc,
        )
        gtl_sensitivity_mass = integrated_mass_msun(
            gtl_sigma,
            pixel_scale=scale,
            distance_kpc=distance_kpc,
        )
        background_sensitivity.append(
            {
                "background_level_mjy_sr": float(level),
                "bt12_mass_msun": bt12_sensitivity_mass,
                "gtl_mass_msun": gtl_sensitivity_mass,
                "gtl_to_bt12_mass_ratio": (
                    gtl_sensitivity_mass / bt12_sensitivity_mass
                ),
                "gtl_invalid_background_pixels": int(
                    np.count_nonzero(invalid_background)
                ),
            }
        )

    shared_metadata = {
        "target": "Sgr C E NIRCam",
        "x_range": str(x_range),
        "y_range": str(y_range),
        "bg_level": background_level,
        "bg_boxes": "N,NE,E,SE,S touching target",
        "distance_kpc": distance_kpc,
    }
    bt12_result.metadata.update(shared_metadata)
    gtl_lower_limits.metadata.update(shared_metadata)
    bt12_result.write(output_dir / "sgrc_f480m_bt12.fits", overwrite=True)
    gtl_lower_limits.write(
        output_dir / "sgrc_f480m_gtl_lower_limits.fits",
        overwrite=True,
    )

    summary: dict[str, object] = {
        "input": str(input_path),
        "x_range": list(x_range),
        "y_range": list(y_range),
        "distance_kpc": distance_kpc,
        "noise_sigma_mjy_sr": noise_sigma,
        "background_level_mjy_sr": background_level,
        "background_boxes_row_start_row_stop_col_start_col_stop": [
            list(box) for box in boxes
        ],
        "background_box_medians_mjy_sr": background_diagnostics[
            "box_medians"
        ],
        "background_box_scatter_mjy_sr": background_diagnostics[
            "box_median_scatter"
        ],
        "background_box_standard_error_mjy_sr": background_std,
        "bt12_foreground_mjy_sr": float(bt12_foreground.values[0, 0]),
        "gtl_foreground_min_median_max_mjy_sr": [
            float(np.min(gtl_foreground.values)),
            float(np.median(gtl_foreground.values)),
            float(np.max(gtl_foreground.values)),
        ],
        "bt12_saturated_pixels": int(np.count_nonzero(bt12_saturated)),
        "bt12_local_saturated_pixels": int(
            np.count_nonzero(bt12_local_saturated)
        ),
        "bt12_independent_saturated_pixels": int(
            bt12_foreground.diagnostics[
                "independent_saturated_pixel_count"
            ]
        ),
        "gtl_local_saturated_pixels": int(
            np.count_nonzero(gtl_local_saturated)
        ),
        "gtl_strict_lower_limit_pixels": int(
            np.count_nonzero(gtl_lower_limits.saturated_mask)
        ),
        "gtl_raw_sample_count": samples.raw_detection_count,
        "gtl_unique_sample_count": len(samples),
        "conservative_blend_factor": float(
            gtl_foreground.diagnostics["blend_factor"]
        ),
        "conservative_unregularized_local_saturation_pixels": int(
            gtl_foreground.diagnostics[
                "unregularized_local_saturation_count"
            ]
        ),
        "conservative_unregularized_strict_saturation_pixels": int(
            gtl_foreground.diagnostics[
                "unregularized_strict_saturation_count"
            ]
        ),
        "bt12_mass_msun": bt12_mass,
        "gtl_detected_only_mass_msun": gtl_detected_mass,
        "gtl_lower_limit_mass_msun": gtl_lower_limit_mass,
        "gtl_to_bt12_detected_mass_ratio": gtl_detected_mass / bt12_mass,
        "gtl_to_bt12_lower_limit_mass_ratio": (
            gtl_lower_limit_mass / bt12_mass
        ),
        "grid_sensitivity": grid_sensitivity,
        "opacity_systematic_fraction": 0.30,
        "bt12_opacity_systematic_msun": 0.30 * bt12_mass,
        "gtl_opacity_systematic_msun": 0.30 * gtl_lower_limit_mass,
        "source_mask_dilation_sensitivity": source_mask_sensitivity,
        "background_sensitivity": background_sensitivity,
        "background_box_note": (
            "Five touching comparison boxes; replace with Ruben's exact "
            "machine-readable background regions when available."
        ),
    }
    (output_dir / "sgrc_f480m_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return summary

    bt12_sigma = bt12_result.surface_density.filled(np.nan)
    gtl_sigma = gtl_lower_limits.surface_density.filled(np.nan)
    common_vmax = float(np.nanpercentile(gtl_sigma, 99.0))
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(12.0, 7.5),
        constrained_layout=True,
    )
    panels = [
        (observed, "F480M intensity", "viridis", 2.0, 15.0),
        (
            bt12_sigma,
            f"BT12 Sigma ({bt12_mass:.1f} Msun)",
            "magma",
            0.0,
            common_vmax,
        ),
        (
            gtl_sigma,
            f"Conservative GTL Sigma ({gtl_lower_limit_mass:.1f} Msun)",
            "magma",
            0.0,
            common_vmax,
        ),
        (
            bt12_local_saturated,
            f"BT12 within 2 sigma: {bt12_local_saturated.sum()}",
            "gray_r",
            0,
            1,
        ),
        (
            gtl_local_saturated,
            f"GTL locally consistent: {gtl_local_saturated.sum()}",
            "gray_r",
            0,
            1,
        ),
        (
            gtl_sigma - bt12_sigma,
            "GTL - BT12 Sigma",
            "Reds",
            0.0,
            common_vmax / 2.0,
        ),
    ]
    for axis, (values, title, cmap, vmin, vmax) in zip(
        axes.flat,
        panels,
        strict=True,
    ):
        image = axis.imshow(
            values,
            origin="lower",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        axis.set_title(title)
        axis.set_xlabel("x pixel in Rubén box")
        axis.set_ylabel("y pixel in Rubén box")
        figure.colorbar(image, ax=axis, fraction=0.046)
    figure.suptitle(
        "Sgr C F480M: BT12 vs BT12-floored spatial trend",
        fontsize=14,
    )
    figure.savefig(output_dir / "sgrc_f480m_comparison.png", dpi=200)
    plt.close(figure)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation/sgrc_f480m"),
    )
    parser.add_argument("--distance-kpc", type=float, default=8.15)
    arguments = parser.parse_args()
    summary = run_comparison(
        arguments.input,
        arguments.output_dir,
        distance_kpc=arguments.distance_kpc,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
