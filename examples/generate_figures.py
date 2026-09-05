"""Reproduce the figures and machine-readable tables for the GTLMapping paper.

The preparation step is deliberately one cloud and one profile at a time.
Cloud F is large enough that isolating profiles keeps peak memory use modest.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits

from gtlmapping import GTLMapper
from gtlmapping.geometry import ellipse_mask, pixel_scale_arcsec


TEAL = "#1FA89C"
TERRACOTTA = "#C65A2E"
PLUM = "#6A4C93"
STEEL = "#3D6B9E"
SAFFRON = "#E0A23C"
OFF_WHITE = "#F7F3EB"
PROFILE_COLORS = {
    "conservative": TEAL,
    "moderate": SAFFRON,
    "liberal": TERRACOTTA,
}
PROFILE_MARKERS = {
    "conservative": "o",
    "moderate": "s",
    "liberal": "^",
}
PDF_PROFILE_COLORS = {
    "bt12": "#2B2730",
    "conservative": TEAL,
    "moderate": SAFFRON,
    "liberal": TERRACOTTA,
}
PDF_PROFILE_LINESTYLES = {
    "bt12": "-",
    "conservative": (0, (5, 2)),
    "moderate": (0, (3, 1.5)),
    "liberal": (0, (1, 1.4)),
}


@dataclass(frozen=True)
class CloudCase:
    label: str
    catalog_name: str
    image_path: Path
    background_path: Path | None


def crop_bounds(mask: np.ndarray) -> tuple[slice, slice]:
    rows, cols = np.where(mask)
    if rows.size == 0:
        raise ValueError("The catalog ellipse does not overlap the image.")
    height = int(rows.max() - rows.min() + 1)
    width = int(cols.max() - cols.min() + 1)
    pad_y = max(10, int(round(0.04 * height)))
    pad_x = max(10, int(round(0.04 * width)))
    return (
        slice(max(0, int(rows.min()) - pad_y), min(mask.shape[0], int(rows.max()) + pad_y + 1)),
        slice(max(0, int(cols.min()) - pad_x), min(mask.shape[1], int(cols.max()) + pad_x + 1)),
    )


def positive_sum(values: np.ma.MaskedArray, mask: np.ndarray) -> float:
    data = values.filled(np.nan)
    selected = data[mask & np.isfinite(data)]
    return float(np.sum(np.maximum(selected, 0.0)))


def prepare_cloud(
    case: CloudCase,
    method: str,
    *,
    catalog_path: Path,
    output_dir: Path,
    noise_sigma: float,
) -> Path:
    if method not in {"bt12", "conservative", "moderate", "liberal"}:
        raise ValueError(
            "Paper map products are prepared for bt12, conservative, "
            "moderate, or liberal."
        )

    print(f"Loading Cloud {case.label} ({case.image_path.name})", flush=True)
    base = GTLMapper.from_fits(case.image_path)
    cloud = base.select_cloud(catalog_path, case.catalog_name)
    cloud_mask = ellipse_mask(base.observed.shape, base.wcs, cloud)
    samples = base.detect_foreground(noise_sigma=noise_sigma)
    if case.background_path is None:
        print("Estimating the SMF background", flush=True)
        background = base.estimate_background(method="smf", cloud=cloud).values
    else:
        print(f"Loading {case.background_path.name}", flush=True)
        background = base.set_background_from_fits(case.background_path).values

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
        result = mapper.compute_moderate(
            kappa_cm2_g=7.5,
            bright_pixel_policy="zero",
        )
    elif method == "liberal":
        result = mapper.compute_liberal(
            kappa_cm2_g=7.5,
            bright_pixel_policy="zero",
        )
    else:
        result = mapper.compute(
            kappa_cm2_g=7.5,
            bright_pixel_policy="zero",
        )

    crop = crop_bounds(cloud_mask)
    y0 = int(crop[0].start or 0)
    x0 = int(crop[1].start or 0)
    sample_rows = np.asarray(samples.rows - y0)
    sample_cols = np.asarray(samples.cols - x0)
    height, width = base.observed[crop].shape
    inside_crop = (
        (sample_rows >= 0)
        & (sample_rows < height)
        & (sample_cols >= 0)
        & (sample_cols < width)
    )
    strict = np.asarray(result.saturated_mask, dtype=bool) & cloud_mask
    foreground = np.asarray(result.foreground, dtype=float)
    local = cloud_mask & np.isfinite(base.observed) & (
        base.observed <= foreground + 2.0 * noise_sigma
    )
    metadata = {
        "cloud": case.label,
        "catalog_name": case.catalog_name,
        "method": method,
        "image": case.image_path.name,
        "background": (
            case.background_path.name
            if case.background_path is not None
            else "GTLMapping SMF estimate"
        ),
        "noise_sigma_mjy_sr": noise_sigma,
        "accepted_windows": int(samples.accepted_windows),
        "accepted_unique_sites": int(len(samples)),
        "surface_density_sum": positive_sum(result.surface_density, cloud_mask),
        "local_saturation_pixels": int(np.count_nonzero(local)),
        "strict_lower_limit_pixels": int(np.count_nonzero(strict)),
        "unresolved_pixels": int(np.count_nonzero(result.unresolved_mask & cloud_mask)),
        "detected_surface_density_sum": positive_sum(result.surface_density, cloud_mask & result.detection_mask),
        "limit_surface_density_sum": positive_sum(result.surface_density, cloud_mask & result.unresolved_mask),
        "samples_inside_aperture": int(np.count_nonzero(cloud_mask[samples.rows, samples.cols])),
        "fit_diagnostics": mapper.foreground_result.diagnostics,
        "foreground_min": float(np.nanmin(foreground[cloud_mask])),
        "foreground_median": float(np.nanmedian(foreground[cloud_mask])),
        "foreground_max": float(np.nanmax(foreground[cloud_mask])),
        "pixel_scale_arcsec": float(pixel_scale_arcsec(base.wcs)),
        "crop_x0": x0,
        "crop_y0": y0,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"cloud_{case.label.lower()}_{method}.npz"
    np.savez_compressed(
        output_path,
        observed=np.asarray(base.observed[crop], dtype=np.float64),
        background=np.asarray(background[crop], dtype=np.float64),
        cloud_mask=np.asarray(cloud_mask[crop], dtype=bool),
        foreground=np.asarray(foreground[crop], dtype=np.float64),
        surface_density=np.asarray(
            result.surface_density[crop].filled(np.nan), dtype=np.float64
        ),
        saturated=np.asarray(strict[crop], dtype=bool),
        unresolved=np.asarray(result.unresolved_mask[crop] & cloud_mask[crop], dtype=bool),
        bright=np.asarray(result.bright_mask[crop], dtype=bool),
        invalid_background=np.asarray(
            result.invalid_background_mask[crop], dtype=bool
        ),
        sample_rows=np.asarray(sample_rows[inside_crop], dtype=np.int32),
        sample_cols=np.asarray(sample_cols[inside_crop], dtype=np.int32),
        metadata=json.dumps(metadata),
    )
    print(json.dumps(metadata, indent=2), flush=True)
    print(f"Wrote {output_path}", flush=True)
    return output_path


def apply_plot_style(style_path: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    font_dir = style_path.resolve().parent / "fonts"
    font_path = font_dir / "AVHersheyComplexMedium.ttf"
    if not font_path.exists():
        raise FileNotFoundError(
            f"AVHershey Complex Medium was not found at {font_path}."
        )
    font_manager.fontManager.addfont(font_path)
    font_family = font_manager.FontProperties(fname=font_path).get_name()

    plt.style.use(style_path)
    mpl.rcParams.update(
        {
            "font.family": font_family,
            "font.sans-serif": [font_family],
            "mathtext.fontset": "custom",
            "mathtext.rm": font_family,
            "mathtext.it": font_family,
            "mathtext.bf": font_family,
            "mathtext.cal": font_family,
            "mathtext.sf": font_family,
            "mathtext.tt": font_family,
            "mathtext.fallback": "stix",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": False,
        }
    )


def save_figure(figure, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".png"), dpi=600, bbox_inches="tight")


def load_product(path: Path) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    with np.load(path) as product:
        arrays = {name: product[name] for name in product.files if name != "metadata"}
        metadata = json.loads(str(product["metadata"]))
    return arrays, metadata


def offset_extent(shape: tuple[int, int], scale: float) -> tuple[float, float, float, float]:
    height, width = shape
    return (
        -0.5 * width * scale,
        0.5 * width * scale,
        -0.5 * height * scale,
        0.5 * height * scale,
    )


def label_panel(
    axis,
    label: str,
    *,
    fontsize: float = 14,
    x: float = 0.025,
    y: float = 0.965,
) -> None:
    axis.text(
        x,
        y,
        label,
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=fontsize,
        fontweight="bold",
        color="black",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2.0},
        zorder=20,
    )


def plot_cloud(
    label: str,
    *,
    build_dir: Path,
    figure_dir: Path,
    style_path: Path,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.colors import AsinhNorm, LinearSegmentedColormap, Normalize, TwoSlopeNorm

    apply_plot_style(style_path)
    mpl.rcParams.update(
        {
            "axes.titlesize": 13,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 9,
            "lines.linewidth": 2.2,
        }
    )
    key = label.lower()
    bt12, bt12_meta = load_product(build_dir / f"cloud_{key}_bt12.npz")
    moderate, moderate_meta = load_product(build_dir / f"cloud_{key}_moderate.npz")
    if bt12["observed"].shape != moderate["observed"].shape:
        raise ValueError("BT12 and moderate products have different crops.")

    observed = bt12["observed"]
    mask = bt12["cloud_mask"].astype(bool)
    scale = float(bt12_meta["pixel_scale_arcsec"])
    extent = offset_extent(observed.shape, scale)
    sample_x = (bt12["sample_cols"] - 0.5 * observed.shape[1]) * scale
    sample_y = (bt12["sample_rows"] - 0.5 * observed.shape[0]) * scale
    change = np.where(mask, moderate["foreground"] - bt12["foreground"], np.nan)
    bt12_sigma = np.where(mask, bt12["surface_density"], np.nan)
    moderate_sigma = np.where(mask, moderate["surface_density"], np.nan)

    finite_observed = observed[np.isfinite(observed)]
    obs_low, obs_high = np.nanpercentile(finite_observed, (1.0, 99.5))
    observed_norm = Normalize(vmin=float(obs_low), vmax=float(obs_high))
    finite_change = np.abs(change[np.isfinite(change)])
    change_limit = max(float(np.nanpercentile(finite_change, 99.0)), 1.0e-6)
    change_norm = TwoSlopeNorm(vmin=-change_limit, vcenter=0.0, vmax=change_limit)
    diverging = LinearSegmentedColormap.from_list(
        "gtl_diverging", [STEEL, OFF_WHITE, TERRACOTTA]
    )
    density_values = np.concatenate(
        [
            bt12_sigma[np.isfinite(bt12_sigma)],
            moderate_sigma[np.isfinite(moderate_sigma)],
        ]
    )
    density_high = max(float(np.nanpercentile(density_values, 99.5)), 0.05)
    density_norm = AsinhNorm(
        linear_width=max(0.01 * density_high, 1.0e-4),
        vmin=0.0,
        vmax=density_high,
    )

    panel_aspect = observed.shape[1] / observed.shape[0]
    figure_height = max(4.6, min(7.4, 1.3 + 6.1 / panel_aspect))
    figure, axes = plt.subplots(
        2, 2, figsize=(9.0, figure_height), constrained_layout=True
    )
    observed_image = axes[0, 0].imshow(
        observed,
        origin="lower",
        extent=extent,
        cmap="gray",
        norm=observed_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[0, 0].contour(mask, levels=[0.5], colors="white", linewidths=1.0, extent=extent)
    axes[0, 0].scatter(
        sample_x,
        sample_y,
        s=18,
        facecolors="none",
        edgecolors=TEAL,
        linewidths=1.0,
        label=f"{len(sample_x)} sites shown",
    )
    axes[0, 0].legend(
        loc="upper right",
        fontsize=9,
        frameon=True,
        facecolor="white",
        framealpha=0.78,
    )
    axes[0, 0].set_title(r"Observed $I_{\rm obs}$")
    figure.colorbar(observed_image, ax=axes[0, 0], label=r"MJy sr$^{-1}$")

    change_image = axes[0, 1].imshow(
        change,
        origin="lower",
        extent=extent,
        cmap=diverging,
        norm=change_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[0, 1].contour(mask, levels=[0.5], colors="black", linewidths=0.8, extent=extent)
    axes[0, 1].set_title(r"$I_{\rm fg}^{\rm moderate}-I_{\rm fg}^{\rm BT12}$")
    figure.colorbar(change_image, ax=axes[0, 1], label=r"MJy sr$^{-1}$")

    bt12_image = axes[1, 0].imshow(
        bt12_sigma,
        origin="lower",
        extent=extent,
        cmap="inferno",
        norm=density_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[1, 0].contour(mask, levels=[0.5], colors="white", linewidths=0.8, extent=extent)
    limit_rows, limit_cols = np.where(bt12["unresolved"].astype(bool) & mask)
    axes[1, 0].scatter((limit_cols - 0.5*observed.shape[1])*scale,
                       (limit_rows - 0.5*observed.shape[0])*scale,
                       marker="x", s=9, color=SAFFRON, linewidths=0.6)
    axes[1, 0].set_title(r"BT12 $\Sigma$")
    figure.colorbar(bt12_image, ax=axes[1, 0], label=r"$\Sigma$ (g cm$^{-2}$)")

    moderate_image = axes[1, 1].imshow(
        moderate_sigma,
        origin="lower",
        extent=extent,
        cmap="inferno",
        norm=density_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[1, 1].contour(mask, levels=[0.5], colors="white", linewidths=0.8, extent=extent)
    sat_rows, sat_cols = np.where(moderate["unresolved"].astype(bool) & mask)
    if sat_rows.size:
        sat_x = (sat_cols - 0.5 * observed.shape[1]) * scale
        sat_y = (sat_rows - 0.5 * observed.shape[0]) * scale
        axes[1, 1].scatter(
            sat_x,
            sat_y,
            s=10,
            marker="x",
            color=SAFFRON,
            linewidths=0.8,
            label="Unresolved transmission",
        )
        axes[1, 1].legend(
            loc="upper right",
            fontsize=9,
            frameon=True,
            facecolor="white",
            framealpha=0.78,
        )
    axes[1, 1].set_title(r"Moderate GTL $\Sigma$")
    figure.colorbar(moderate_image, ax=axes[1, 1], label=r"$\Sigma$ (g cm$^{-2}$)")

    for index, (axis, panel) in enumerate(
        zip(axes.flat, ("(a)", "(b)", "(c)", "(d)"), strict=True)
    ):
        label_panel(axis, panel)
        if index >= 2:
            axis.set_xlabel(r"Offset (arcsec)")
        if index in (0, 2):
            axis.set_ylabel(r"Offset (arcsec)")
        axis.set_aspect("equal")
    save_figure(figure, figure_dir / f"fig_cloud_{key}")
    plt.close(figure)

    ratio = float(moderate_meta["surface_density_sum"]) / float(
        bt12_meta["surface_density_sum"]
    )
    print(
        f"Cloud {label}: moderate/BT12 sum={ratio:.8f}, "
        f"lower limits={int(moderate_meta['strict_lower_limit_pixels'])}"
    )


def write_profile_table(summary: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "cloud",
                "catalog_name",
                "profile",
                "surface_density_sum",
                "ratio_to_bt12",
                "local_saturation_pixels",
                "strict_lower_limit_pixels",
                "foreground_min_mjy_sr",
                "foreground_median_mjy_sr",
                "foreground_max_mjy_sr",
            ]
        )
        for cloud in summary:
            for profile, values in cloud["profiles"].items():
                writer.writerow(
                    [
                        cloud["cloud"],
                        cloud["catalog_name"],
                        profile,
                        values["surface_density_sum"],
                        values["ratio_to_bt12"],
                        values["local_saturation_pixels"],
                        values["strict_lower_limit_pixels"],
                        values["foreground_min"],
                        values["foreground_median"],
                        values["foreground_max"],
                    ]
                )


def log_probability_density(
    values: np.ndarray,
    bins: np.ndarray,
    *,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    histogram, _ = np.histogram(values, bins=bins, weights=weights)
    normalization = float(np.sum(histogram))
    if normalization <= 0.0:
        return np.zeros(len(bins) - 1, dtype=float)
    return histogram / (normalization * np.diff(np.log(bins)))


def plot_sigma_pdfs(
    *,
    summary_path: Path,
    build_dir: Path,
    data_dir: Path,
    figure_dir: Path,
    style_path: Path,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    apply_plot_style(style_path)
    mpl.rcParams.update(
        {
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "lines.linewidth": 2.0,
        }
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    write_profile_table(summary, data_dir / "profile_comparison.csv")

    profiles = ("bt12", "conservative", "moderate", "liberal")
    display_minimum = 0.01
    cloud_products: list[dict[str, object]] = []
    all_common_values: list[np.ndarray] = []
    for cloud in ("C", "F", "H"):
        products = {
            profile: load_product(
                build_dir / f"cloud_{cloud.lower()}_{profile}.npz"
            )[0]
            for profile in profiles
        }
        shapes = {products[profile]["surface_density"].shape for profile in profiles}
        if len(shapes) != 1:
            raise ValueError(f"Cloud {cloud} products have different crops.")

        cloud_mask = products["bt12"]["cloud_mask"].astype(bool)
        common = cloud_mask.copy()
        for profile in profiles:
            product = products[profile]
            sigma = np.asarray(product["surface_density"], dtype=float)
            common &= np.isfinite(sigma) & (sigma >= display_minimum)
            common &= ~product["unresolved"].astype(bool)
            common &= product["observed"] - product["foreground"] > 1.8
            common &= ~product["bright"].astype(bool)
            common &= ~product["invalid_background"].astype(bool)
        if not np.any(common):
            raise ValueError(f"Cloud {cloud} has no common PDF sample.")
        for profile in profiles:
            all_common_values.append(
                np.asarray(products[profile]["surface_density"], dtype=float)[common]
            )
        cloud_products.append(
            {"cloud": cloud, "products": products, "common": common}
        )

    maximum = max(float(np.nanmax(values)) for values in all_common_values)
    display_maximum = 10.0 ** np.ceil(np.log10(maximum * 1.000001))
    bins = np.geomspace(display_minimum, display_maximum, 33)
    centers = np.sqrt(bins[:-1] * bins[1:])

    data_dir.mkdir(parents=True, exist_ok=True)
    table_path = data_dir / "sigma_pdf_profiles.csv"
    figure, axes = plt.subplots(
        3,
        2,
        figsize=(8.2, 8.7),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    csv_rows: list[list[object]] = []
    panels = iter(("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
    for row, cloud_data in enumerate(cloud_products):
        cloud = str(cloud_data["cloud"])
        products = cloud_data["products"]
        common = np.asarray(cloud_data["common"], dtype=bool)
        common_count = int(np.count_nonzero(common))
        for column, weighting in enumerate(("area", "mass")):
            axis = axes[row, column]
            for profile in profiles:
                product = products[profile]
                sigma = np.asarray(product["surface_density"], dtype=float)[common]
                weights = None if weighting == "area" else sigma
                density = log_probability_density(sigma, bins, weights=weights)
                normalized_density = density / np.max(density)
                axis.stairs(
                    normalized_density,
                    bins,
                    color=PDF_PROFILE_COLORS[profile],
                    linestyle=PDF_PROFILE_LINESTYLES[profile],
                    linewidth=2.1 if profile == "bt12" else 1.9,
                    label="BT12" if profile == "bt12" else profile.capitalize(),
                )
                censored_count = int(
                    np.count_nonzero(
                        product["unresolved"].astype(bool)
                        & product["cloud_mask"].astype(bool)
                    )
                )
                for index, density_value in enumerate(density):
                    csv_rows.append(
                        [
                            cloud,
                            profile,
                            weighting,
                            bins[index],
                            bins[index + 1],
                            centers[index],
                            density_value,
                            normalized_density[index],
                            common_count,
                            censored_count,
                            display_minimum,
                        ]
                    )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_ylim(1.0e-3, 1.2)
            axis.grid(True, which="major", color="0.90", linewidth=0.65)
            label_panel(axis, next(panels), fontsize=11)
            axis.text(
                0.975,
                0.94,
                f"Cloud {cloud}",
                transform=axis.transAxes,
                ha="right",
                va="top",
                fontsize=11,
                fontweight="bold",
            )
            if column == 0:
                axis.set_ylabel(r"$p_A(\ln\Sigma)$ (normalized)")
            else:
                axis.set_ylabel(r"$p_M(\ln\Sigma)$ (normalized)")
            if row == 2:
                axis.set_xlabel(r"$\Sigma$ (g cm$^{-2}$)")
    axes[0, 0].set_title("Area weighted")
    axes[0, 1].set_title("Mass weighted")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="outside upper center",
        ncol=4,
        frameon=False,
    )
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "cloud",
                "profile",
                "weighting",
                "sigma_bin_left_g_cm2",
                "sigma_bin_right_g_cm2",
                "sigma_bin_center_g_cm2",
                "probability_density_per_ln_sigma",
                "probability_density_normalized_to_peak",
                "common_pixel_count",
                "profile_censored_pixel_count",
                "shared_display_minimum_g_cm2",
            ]
        )
        writer.writerows(csv_rows)
    save_figure(figure, figure_dir / "fig_sigma_pdfs")
    plt.close(figure)
    print(
        "Wrote common-sample PDFs with "
        + ", ".join(
            f"Cloud {item['cloud']}: "
            f"{int(np.count_nonzero(item['common'])):,} pixels"
            for item in cloud_products
        )
    )


def write_sgrc_tables(summary: dict[str, object], data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    with (data_dir / "sgrc_grid_sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        rows = summary["grid_sensitivity"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (data_dir / "sgrc_background_sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        rows = summary["background_sensitivity"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_sgrc(
    *,
    input_path: Path,
    bt12_path: Path,
    gtl_path: Path,
    summary_path: Path,
    data_dir: Path,
    figure_dir: Path,
    style_path: Path,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.colors import AsinhNorm, LinearSegmentedColormap, Normalize

    apply_plot_style(style_path)
    mpl.rcParams.update(
        {
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "lines.linewidth": 1.8,
            "lines.markersize": 6,
        }
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    write_sgrc_tables(summary, data_dir)
    x0, x1 = summary["x_range"]
    y0, y1 = summary["y_range"]
    with fits.open(input_path, memmap=True) as hdul:
        observed = np.asarray(hdul["SCI"].data[y0:y1, x0:x1], dtype=float)
    with fits.open(bt12_path, memmap=True) as hdul:
        bt12_sigma = np.asarray(hdul["SIGMA"].data, dtype=float)
        bt12_limits = np.asarray(hdul["UNRESOLVED"].data, dtype=bool)
    with fits.open(gtl_path, memmap=True) as hdul:
        gtl_sigma = np.asarray(hdul["SIGMA"].data, dtype=float)
        gtl_limits = np.asarray(hdul["UNRESOLVED"].data, dtype=bool)
    difference = np.where(bt12_limits | gtl_limits, np.nan, gtl_sigma - bt12_sigma)

    observed_norm = Normalize(vmin=2.0, vmax=15.0)
    density_high = max(float(np.nanpercentile(gtl_sigma, 99.0)), 0.05)
    density_norm = AsinhNorm(
        linear_width=max(0.01 * density_high, 1.0e-4), vmin=0.0, vmax=density_high
    )
    diff_high = max(float(np.nanpercentile(difference, 99.5)), 1.0e-5)
    sequential = LinearSegmentedColormap.from_list(
        "gtl_sequential", [OFF_WHITE, SAFFRON, TERRACOTTA]
    )

    figure, axes = plt.subplots(2, 3, figsize=(10.5, 6.8), constrained_layout=True)
    image = axes[0, 0].imshow(
        observed,
        origin="lower",
        cmap="gray",
        norm=observed_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[0, 0].set_title("F480M intensity")
    figure.colorbar(image, ax=axes[0, 0], label=r"MJy sr$^{-1}$")

    image = axes[0, 1].imshow(
        bt12_sigma,
        origin="lower",
        cmap="inferno",
        norm=density_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[0, 1].set_title(f"BT12: {summary['bt12_mass_msun']:.2f} " + r"$M_\odot$")
    rows, cols = np.where(bt12_limits)
    axes[0, 1].scatter(cols, rows, marker="x", s=3, color=SAFFRON, linewidths=0.5)
    figure.colorbar(image, ax=axes[0, 1], label=r"g cm$^{-2}$")

    image = axes[0, 2].imshow(
        gtl_sigma,
        origin="lower",
        cmap="inferno",
        norm=density_norm,
        interpolation="nearest",
        rasterized=True,
    )
    axes[0, 2].set_title(f"Conservative: {summary['gtl_lower_limit_mass_msun']:.2f} " + r"$M_\odot$")
    rows, cols = np.where(gtl_limits)
    axes[0, 2].scatter(cols, rows, marker="x", s=3, color=SAFFRON, linewidths=0.5)
    figure.colorbar(image, ax=axes[0, 2], label=r"g cm$^{-2}$")

    image = axes[1, 0].imshow(
        difference,
        origin="lower",
        cmap=sequential,
        vmin=0.0,
        vmax=diff_high,
        interpolation="nearest",
        rasterized=True,
    )
    axes[1, 0].set_title(r"Resolved $\Sigma_{\rm GTL}-\Sigma_{\rm BT12}$")
    figure.colorbar(image, ax=axes[1, 0], label=r"g cm$^{-2}$")

    grid = summary["grid_sensitivity"]
    grid_n = np.asarray([row["grid_n"] for row in grid])
    grid_ratio = np.asarray([row["gtl_to_bt12_mass_ratio"] for row in grid])
    axes[1, 1].plot(grid_n, grid_ratio, marker="o", color=TEAL)
    axes[1, 1].axhline(
        summary["gtl_to_bt12_lower_limit_mass_ratio"],
        color=TERRACOTTA,
        linestyle="--",
        linewidth=1.5,
        label="Fiducial",
    )
    axes[1, 1].set_xlabel("Scan grid, $n$")
    axes[1, 1].set_ylabel(r"$M_{\rm GTL}/M_{\rm BT12}$")
    axes[1, 1].set_title("Grid sensitivity")
    axes[1, 1].legend(loc="best", fontsize=9)

    background = summary["background_sensitivity"]
    background_level = np.asarray([row["background_level_mjy_sr"] for row in background])
    bt12_mass = np.asarray([row["bt12_mass_msun"] for row in background])
    gtl_mass = np.asarray([row["gtl_mass_msun"] for row in background])
    axes[1, 2].plot(
        background_level,
        bt12_mass,
        marker="s",
        color=PLUM,
        label="BT12",
    )
    axes[1, 2].plot(
        background_level,
        gtl_mass,
        marker="o",
        color=SAFFRON,
        label="Conservative GTL",
    )
    axes[1, 2].set_xlabel(r"$I_{\rm bg}$ (MJy sr$^{-1}$)")
    axes[1, 2].set_ylabel(r"Mass ($M_\odot$)")
    axes[1, 2].set_title("Background sensitivity")
    axes[1, 2].legend(loc="best", fontsize=9)

    for axis, panel in zip(
        axes.flat, ("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), strict=True
    ):
        label_panel(axis, panel, fontsize=11)
    for index, axis in enumerate(axes.flat[:4]):
        if index == 3:
            axis.set_xlabel("Image column")
        if index in (0, 3):
            axis.set_ylabel("Image row")
        axis.set_aspect("equal")
    save_figure(figure, figure_dir / "fig_sgrc_benchmark")
    plt.close(figure)


def parse_case(args: argparse.Namespace) -> CloudCase:
    background = None if args.background is None else Path(args.background)
    return CloudCase(
        label=args.cloud,
        catalog_name=args.catalog_name,
        image_path=Path(args.image),
        background_path=background,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare-cloud")
    prepare.add_argument("--cloud", required=True, choices=("C", "F", "H"))
    prepare.add_argument("--catalog-name", required=True)
    prepare.add_argument("--image", required=True, type=Path)
    prepare.add_argument("--background", type=Path)
    prepare.add_argument("--catalog", required=True, type=Path)
    prepare.add_argument(
        "--method",
        required=True,
        choices=("bt12", "conservative", "moderate", "liberal"),
    )
    prepare.add_argument("--noise-sigma", type=float, default=0.6)
    prepare.add_argument("--output-dir", required=True, type=Path)

    cloud_plot = subparsers.add_parser("plot-cloud")
    cloud_plot.add_argument("--cloud", required=True, choices=("C", "F", "H"))
    cloud_plot.add_argument("--build-dir", required=True, type=Path)
    cloud_plot.add_argument("--figure-dir", required=True, type=Path)
    cloud_plot.add_argument("--style", required=True, type=Path)

    summary_plot = subparsers.add_parser("plot-summary")
    summary_plot.add_argument("--summary", required=True, type=Path)
    summary_plot.add_argument("--build-dir", required=True, type=Path)
    summary_plot.add_argument("--data-dir", required=True, type=Path)
    summary_plot.add_argument("--figure-dir", required=True, type=Path)
    summary_plot.add_argument("--style", required=True, type=Path)

    sgrc_plot = subparsers.add_parser("plot-sgrc")
    sgrc_plot.add_argument("--input", required=True, type=Path)
    sgrc_plot.add_argument("--bt12", required=True, type=Path)
    sgrc_plot.add_argument("--gtl", required=True, type=Path)
    sgrc_plot.add_argument("--summary", required=True, type=Path)
    sgrc_plot.add_argument("--data-dir", required=True, type=Path)
    sgrc_plot.add_argument("--figure-dir", required=True, type=Path)
    sgrc_plot.add_argument("--style", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "prepare-cloud":
        prepare_cloud(
            parse_case(args),
            args.method,
            catalog_path=args.catalog,
            output_dir=args.output_dir,
            noise_sigma=args.noise_sigma,
        )
    elif args.command == "plot-cloud":
        plot_cloud(
            args.cloud,
            build_dir=args.build_dir,
            figure_dir=args.figure_dir,
            style_path=args.style,
        )
    elif args.command == "plot-summary":
        plot_sigma_pdfs(
            summary_path=args.summary,
            build_dir=args.build_dir,
            data_dir=args.data_dir,
            figure_dir=args.figure_dir,
            style_path=args.style,
        )
    elif args.command == "plot-sgrc":
        plot_sgrc(
            input_path=args.input,
            bt12_path=args.bt12,
            gtl_path=args.gtl,
            summary_path=args.summary,
            data_dir=args.data_dir,
            figure_dir=args.figure_dir,
            style_path=args.style,
        )


if __name__ == "__main__":
    main()
