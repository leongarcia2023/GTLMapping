"""Command-line entry point for a reproducible GTLMapping run."""

from __future__ import annotations

import argparse
from pathlib import Path

from .mapper import GTLMapper


def _hdu(value: str) -> int | str:
    try:
        return int(value)
    except ValueError:
        return value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gtlmapping",
        description="Create a MIREX mass-surface-density map.",
    )
    parser.add_argument("observed", type=Path, help="Observed MIR FITS image")
    parser.add_argument("background", type=Path, help="Aligned background FITS image")
    parser.add_argument("output", type=Path, help="Output multi-extension FITS file")
    parser.add_argument("--catalog", type=Path, required=True, help="Simon catalog.dat")
    parser.add_argument("--cloud", required=True, help="Cloud name, e.g. G028.37+00.07")
    parser.add_argument(
        "--foreground-method",
        choices=[
            "conservative",
            "moderate",
            "liberal",
            "kriging",
            "rbf",
            "spline",
            "gaussian",
            "cauchy",
            "flat",
            "bt12",
        ],
        default="conservative",
    )
    parser.add_argument("--noise-sigma", type=float, default=0.6)
    parser.add_argument("--foreground-margin", type=float, default=None)
    parser.add_argument(
        "--target-local-saturation-fraction",
        type=float,
        default=None,
        help="Override the moderate/liberal near-saturation budget",
    )
    parser.add_argument(
        "--maximum-strict-saturation-fraction",
        type=float,
        default=None,
        help="Override the moderate/liberal censored-pixel budget",
    )
    parser.add_argument(
        "--bt12-anchor-weight",
        type=float,
        default=None,
        help="Override the moderate/liberal soft BT12 weight",
    )
    parser.add_argument("--observed-hdu", type=_hdu, default=0)
    parser.add_argument("--uncertainty-hdu", type=_hdu)
    opacity = parser.add_mutually_exclusive_group()
    opacity.add_argument("--kappa", type=float)
    opacity.add_argument(
        "--filter",
        choices=["IRAC2", "F480M", "IRAC4", "F770W", "F2100W"],
    )
    parser.add_argument("--gas-to-dust-ratio", type=float, default=156.0)
    parser.add_argument("--kappa-std", type=float, default=0.0)
    parser.add_argument(
        "--bright-pixel-policy",
        choices=["allow", "zero", "mask"],
        default="allow",
    )
    parser.add_argument(
        "--saturation-policy",
        choices=["mask", "lower_limit"],
        default=None,
        help="Defaults to lower_limit for moderate/liberal and mask otherwise",
    )
    parser.add_argument("--intensity-floor", type=float)
    parser.add_argument(
        "--constrain-foreground",
        action="store_true",
        help="Constrain foreground below background by --intensity-floor",
    )
    parser.add_argument(
        "--background-fractional-uncertainty",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--align-background",
        action="store_true",
        help="Reproject the background (requires GTLMapping[align])",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line workflow."""

    args = _parser().parse_args(argv)
    mapper = GTLMapper.from_fits(
        args.observed,
        hdu=args.observed_hdu,
        uncertainty_hdu=args.uncertainty_hdu,
    )
    cloud = mapper.select_cloud(args.catalog, args.cloud)
    margin = (
        2.0 * args.noise_sigma
        if args.foreground_margin is None
        else args.foreground_margin
    )
    if args.foreground_method == "bt12":
        foreground = mapper.fit_foreground(
            method="bt12",
            cloud=cloud,
            noise_sigma=args.noise_sigma,
            foreground_margin=margin,
        )
        print(
            "BT12 saturated pixels: "
            f"{foreground.diagnostics['saturated_pixel_count']}; "
            "constant foreground: "
            f"{foreground.values[0, 0]:.6g}."
        )
    else:
        samples = mapper.detect_foreground(noise_sigma=args.noise_sigma)
        print(
            f"Foreground samples: {samples.raw_detection_count} accepted windows, "
            f"{len(samples)} unique coordinates."
        )
        if args.foreground_method == "conservative":
            mapper.fit_foreground(
                method="conservative",
                samples=samples,
                cloud=cloud,
                noise_sigma=args.noise_sigma,
                foreground_margin=margin,
            )
        elif args.foreground_method in {"moderate", "liberal"}:
            defaults = (
                (0.005, 0.0001, 0.5)
                if args.foreground_method == "moderate"
                else (0.01, 0.001, 0.0)
            )
            local_fraction = (
                defaults[0]
                if args.target_local_saturation_fraction is None
                else args.target_local_saturation_fraction
            )
            strict_fraction = (
                defaults[1]
                if args.maximum_strict_saturation_fraction is None
                else args.maximum_strict_saturation_fraction
            )
            anchor_weight = (
                defaults[2]
                if args.bt12_anchor_weight is None
                else args.bt12_anchor_weight
            )
            foreground = mapper.fit_foreground(
                method=args.foreground_method,
                samples=samples,
                cloud=cloud,
                noise_sigma=args.noise_sigma,
                foreground_margin=margin,
                target_local_saturation_fraction=local_fraction,
                maximum_strict_saturation_fraction=strict_fraction,
                bt12_anchor_weight=anchor_weight,
            )
            print(
                f"{args.foreground_method.title()} near-saturated pixels: "
                f"{foreground.diagnostics['local_saturation_count']} "
                "(BT12 comparison: "
                f"{foreground.diagnostics['reference_local_saturation_count']}); "
                "strict lower limits: "
                f"{foreground.diagnostics['strict_saturation_count']}."
            )
        else:
            mapper.fit_foreground(
                method=args.foreground_method,
                samples=samples,
                foreground_margin=margin,
            )
    mapper.set_background_from_fits(
        args.background,
        align=args.align_background,
    )
    saturation_policy = args.saturation_policy
    if saturation_policy is None:
        saturation_policy = (
            "lower_limit"
            if args.foreground_method in {"moderate", "liberal"}
            else "mask"
        )
    intensity_floor = args.intensity_floor
    if (
        args.foreground_method in {"moderate", "liberal"}
        and intensity_floor is None
    ):
        intensity_floor = 2.0 * args.noise_sigma
    if args.foreground_method in {"moderate", "liberal"}:
        mapper.constrain_foreground(
            minimum_transmitted_intensity=intensity_floor,
            minimum_foreground=0.0,
        )
    elif args.constrain_foreground:
        if args.intensity_floor is None:
            raise ValueError(
                "--constrain-foreground requires --intensity-floor."
            )
        mapper.constrain_foreground(
            minimum_transmitted_intensity=args.intensity_floor,
        )
    background_std = (
        None
        if args.background_fractional_uncertainty == 0
        else (
            args.background_fractional_uncertainty
            * mapper.background_result.values
        )
    )
    result = mapper.compute(
        kappa_cm2_g=args.kappa,
        filter_name=args.filter,
        gas_to_dust_ratio=args.gas_to_dust_ratio,
        kappa_std_cm2_g=args.kappa_std,
        bright_pixel_policy=args.bright_pixel_policy,
        saturation_policy=saturation_policy,
        intensity_floor=intensity_floor,
        background_std=background_std,
    )
    result.write(args.output, overwrite=args.overwrite)
    valid_fraction = float(
        1.0 - result.surface_density.mask.mean()
    )
    print(f"Wrote {args.output} ({valid_fraction:.2%} valid pixels).")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
