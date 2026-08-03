"""Run the staged GTLMapping workflow on a Cloud C image."""

from __future__ import annotations

import argparse
from pathlib import Path

from gtlmapping import GTLMapper


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("observed", type=Path)
    parser.add_argument("catalog", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    mapper = GTLMapper.from_fits(args.observed)
    cloud = mapper.select_cloud(args.catalog, "G028.37+00.07")
    samples = mapper.detect_foreground()
    print(
        f"Accepted {samples.raw_detection_count} windows at "
        f"{len(samples)} unique coordinates."
    )
    foreground = mapper.fit_foreground(
        method="conservative",
        cloud=cloud,
        noise_sigma=0.6,
    )
    print(
        "Conservative spatial-trend blend factor: "
        f"{foreground.diagnostics['blend_factor']:.3f}."
    )
    mapper.estimate_background(method="smf", cloud=cloud)
    result = mapper.compute(kappa_cm2_g=7.5)
    result.write(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
