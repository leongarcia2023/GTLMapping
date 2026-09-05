"""Compare foreground profiles with shared data, apertures and sensitivity."""
from pathlib import Path
import argparse
import gc
import json
from shutil import copyfile
import generate_figures as figures

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("catalog", "cloud-c", "cloud-c-background", "cloud-f", "cloud-h", "output-dir"):
        parser.add_argument("--" + flag, required=True, type=Path)
    parser.add_argument("--noise-sigma", type=float, default=0.6)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    build = args.output_dir / "build"
    style = Path(__file__).with_name("plotstyle.mplstyle")
    cases = [
        figures.CloudCase("C", "G028.37+00.07", args.cloud_c, args.cloud_c_background),
        figures.CloudCase("F", "G034.43+00.24", args.cloud_f, None),
        figures.CloudCase("H", "G035.39-00.33", args.cloud_h, None),
    ]
    rows = []
    for case in cases:
        row = {"cloud": case.label, "catalog_name": case.catalog_name,
               "image": case.image_path.name,
               "background": case.background_path.name if case.background_path else "GTLMapping SMF estimate",
               "profiles": {}}
        for profile in ("bt12", "conservative", "moderate", "liberal"):
            product = figures.prepare_cloud(case, profile, catalog_path=args.catalog,
                                            output_dir=build, noise_sigma=args.noise_sigma)
            _, meta = figures.load_product(product)
            row["profiles"][profile] = meta
            gc.collect()
        baseline = row["profiles"]["bt12"]["surface_density_sum"]
        for meta in row["profiles"].values():
            meta["ratio_to_bt12"] = meta["surface_density_sum"] / baseline
        for key in ("accepted_windows", "accepted_unique_sites", "noise_sigma_mjy_sr", "samples_inside_aperture"):
            row[key] = row["profiles"]["bt12"][key]
        rows.append(row)
        figures.plot_cloud(case.label, build_dir=build, figure_dir=args.output_dir, style_path=style)
        copyfile(args.output_dir / f"fig_cloud_{case.label.lower()}.png",
                 args.output_dir / f"cloud_{case.label.lower()}_method.png")
    summary = args.output_dir / "cloud_method_summary.json"
    summary.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    figures.plot_sigma_pdfs(summary_path=summary, build_dir=build, data_dir=args.output_dir,
                           figure_dir=args.output_dir, style_path=style)

if __name__ == "__main__":
    main()
