"""Render the Recipes gallery from reproducible synthetic data.

Run ``python examples/recipes.py`` for PNGs and their settings sidecars.
Add ``--exports`` to also write a transparent PNG, 24-second GIF, interactive
HTML, numerical NPZ and VTP mesh to outputs/recipes. HTML needs the notebook extra.
"""

import argparse
from pathlib import Path

import scalps as sca

VIEW = dict(height=0.16, skirt=False, contours=8, window_size=(1200, 900))

# Aggregation and display settings mirror the examples in docs/recipes.md.
RECIPES = {
    "cell-density": ({}, dict(preset="glacier", title="Cell density")),
    "gene-expression": (
        dict(value="MKI67", layer="counts"),
        dict(preset="ember", title="MKI67 expression"),
    ),
    "score": (dict(value="prolif_score"), dict(transform="sqrt", title="Proliferation score")),
    "signed-score": (
        dict(value="immune_balance"),
        dict(cmap="RdBu_r", title="Signed immune balance"),
    ),
    "tcell-density": (
        dict(groupby="cell_type", groups="T cell"),
        dict(preset="glacier", title="T-cell density"),
    ),
    "tcell-fraction": (
        dict(groupby="cell_type", groups="T cell", statistic="fraction"),
        dict(preset="glacier", clim=(0, 1), vmax=1, title="Local T-cell fraction"),
    ),
    "gene-sum": (
        dict(value="MKI67", layer="counts", statistic="sum"),
        dict(preset="ember", title="MKI67 sum", scalar_bar_title="Smoothed MKI67 sum per grid bin"),
    ),
    "smoothing-sharp": (
        dict(value="MKI67", layer="counts", resolution=250, smooth=1),
        dict(preset="ember", vmax=20, clim=(0, 20), title="MKI67 · Gaussian sigma = 1 pixel"),
    ),
    "smoothing-broad": (
        dict(value="MKI67", layer="counts", resolution=250, smooth=4),
        dict(preset="ember", vmax=20, clim=(0, 20), title="MKI67 · Gaussian sigma = 4 pixels"),
    ),
    "height-color": (
        dict(groupby="cell_type", groups="T cell", color="activation_score"),
        dict(
            cmap="viridis",
            clim=(0, 1),
            title="T-cell density and activation",
            height_label="T-cell density (cells/mm²)",
            scalar_bar_title="Mean activation score in T cells (synthetic)",
        ),
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs/assets/recipes",
    )
    parser.add_argument("--exports", action="store_true")
    args = parser.parse_args()
    adata = sca.demo(seed=7)
    for name, (aggregation, display) in RECIPES.items():
        mountain = sca.terrain(adata, **aggregation)
        visual = dict(VIEW, **display)
        visual["title"] += " · Synthetic tissue"
        print(mountain.save(args.output / f"{name}.png", **visual), flush=True)
    if args.exports:
        mountain = sca.terrain(adata, "MKI67", layer="counts")
        out = Path("outputs/recipes")
        visual = dict(VIEW, preset="ember", title="MKI67 · Synthetic tissue")
        print(mountain.save(out / "MKI67.npz"), flush=True)
        print(mountain.save(out / "MKI67.vtp", height=VIEW["height"]), flush=True)
        print(
            mountain.save(out / "MKI67-transparent.png", transparent_background=True, **visual),
            flush=True,
        )
        print(mountain.save(out / "MKI67.html", **visual), flush=True)
        print(mountain.save(out / "MKI67.gif", **dict(visual, window_size=(700, 500))), flush=True)


if __name__ == "__main__":
    main()
