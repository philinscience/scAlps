"""Render a real Xenium h5ad without loading unrelated embeddings or all gene layers.

Example:
    python examples/xenium.py /path/to/GS52_spatial.h5ad --output outputs/GS52

Uses the supplied log1p_norm layer as-is. Native coordinate units are the default
because an obsm key alone does not establish a physical calibration.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from anndata import AnnData
from anndata.io import read_elem, sparse_dataset

import scalps as sca


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("outputs/xenium"))
    parser.add_argument("--genes", nargs="+", default=["Mki67", "Krt19", "Cd3e"])
    parser.add_argument("--layer", default="log1p_norm")
    parser.add_argument("--unit", default="native coordinate units")
    parser.add_argument("--density-scale", type=float, default=1)
    parser.add_argument("--density-percentile", type=float, default=1)
    parser.add_argument(
        "--gif", action="store_true", help="Also render a slow orbit of the first gene"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.input, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        indices = var.index.get_indexer(args.genes)
        if (indices < 0).any():
            raise ValueError(f"Genes not found: {np.asarray(args.genes)[indices < 0].tolist()}")
        matrix = sparse_dataset(handle[f"layers/{args.layer}"])[:, indices]
        adata = AnnData(matrix, obs=obs, var=var.iloc[indices].copy())
        adata.obsm["spatial"] = handle["obsm/spatial"][:]
        adata.layers[args.layer] = matrix
    sample = str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs else args.input.stem
    settings = dict(
        resolution=400,
        smooth=2.5,
        unit=args.unit,
        density_scale=args.density_scale,
        density_percentile=args.density_percentile,
    )
    visual = dict(height=0.16, skirt=False, contours=8, window_size=(1600, 1200))
    density_unit = (
        "cells/mm²"
        if args.unit == "µm" and args.density_scale == 1e6
        else (
            "cells / native coordinate unit²"
            if args.unit == "native coordinate units"
            else f"cells/{args.unit}² × {args.density_scale:g}"
        )
    )
    report = dict(
        source=str(args.input.resolve()),
        n_cells=adata.n_obs,
        genes=args.genes,
        expression_layer=args.layer,
        coordinate_unit=args.unit,
        expression_note="Supplied layer used as-is; no extra normalization or imputation.",
        plots=[],
    )
    for gene in args.genes:
        mountain = sca.terrain(adata, gene, layer=args.layer, **settings)
        path = mountain.save(
            args.output / f"{gene}.png",
            preset="ember",
            title=f"{sample} · {gene}",
            scalar_bar_title=f"Mean {gene} ({args.layer})",
            **visual,
        )
        report["plots"].append(dict(file=path.name, **mountain.metadata))
        print(path, flush=True)
        if gene == args.genes[0]:
            mountain.save(
                args.output / f"{gene}_transparent.png",
                preset="ember",
                title=f"{sample} · {gene}",
                transparent_background=True,
                **visual,
            )
            if args.gif:
                path = mountain.save(
                    args.output / f"{gene}.gif",
                    preset="ember",
                    title=f"{sample} · {gene}",
                    **dict(visual, window_size=(800, 600)),
                )
                print(path, flush=True)
    mountain = sca.terrain(adata, **settings)
    path = mountain.save(
        args.output / "cell_density.png",
        preset="glacier",
        title=f"{sample} · All-cell density",
        scalar_bar_title=f"Cell density ({density_unit})",
        **visual,
    )
    report["plots"].append(dict(file=path.name, **mountain.metadata))
    print(path, flush=True)
    if "Level_2" in adata.obs and "Lymphoid - T cells" in set(adata.obs["Level_2"]):
        mountain = sca.terrain(adata, groupby="Level_2", groups="Lymphoid - T cells", **settings)
        path = mountain.save(
            args.output / "T_cell_density.png",
            preset="glacier",
            title=f"{sample} · T-cell density",
            scalar_bar_title=f"T-cell density ({density_unit})",
            **visual,
        )
        report["plots"].append(dict(file=path.name, **mountain.metadata))
        print(path, flush=True)
    (args.output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
