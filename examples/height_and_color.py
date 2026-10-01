"""T-cell density for height, a separate score or gene for color.

Without --input, renders a synthetic activation example. With --input, reads
only Cd8a (or --color-gene) from a Xenium h5ad's log1p_norm layer. Coordinates
in that file remain in native units. No activation score is inferred.
"""

import argparse
from pathlib import Path

import scalps as sca


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--color-gene", default="Cd8a")
    parser.add_argument("--layer", default="log1p_norm")
    parser.add_argument("--groupby", default="Level_2")
    parser.add_argument("--group", default="Lymphoid - T cells")
    parser.add_argument("--output", type=Path, default=Path("outputs/height-and-color"))
    parser.add_argument("--gif", action="store_true")
    args = parser.parse_args()
    if args.input:
        import h5py
        from anndata import AnnData
        from anndata.io import read_elem, sparse_dataset

        with h5py.File(args.input, "r") as handle:
            var = read_elem(handle["var"])
            indices = var.index.get_indexer([args.color_gene])
            if (indices < 0).any():
                raise ValueError(f"Gene not found: {args.color_gene!r}")
            matrix = sparse_dataset(handle[f"layers/{args.layer}"])[:, indices]
            adata = AnnData(matrix, obs=read_elem(handle["obs"]), var=var.iloc[indices].copy())
            adata.layers[args.layer] = matrix
            adata.obsm["spatial"] = handle["obsm/spatial"][:]
        options = dict(
            color=f"gene:{args.color_gene}",
            color_layer=args.layer,
            groupby=args.groupby,
            groups=args.group,
            unit="native coordinate units",
            density_scale=1,
        )
        sample = str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs else args.input.stem
        title = f"{sample} · T-cell density and {args.color_gene}"
        height_label = "T-cell density (cells / native coordinate unit²)"
        color_label = f"Mean {args.color_gene} in T cells ({args.layer})"
    else:
        adata = sca.demo()
        options = dict(color="activation_score", groupby="cell_type", groups="T cell")
        title = "T-cell density and activation · Synthetic tissue"
        height_label = "T-cell density (cells/mm²)"
        color_label = "Mean activation score in T cells (synthetic)"
    mountain = sca.terrain(adata, resolution=350, smooth=2.5, **options)
    if args.input:
        mountain.metadata["source"] = str(args.input.resolve())
    visual = dict(
        cmap="viridis",
        height=0.16,
        skirt=False,
        contours=8,
        title=title,
        height_label=height_label,
        scalar_bar_title=color_label,
    )
    if not args.input:
        visual["clim"] = (0, 1)
    print(
        mountain.save(args.output / "terrain.png", window_size=(1600, 1200), **visual), flush=True
    )
    mountain.save(args.output / "terrain.npz")
    if args.gif:
        print(
            mountain.save(args.output / "terrain.gif", window_size=(900, 700), **visual), flush=True
        )


if __name__ == "__main__":
    main()
