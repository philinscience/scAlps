"""A small command line entry point for demos and h5ad files."""

import argparse

from .core import terrain
from .datasets import demo


def main():
    parser = argparse.ArgumentParser(description="scAlps · turn spatial cells into landscapes")
    parser.add_argument("input", help="Path to .h5ad, or 'demo' for synthetic tissue")
    parser.add_argument("--value", help="Gene or numeric obs column; omit for density")
    parser.add_argument("--groupby")
    parser.add_argument("--groups", nargs="+")
    parser.add_argument("--layer")
    parser.add_argument("--resolution", type=int, default=250)
    parser.add_argument("--smooth", type=float, default=2)
    parser.add_argument("--preset", choices=["alpine", "ember", "glacier"], default="alpine")
    parser.add_argument("--output", "-o", help=".png, .html, .gif, .vtp or .npz; omit to explore")
    args = parser.parse_args()
    if args.input == "demo":
        adata = demo()
    else:
        from anndata import read_h5ad

        adata = read_h5ad(args.input, backed="r")
    try:
        landscape = terrain(
            adata,
            args.value,
            groupby=args.groupby,
            groups=args.groups,
            layer=args.layer,
            resolution=args.resolution,
            smooth=args.smooth,
        )
        if args.output:
            options = {} if args.output.endswith((".npz", ".vtp")) else {"preset": args.preset}
            print(landscape.save(args.output, **options))
        else:
            landscape.plot(preset=args.preset)
    finally:
        if adata.isbacked:
            adata.file.close()
