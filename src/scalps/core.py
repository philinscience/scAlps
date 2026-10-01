"""Numerical terrain construction, independent of the renderer."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from ._input import resolve, values


@dataclass
class Terrain:
    """Gridded data in original coordinate units, with auditable support and metadata."""

    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    mask: np.ndarray
    counts: np.ndarray
    label: str
    metadata: dict = field(default_factory=dict)

    def mesh(self, *, height=0.22, vmax=None, transform="linear"):
        """Create a PyVista surface. Height affects geometry; scalar colors remain raw."""
        import pyvista as pv

        if not np.isfinite(height) or height < 0:
            raise ValueError("height must be finite and nonnegative.")
        if transform not in {"linear", "sqrt", "log1p"}:
            raise ValueError("transform must be 'linear', 'sqrt' or 'log1p'.")
        raw = np.where(self.mask, self.values, 0.0)
        scale = np.percentile(np.abs(raw[self.mask]), 99) if vmax is None else float(vmax)
        if not np.isfinite(scale) or scale < 0 or (vmax is not None and scale == 0):
            raise ValueError("vmax must be finite and positive.")
        scale = scale or 1.0
        z = np.clip(raw / scale, -1, 1)
        if transform == "sqrt":
            z = np.sign(z) * np.sqrt(np.abs(z))
        elif transform == "log1p":
            z = np.sign(z) * np.log1p(9 * np.abs(z)) / np.log(10)
        z *= max(np.ptp(self.x), np.ptp(self.y)) * height
        xx, yy = np.meshgrid(self.x, self.y, indexing="ij")
        grid = pv.StructuredGrid(xx, yy, z)
        grid.point_data["value"] = raw.ravel(order="F")
        grid.point_data["support"] = self.mask.astype(float).ravel(order="F")
        # Drop quads touching unsupported vertices; no bridges across large tissue gaps.
        surface = grid.threshold(0.5, scalars="support", preference="point", all_scalars=True)
        if not surface.n_cells:
            raise ValueError("No connected terrain; increase smooth or reduce resolution/support.")
        return surface.extract_surface(algorithm="dataset_surface").triangulate()

    def plot(self, **kwargs):
        """Build a styled PyVista plotter; see scalps.plot."""
        from .render import plot

        return plot(self, **kwargs)

    def save(self, path, **kwargs):
        """Export PNG, interactive HTML, orbit GIF, VTP mesh or numerical NPZ."""
        from .render import save

        return save(self, Path(path), **kwargs)


def terrain(
    data,
    value=None,
    *,
    groupby=None,
    groups=None,
    statistic=None,
    layer=None,
    spatial_key="spatial",
    table=None,
    element=None,
    coordinate_system="global",
    resolution=250,
    smooth=2.0,
    support=0.02,
    unit="µm",
    density_scale=1_000_000.0,
):
    """Aggregate cells into a terrain with Gaussian smoothing and explicit tissue support.

    ``value`` is an obs column, gene, numeric array, or None (density).
    ``resolution`` sets the longest grid dimension; ``smooth`` is in grid pixels.
    Density defaults to cells/mm² when coordinates are micrometers. Means use
    smoothed sums / smoothed valid counts, rather than means of occupied bins.
    ``groups`` selects cells but retains the full tissue footprint.
    """
    if isinstance(resolution, bool) or not isinstance(resolution, (int, np.integer)):
        raise ValueError("resolution must be an integer between 16 and 2000.")
    if not 16 <= resolution <= 2000:
        raise ValueError("resolution must be between 16 and 2000.")
    if not np.isfinite(smooth) or smooth < 0:
        raise ValueError("smooth must be finite and nonnegative.")
    if not np.isfinite(support) or not 0 < support <= 1:
        raise ValueError("support must be in (0, 1].")
    if not np.isfinite(density_scale) or density_scale <= 0:
        raise ValueError("density_scale must be finite and positive.")
    adata, xy = resolve(data, spatial_key, table, element, coordinate_system)
    val, label = values(adata, value, layer)
    statistic = statistic or ("density" if val is None else "mean")
    if statistic not in {"density", "mean", "sum", "fraction"}:
        raise ValueError("statistic must be density, mean, sum or fraction.")
    if (statistic in {"density", "fraction"}) != (val is None):
        raise ValueError("density/fraction require value=None; mean/sum require a numeric value.")
    selected = np.ones(len(xy), dtype=bool)
    if (groupby is None) != (groups is None):
        raise ValueError("Provide both groupby= and groups=.")
    if groupby is not None:
        groups = [groups] if np.isscalar(groups) else list(groups)
        missing = set(groups) - set(adata.obs[groupby].dropna())
        if missing:
            raise ValueError(f"Unknown groups: {missing}")
        selected = adata.obs[groupby].isin(groups).to_numpy()
        label = ", ".join(map(str, groups)) + (f" · {label}" if val is not None else "")
    if not selected.any():
        raise ValueError("No cells selected.")
    if statistic == "fraction" and groupby is None:
        raise ValueError("fraction requires groupby= and groups=.")
    if val is not None:
        selected &= np.isfinite(val)
        if not selected.any():
            raise ValueError("No finite values in the selected cells.")
    span = np.ptp(xy, axis=0)
    spacing = span.max() / (resolution - 1)
    pad = int(np.ceil(3 * smooth)) + 1
    shape = np.maximum(2, np.ceil(span / spacing).astype(int) + 1)
    edges = [
        xy[:, i].min() + (np.arange(shape[i] + 2 * pad + 1) - pad - 0.5) * spacing for i in range(2)
    ]
    hist = lambda points, weights=None: np.histogram2d(  # noqa: E731
        points[:, 0], points[:, 1], bins=edges, weights=weights
    )[0]
    blur = lambda a: gaussian_filter(a, smooth, mode="constant")  # noqa: E731
    counts = hist(xy)
    occupancy = blur((counts > 0).astype(float))
    mask = occupancy >= support
    denominator = blur(hist(xy[selected]))
    if statistic == "density":
        result = denominator * density_scale / spacing**2
        label += " · cells/mm²" if unit == "µm" and density_scale == 1e6 else " · density"
    elif statistic == "fraction":
        total = blur(counts)
        result = np.divide(denominator, total, out=np.zeros_like(total), where=total > 1e-12)
        label += " · fraction"
    else:
        numerator = blur(hist(xy[selected], val[selected]))
        if statistic == "mean":
            mask &= denominator > 1e-12
            result = np.divide(
                numerator, denominator, out=np.zeros_like(numerator), where=denominator > 1e-12
            )
        else:
            result = numerator
    if not mask.any():
        raise ValueError(
            "No tissue support; reduce support, reduce resolution, or increase smooth."
        )
    result[~mask] = np.nan
    return Terrain(
        (edges[0][:-1] + edges[0][1:]) / 2,
        (edges[1][:-1] + edges[1][1:]) / 2,
        result,
        mask,
        counts,
        label,
        dict(
            statistic=statistic,
            smooth=float(smooth),
            support=float(support),
            resolution=int(resolution),
            spacing=float(spacing),
            unit=unit,
            density_scale=float(density_scale),
            n_cells=len(xy),
            n_selected=int(selected.sum()),
            spatial_key=spatial_key,
            table=table,
            element=element,
            coordinate_system=coordinate_system if element is not None else None,
            value=value if isinstance(value, str) else ("array" if value is not None else None),
            layer=layer,
            groupby=groupby,
            groups=None if groups is None else list(map(str, groups)),
        ),
    )
