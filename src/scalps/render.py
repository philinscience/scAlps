"""PyVista styling and exports, kept separate from data aggregation."""

import inspect
import json
from pathlib import Path

import numpy as np

from .core import Terrain, terrain

PRESETS = {
    "alpine": dict(cmap="terrain", background="#e8edf1", ink="#243944", rock="#9aa7aa"),
    "ember": dict(cmap="magma", background="#101420", ink="#e9e7ef", rock="#222735"),
    "glacier": dict(cmap="viridis", background="#091f2c", ink="#d6f1f4", rock="#183d4b"),
}


def _colormap(preset):
    from matplotlib.colors import LinearSegmentedColormap

    palettes = {
        "alpine": [
            "#26485a",
            "#457986",
            "#7ca79b",
            "#bbc696",
            "#cbb791",
            "#9a8070",
            "#ece6da",
            "#ffffff",
        ],
        "glacier": ["#0c263c", "#145776", "#278e9c", "#78d8c2", "#e9ffef"],
    }
    if preset in palettes:
        return LinearSegmentedColormap.from_list(f"scalps_{preset}", palettes[preset])
    return PRESETS[preset]["cmap"]


def plot(
    data,
    value=None,
    *,
    preset="alpine",
    height=0.22,
    transform="linear",
    vmax=None,
    clim=None,
    cmap=None,
    contours=12,
    skirt=True,
    title=None,
    axes=False,
    scalar_bar=True,
    image_coordinates=True,
    off_screen=None,
    window_size=(1400, 1000),
    show=True,
    **terrain_kwargs,
):
    """Plot AnnData, SpatialData or a Terrain; return the editable PyVista Plotter.

    Use ``show=False`` to customize it, or ``off_screen=True, show=False`` on a server.
    ``image_coordinates=True`` displays increasing y downward (Xenium convention).
    ``clim`` controls raw color limits; ``vmax`` controls absolute geometric scaling.
    """
    import pyvista as pv

    if preset not in PRESETS:
        raise ValueError(f"Unknown preset {preset!r}; choose from {list(PRESETS)}.")
    if isinstance(data, Terrain):
        if value is not None or terrain_kwargs:
            raise ValueError("Aggregation arguments cannot be applied to an existing Terrain.")
        landscape = data
    else:
        landscape = terrain(data, value, **terrain_kwargs)
    style = PRESETS[preset]
    mesh = landscape.mesh(height=height, vmax=vmax, transform=transform)
    if image_coordinates:
        mesh.points[:, 1] *= -1
        mesh.flip_faces(inplace=True)
    finite = landscape.values[landscape.mask]
    signed = finite.min() < 0
    if clim is None:
        if signed:
            limit = max(float(np.max(np.abs(finite))), 1e-12)
            clim = (-limit, limit)
        else:
            clim = (0, max(float(finite.max()), 1e-12))
    p = pv.Plotter(off_screen=off_screen, window_size=window_size)
    p.set_background(style["background"])
    p.add_mesh(
        mesh,
        scalars="value",
        cmap=cmap or ("RdBu_r" if signed else _colormap(preset)),
        clim=clim,
        smooth_shading=True,
        ambient=0.25,
        diffuse=0.8,
        specular=0.12,
        specular_power=30,
        show_scalar_bar=scalar_bar,
        scalar_bar_args=dict(
            title=landscape.label,
            color=style["ink"],
            title_font_size=15,
            label_font_size=12,
            fmt="%.2g",
            position_x=0.27,
            position_y=0.045,
            width=0.46,
            height=0.085,
        ),
    )
    span = max(np.ptp(landscape.x), np.ptp(landscape.y))
    if skirt:
        # Build a vertical curtain from each boundary edge to a common base plane.
        boundary = mesh.extract_feature_edges(
            boundary_edges=True, feature_edges=False, manifold_edges=False, non_manifold_edges=False
        )
        if boundary.n_cells:
            top = boundary.points.copy()
            bottom = top.copy()
            bottom[:, 2] = min(0, mesh.bounds[4]) - span * 0.035
            n = len(top)
            segments = boundary.lines.reshape(-1, 3)[:, 1:]
            faces = np.column_stack(
                [
                    np.full(len(segments), 4),
                    segments[:, 0],
                    segments[:, 1],
                    segments[:, 1] + n,
                    segments[:, 0] + n,
                ]
            )
            curtain = pv.PolyData(np.vstack([top, bottom]), faces.ravel())
            p.add_mesh(curtain, color=style["rock"], smooth_shading=False, ambient=0.3)
    if contours and finite.max() > finite.min():
        levels = np.linspace(finite.min(), finite.max(), int(contours) + 2)[1:-1]
        lines = mesh.contour(levels, scalars="value")
        if lines.n_points:
            lines.points[:, 2] += span * 0.0004
            p.add_mesh(lines, color=style["ink"], opacity=0.22, line_width=1, show_scalar_bar=False)
    p.add_text(
        title or landscape.label,
        position=(0.045, 0.91),
        viewport=True,
        font_size=20,
        color=style["ink"],
        font="arial",
    )
    p.add_text(
        "scAlps  /  SPATIAL LANDSCAPES",
        position=(0.045, 0.035),
        viewport=True,
        font_size=9,
        color=style["ink"],
    )
    if axes:
        p.show_bounds(
            xtitle=f"x ({landscape.metadata['unit']})",
            ytitle=("−y" if image_coordinates else "y") + f" ({landscape.metadata['unit']})",
            ztitle="visual height",
            color=style["ink"],
            grid=False,
        )
    p.view_isometric()
    p.enable_parallel_projection()
    p.camera.elevation = 12
    p.camera.zoom(1.25)
    p.enable_anti_aliasing("ssaa")
    if show:
        p.show()
    return p


def save(landscape, path, *, frames=90, fps=24, **kwargs):
    """Write an export and a JSON sidecar recording scientific and visual settings."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in {".npz", ".vtp", ".png", ".html", ".gif"}:
        raise ValueError("Choose .png, .html, .gif, .vtp or .npz.")
    if suffix == ".gif" and (frames < 2 or fps <= 0):
        raise ValueError("GIF exports need frames >= 2 and fps > 0.")
    path.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".npz":
        np.savez_compressed(
            path,
            x=landscape.x,
            y=landscape.y,
            values=landscape.values,
            mask=landscape.mask,
            counts=landscape.counts,
        )
    elif suffix == ".vtp":
        landscape.mesh(**kwargs).save(path)
    else:
        p = plot(landscape, off_screen=True, show=False, **kwargs)
        try:
            if suffix == ".png":
                p.screenshot(path)
            elif suffix == ".html":
                try:
                    import trame_pyvista  # noqa: F401

                    p.trame.export_html(path)
                except ImportError as exc:
                    raise ImportError("HTML export needs: pip install 'scalps[notebook]'") from exc
            else:
                p.open_gif(path, fps=fps)
                for _ in range(frames):
                    p.write_frame()
                    p.camera.azimuth += 360 / frames
        finally:
            p.close()
    render_options = dict(kwargs)
    if suffix in {".png", ".html", ".gif"}:
        render_options = {
            name: param.default
            for name, param in inspect.signature(plot).parameters.items()
            if param.default is not inspect.Parameter.empty and name != "value"
        }
        render_options.update(kwargs, off_screen=True, show=False)
    metadata = dict(
        landscape.metadata,
        label=landscape.label,
        render=render_options,
        scalps_version="0.1.0",
        value_range=[float(np.nanmin(landscape.values)), float(np.nanmax(landscape.values))],
    )
    if suffix == ".gif":
        metadata.update(frames=frames, fps=fps)
    path.with_suffix(path.suffix + ".json").write_text(
        json.dumps(metadata, indent=2, default=str) + "\n"
    )
    return path
