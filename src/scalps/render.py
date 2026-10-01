"""PyVista styling and exports, kept separate from data aggregation."""

import inspect
import json
from pathlib import Path

import numpy as np

from .core import Terrain, terrain

PRESETS = {
    "alpine": dict(cmap="terrain", rock="#9aa7aa"),
    "ember": dict(cmap="magma", rock="#8b8991"),
    "glacier": dict(cmap="viridis", rock="#78969e"),
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
    footprint=True,
    footprint_color="#d6d6d6",
    footprint_gap=0.06,
    title=None,
    background="white",
    axes=False,
    scalar_bar=True,
    scalar_bar_title=None,
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
    All presets use a white background by default. The legend has a dedicated
    viewport and uses regular DejaVu Sans, independent of the global PyVista theme.
    """
    import pyvista as pv
    from matplotlib.font_manager import findfont

    if preset not in PRESETS:
        raise ValueError(f"Unknown preset {preset!r}; choose from {list(PRESETS)}.")
    if not np.isfinite(footprint_gap) or footprint_gap <= 0:
        raise ValueError("footprint_gap must be finite and positive.")
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
    # Keep the legend outside the 3D viewport, even when the camera is rotated.
    layout = dict(shape=(2, 1), row_weights=(0.84, 0.16)) if scalar_bar else {}
    p = pv.Plotter(off_screen=off_screen, window_size=window_size, border=False, **layout)
    p.set_background(background, all_renderers=True)
    font_file = findfont("DejaVu Sans")
    ink = "#202020"
    font_scale = min(window_size) / 1000
    label_size = max(9, round(20 * font_scale))
    title_size = max(10, round(22 * font_scale))
    actor = p.add_mesh(
        mesh,
        scalars="value",
        cmap=cmap or ("RdBu_r" if signed else _colormap(preset)),
        clim=clim,
        smooth_shading=True,
        ambient=0.25,
        diffuse=0.8,
        specular=0.12,
        specular_power=30,
        show_scalar_bar=False,
    )
    span = max(np.ptp(landscape.x), np.ptp(landscape.y))
    if footprint:
        # Below both valleys and any decorative skirt; colors encode no extra data.
        floor_z = min(0, mesh.bounds[4]) - span * (footprint_gap + (0.035 if skirt else 0))
        floor = landscape.footprint(z=floor_z)
        if image_coordinates:
            floor.points[:, 1] *= -1
            floor.flip_faces(inplace=True)
        p.add_mesh(
            floor,
            color=footprint_color,
            lighting=False,
            show_scalar_bar=False,
            name="tissue-footprint",
        )
        outline = floor.extract_feature_edges(
            boundary_edges=True, feature_edges=False, manifold_edges=False, non_manifold_edges=False
        )
        if outline.n_points:
            p.add_mesh(
                outline,
                color="#999999",
                line_width=1,
                lighting=False,
                show_scalar_bar=False,
                name="tissue-footprint-outline",
            )
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
            p.add_mesh(lines, color=ink, opacity=0.22, line_width=1, show_scalar_bar=False)
    p.add_text(
        landscape.label if title is None else title,
        position=(0.045, 0.91),
        viewport=True,
        font_size=title_size,
        color=ink,
        font_file=font_file,
    )
    if axes:
        p.show_bounds(
            xtitle=f"x ({landscape.metadata['unit']})",
            ytitle=("−y" if image_coordinates else "y") + f" ({landscape.metadata['unit']})",
            ztitle="visual height",
            color=ink,
            grid=False,
        )
    p.view_isometric()
    p.enable_parallel_projection()
    p.camera.elevation = 12
    p.camera.zoom(1.25)
    p.enable_anti_aliasing("ssaa")
    if scalar_bar:
        p.subplot(1, 0)
        statistic = landscape.metadata.get("statistic", "mean")
        if scalar_bar_title is None:
            scalar_bar_title = landscape.label
            if statistic in {"mean", "sum"}:
                scalar_bar_title = f"{statistic.capitalize()} {landscape.label}"
                layer = landscape.metadata.get("layer")
                if layer:
                    scalar_bar_title += f" ({layer})"
        p.add_text(
            scalar_bar_title,
            position=(0.23, 0.75),
            viewport=True,
            font_size=label_size,
            font_file=font_file,
            color=ink,
        )
        bar = p.add_scalar_bar(
            title="",
            mapper=actor.mapper,
            color=ink,
            font_family="arial",
            label_font_size=label_size,
            bold=False,
            italic=False,
            shadow=False,
            fmt="%.3g",
            n_labels=5,
            position_x=0.23,
            position_y=0.13,
            width=0.54,
            height=0.42,
            vertical=False,
        )
        # PyVista exposes no font-file keyword on scalar bars. Use its TextProperty
        # wrapper to set the same portable font used by the plot and legend labels.
        bar.SetLabelTextProperty(
            pv.TextProperty(
                font_file=font_file, font_size=label_size, color=ink, bold=False, italic=False
            )
        )
        p.subplot(0, 0)
    if show:
        p.show()
    return p


def save(landscape, path, *, frames=360, fps=15, transparent_background=False, **kwargs):
    """Write an export and a JSON sidecar recording scientific and visual settings."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in {".npz", ".vtp", ".png", ".html", ".gif"}:
        raise ValueError("Choose .png, .html, .gif, .vtp or .npz.")
    if suffix == ".gif" and (
        isinstance(frames, bool)
        or not isinstance(frames, (int, np.integer))
        or frames < 2
        or not np.isfinite(fps)
        or fps <= 0
    ):
        raise ValueError("GIF exports need an integer frames >= 2 and finite fps > 0.")
    if transparent_background and suffix != ".png":
        raise ValueError("Transparent backgrounds are supported for PNG exports only.")
    path.parent.mkdir(parents=True, exist_ok=True)
    if suffix == ".npz":
        np.savez_compressed(
            path,
            x=landscape.x,
            y=landscape.y,
            values=landscape.values,
            mask=landscape.mask,
            counts=landscape.counts,
            tissue_mask=landscape.mask if landscape.tissue_mask is None else landscape.tissue_mask,
        )
    elif suffix == ".vtp":
        landscape.mesh(**kwargs).save(path)
    else:
        p = plot(landscape, off_screen=True, show=False, **kwargs)
        try:
            if suffix == ".png":
                p.screenshot(path, transparent_background=transparent_background)
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
        if suffix == ".png":
            render_options["transparent_background"] = transparent_background
    metadata = dict(
        landscape.metadata,
        label=landscape.label,
        render=render_options,
        scalps_version="0.1.0",
        value_range=[float(np.nanmin(landscape.values)), float(np.nanmax(landscape.values))],
    )
    if suffix == ".gif":
        metadata.update(frames=frames, fps=fps, orbit_seconds=frames / fps)
    path.with_suffix(path.suffix + ".json").write_text(
        json.dumps(metadata, indent=2, default=str) + "\n"
    )
    return path
