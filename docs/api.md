# API

## `scalps.terrain(data, value=None, **options) -> Terrain`

| Option | Default | Meaning |
| --- | --- | --- |
| `data` | required | AnnData or SpatialData |
| `value` | `None` | Gene, obs column, array; omitted for density |
| `statistic` | inferred | `density`, `mean`, `sum`, `fraction` |
| `groupby`, `groups` | `None` | obs column and selected category/categories |
| `layer` | `None` | Expression layer; default uses X |
| `spatial_key` | `"spatial"` | obsm coordinate key |
| `table`, `element` | `None` | SpatialData table and spatial element |
| `coordinate_system` | `"global"` | Target system for explicit element centroids |
| `resolution` | `250` | Longest unpadded grid dimension, 16–2000 |
| `smooth` | `2.0` | Gaussian sigma in pixels; zero disables smoothing |
| `support` | `0.02` | Minimum smoothed occupancy |
| `unit` | `"µm"` | Coordinate label; does not convert coordinates |
| `density_scale` | `1_000_000` | Multiplier for cells per coordinate-unit² |

`Terrain` exposes `x`, `y`, `values`, `mask`, `counts`, `label`, and `metadata`.
Arrays use `(x, y)` indexing. `counts` is the unsmoothed all-cell histogram.
Input objects are never modified.

## `scalps.plot(data, value=None, **options) -> pyvista.Plotter`

Accepts the above data options or a prebuilt `Terrain`.

| Option | Default | Meaning |
| --- | --- | --- |
| `preset` | `"alpine"` | `alpine`, `ember`, `glacier` |
| `height` | `0.22` | Maximum relief relative to longest padded extent |
| `vmax` | automatic | Absolute value mapping to maximum relief |
| `transform` | `"linear"` | `linear`, `sqrt`, `log1p`, signed height mapping |
| `clim`, `cmap` | automatic | Raw color limits and colormap |
| `contours` | `12` | Contour count; 0 disables |
| `skirt` | `True` | Draw decorative cutaway sides |
| `title` | value label | Figure title |
| `axes` | `False` | Show coordinate bounds |
| `scalar_bar` | `True` | Show legend |
| `image_coordinates` | `True` | Reverse y for image convention |
| `off_screen` | PyVista default | Headless rendering |
| `window_size` | `(1400, 1000)` | Render size in pixels |
| `show` | `True` | Display immediately; false returns open plotter |

`Terrain.plot()` accepts the visual options. `Terrain.mesh(height=0.22,
vmax=None, transform="linear")` returns the terrain surface.

## `Terrain.save(path, **options) -> pathlib.Path`

PNG, HTML and GIF accept visual options except `show`/`off_screen` (set internally).
GIF adds `frames=90, fps=24`. HTML requires the notebook extra.
VTP accepts only `mesh()` options. NPZ stores the arrays directly.
Files are overwritten when their path already exists. A settings JSON sidecar
is written next to every successful export.

## `scalps.demo(n_cells=30_000, seed=7) -> AnnData`

Synthetic lobed tissue with sparse genes `MKI67`, `CD3D`, `EPCAM`, cell types
`Tumor`, `T cell`, `Stromal`, and scores `prolif_score`, `immune_balance`.

## CLI

`scalps INPUT [--value KEY] [--groupby COLUMN --groups GROUP ...]
[--layer LAYER] [--resolution N] [--smooth SIGMA] [--preset PRESET]
[-o OUTPUT]`

INPUT is an `.h5ad` file or `demo`. Without output, an interactive plot opens.
