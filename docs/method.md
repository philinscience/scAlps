# Method and interpretation

## Aggregation

Cells are binned on a square grid using bin centers, preserving spatial aspect
ratio. `resolution` controls the longest unpadded grid dimension. Padding adds
three Gaussian standard deviations around the tissue, plus one pixel.

Let `G` denote a Gaussian filter with constant-zero padding. For cell counts
`C`, selected finite-value counts `N`, and selected value sums `S`:

| Statistic | Field |
| --- | --- |
| density | `G(N) / pixel_area * density_scale` |
| mean | `G(S) / G(N)` where the denominator is positive |
| sum | `G(S)` (smoothed total per grid bin) |
| fraction | `G(N) / G(C)` |

For density/fraction, `N` includes the selected cells, without value filtering.
For means and sums, nonfinite values are excluded. Mean calculation weights
individual cells, rather than giving equally weighted votes to bins with
different cell counts. Zero measurements are valid and retained.

## Tissue support

The initial footprint is `G(C > 0) >= support`. By default, it is intersected
with `G(C) >= cutoff`, where `cutoff` is the first percentile of `G(C)` sampled
at the input cell locations. Each cell contributes one sample, including cells
sharing a bin. This avoids letting empty grid area dominate the percentile.
Use `density_percentile=0` to disable this step, or choose another percentile.

The cutoff depends on all cells, independently of genes and selected groups.
It trims display support; it does not delete cells from the input or recompute
the smoothed values after exclusion. The fraction of removed grid area need
not be 1%. Density ties may also mean fewer than 1% of cells fall strictly below
the cutoff. True low-density tissue can be hidden too; this is not cell QC.

Unsupported values are stored as NaN;
quads touching unsupported vertices are removed from the surface. This is a
sampling-based approximation, not a segmentation mask. Tiny holes and nearby
islands may merge with smoothing. Large gaps remain empty. Mean fields also
require local valid-value support; selected-cell means can have a smaller footprint.

Changing resolution, smoothing, support or the density percentile changes this
approximation. A low cell count or an excessively fine grid may need a lower
`support`, a lower `density_percentile`, or more smoothing. Density
does not divide by occupancy; edge density can decrease because there are no
cells outside tissue. Exported masked values are not a mass-conserving estimate.

## Units

Coordinates are not rescaled. The default `unit="µm", density_scale=1_000_000`
reports cells per square millimeter. For coordinates already in millimeters,
use `unit="mm", density_scale=1`. For pixels, choose an appropriate physical
conversion upstream, or use `unit="px", density_scale=1` for cells/pixel².

Smoothing bandwidth is `smooth * metadata["spacing"]` in coordinate units.

## Height and colors

Geometric height is the signed value divided by its 99th percentile absolute
scale, clipped to [-1, 1], optionally transformed, then multiplied by
`height * longest_padded_extent`. Zero maps to the zero plane. Negative values
produce valleys; constant zeros remain flat. `vmax` overrides the scale.
Colors always encode raw aggregated values. `clim` affects color scaling only.

By default, color uses the same field as height. With `color=...`, a second
field is computed on the identical grid and among the same selected groups.
Its default mean is `G(S_color) / G(N_color)`; `color_statistic="sum"` uses
`G(S_color)`. `N_color` counts only finite color measurements, independently
of missing height measurements. Color support requires `G(N_color) > 1e-12`
within height support. Undefined color is NaN and rendered gray; it does not
change the height field or remove geometry. Real zeros contribute to the mean.
The finite Gaussian neighborhood can include very few cells; no minimum-cell
confidence filter or uncertainty estimate is implied by a visible color.
Contours always trace height values, and the color bar describes only color.
Height and color can encode different biological quantities, but neither adds
information beyond the supplied data and chosen spatial aggregation.

Skirts are decorative vertical edges; they are not a watertight solid model.
The optional gray reference layer is a flat copy of the cleaned all-cell tissue
footprint, positioned below the lowest terrain/skirt point. It uses the same
x/y coordinates and preserves holes. For selected-cell means, it can extend
beyond the colored terrain's valid-value support. Its vertical offset is for
visual separation, not a physical distance. VTP export contains only the
colored surface; call `Terrain.footprint()` to obtain the flat reference mesh.
Mesh export preserves the original x/y coordinates. Plotting defaults to
image-style y direction, which is recorded with other render settings.

These are exploratory visualizations, not spatial hypothesis tests. No
expression normalization, batch correction, uncertainty estimation, or
biological inference is performed.
