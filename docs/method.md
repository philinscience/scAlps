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

The footprint is `G(C > 0) >= support`. Unsupported values are stored as NaN;
quads touching unsupported vertices are removed from the surface. This is a
sampling-based approximation, not a segmentation mask. Tiny holes and nearby
islands may merge with smoothing. Large gaps remain empty. Mean fields also
require local valid-value support; selected-cell means can have a smaller footprint.

Changing resolution, smoothing, or support changes this approximation. A low
cell count or an excessively fine grid may need a lower `support`. Density
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

Skirts are decorative vertical edges; they are not a watertight solid model.
Mesh export preserves the original x/y coordinates. Plotting defaults to
image-style y direction, which is recorded with other render settings.

These are exploratory visualizations, not spatial hypothesis tests. No
expression normalization, batch correction, uncertainty estimation, or
biological inference is performed.
