# Recipes

## Cell types and local fractions

```python
sca.plot(adata, groupby="cell_type", groups=["T cell", "B cell"])
sca.plot(
    adata,
    groupby="cell_type",
    groups="T cell",
    statistic="fraction",
    clim=(0, 1),
    vmax=1,
    preset="glacier",
)
```

Density uses only selected cells in its numerator; fraction divides by all
cells locally. Both retain the full tissue footprint.

## Expression and scores

```python
sca.plot(adata, "gene:MKI67", layer="counts", preset="ember")
sca.plot(adata, "obs:prolif_score", transform="sqrt", height=0.3)
sca.plot(adata, "immune_balance", cmap="RdBu_r")
```

The `gene:` and `obs:` prefixes disambiguate duplicate names. `sqrt` and `log1p`
change only the height mapping, preserving signs; they do not normalize expression.
Prepare normalized expression upstream and pass its layer if desired.

## SpatialData with transformed geometry

```python
sca.plot(
    sdata,
    "MKI67",
    table="table",
    element="cell_boundaries",
    coordinate_system="global",
    preset="ember",
)
```

The table must annotate this element using SpatialData's region and instance
metadata. Only rows in the selected element are used. Centroids are transformed
and joined by instance ID. Shapes, labels and points are supported by
[SpatialData's centroid API](https://spatialdata.scverse.org/en/latest/api/operations.html).
For arrays passed as `value`, supply one value per retained table row.

## Custom PyVista composition

```python
mountain = sca.terrain(adata, "MKI67")
p = mountain.plot(show=False, preset="ember", scalar_bar=False, contours=0)
p.add_axes()
p.camera.elevation = 35
p.show()
```

For complete control, `mountain.mesh()` returns a PyVista `PolyData` with
`value` and `support` point arrays in the original x/y coordinate system.
The plotter's default image orientation reverses y for Xenium. Set
`image_coordinates=False` for Cartesian data.

## Comparing slides

Use the same value preprocessing, coordinate units, `vmax`, `clim`, `height`
and transform. Automatic per-slide peak scaling does not convey comparable
absolute height. `smooth` is measured in grid pixels, so equal values do not
imply equal physical bandwidth across slides of different extent. Inspect
`mountain.metadata["spacing"]` and choose resolution/smooth accordingly.

## Large slides

```python
import anndata as ad

adata = ad.read_h5ad("slide.h5ad", backed="r")
try:
    mountain = sca.terrain(adata, "MKI67", resolution=300)
    mountain.save("outputs/slide.png", preset="ember")
finally:
    adata.file.close()
```

Coordinates and one gene vector are loaded; the full expression matrix is not.
Start at resolution 250–400. Doubling resolution roughly quadruples mesh size.
SpatialData label centroids may require scanning the whole segmentation image.
