"""Input adapters; never densify an entire expression matrix."""

import numpy as np
from anndata import AnnData
from scipy import sparse


def resolve(data, spatial_key, table, element, coordinate_system):
    if isinstance(data, AnnData):
        if element is not None or table is not None:
            raise ValueError("table and element are only used with SpatialData.")
        adata = data
        if spatial_key not in adata.obsm:
            raise KeyError(f"Missing adata.obsm[{spatial_key!r}].")
        xy = np.asarray(adata.obsm[spatial_key], dtype=float)
    else:
        try:
            from spatialdata import SpatialData, get_centroids
        except ImportError as exc:
            raise ImportError("SpatialData input needs: pip install 'scalps[spatial]'") from exc
        if not isinstance(data, SpatialData):
            raise TypeError("Expected an AnnData or SpatialData object.")
        if table is None:
            if len(data.tables) != 1:
                raise ValueError(
                    "Choose table= explicitly when SpatialData has multiple/no tables."
                )
            table = next(iter(data.tables))
        adata = data.tables[table]
        if element is None:
            if spatial_key not in adata.obsm:
                raise ValueError(
                    "Specify element= to get coordinates from shapes, labels or points."
                )
            xy = np.asarray(adata.obsm[spatial_key], dtype=float)
        else:
            attrs = adata.uns.get("spatialdata_attrs", {})
            region_key, instance_key = attrs.get("region_key"), attrs.get("instance_key")
            if region_key is None or instance_key is None:
                raise ValueError("The table must annotate the element with region/instance keys.")
            adata = adata[adata.obs[region_key].astype(str).eq(element)]
            centroids = get_centroids(data[element], coordinate_system=coordinate_system).compute()
            ids = adata.obs[instance_key]
            if ids.duplicated().any() or not centroids.index.is_unique:
                raise ValueError("Instance IDs must be unique within the selected element.")
            if not ids.isin(centroids.index).all():
                raise ValueError("Some table instance IDs have no matching spatial centroid.")
            xy = centroids.loc[ids, ["x", "y"]].to_numpy(dtype=float)
    if xy.ndim != 2 or xy.shape != (adata.n_obs, 2):
        raise ValueError("Spatial coordinates must have shape (n_obs, 2); select a 2D plane first.")
    if not len(xy) or not np.isfinite(xy).all():
        raise ValueError("Spatial coordinates must be nonempty and finite.")
    if np.any(np.ptp(xy, axis=0) <= 0):
        raise ValueError("Coordinates must span both x and y axes.")
    return adata, xy


def values(adata, key, layer):
    if key is None:
        if layer is not None:
            raise ValueError("layer= requires a gene value.")
        return None, "Cell density"
    if not isinstance(key, str):
        array = np.asarray(key, dtype=float)
        if array.shape != (adata.n_obs,):
            raise ValueError("An explicit value array must have one entry per selected table row.")
        if layer is not None:
            raise ValueError("layer= requires a gene name.")
        return array, "Value"
    source, sep, name = key.partition(":")
    if not sep or source not in {"obs", "gene"}:
        source, name = "auto", key
    in_obs, in_var = name in adata.obs, name in adata.var_names
    if source == "auto" and in_obs and in_var:
        raise ValueError(f"Ambiguous {name!r}; use 'obs:{name}' or 'gene:{name}'.")
    if source == "obs" or (source == "auto" and in_obs):
        if layer is not None:
            raise ValueError("layer= only applies to genes.")
        try:
            return adata.obs[name].to_numpy(dtype=float), name
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Values must be numeric; use groupby= and groups= for cell types."
            ) from exc
    if not in_var:
        raise KeyError(f"No observation column or gene named {name!r}.")
    if not adata.var_names.is_unique:
        raise ValueError("Gene names must be unique; call adata.var_names_make_unique().")
    # Slice before materializing, including backed sparse datasets.
    matrix = adata[:, [name]].layers[layer] if layer is not None else adata[:, [name]].X
    if matrix is None:
        raise ValueError("The requested expression matrix is empty.")
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float).reshape(-1), name
