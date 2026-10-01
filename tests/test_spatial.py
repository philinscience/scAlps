import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scalps as sca

sd = pytest.importorskip("spatialdata")
from spatialdata.models import PointsModel, TableModel  # noqa: E402
from spatialdata.transformations import Identity, Translation  # noqa: E402


def test_centroids_align_by_instance_and_transform():
    points = PointsModel.parse(
        pd.DataFrame(
            {"x": [0.0, 10.0, 20.0, 30.0], "y": [0.0, 20.0, 10.0, 30.0]}, index=[40, 20, 10, 30]
        ),
        coordinates={"x": "x", "y": "y"},
        transformations={"global": Translation([100.0, 200.0], axes=("x", "y"))},
    )
    table = AnnData(
        np.array([[1.0], [2.0], [3.0], [4.0]]),
        obs=pd.DataFrame(
            {"region": pd.Categorical(["cells"] * 4), "instance": [10, 30, 40, 20]},
            index=["a", "b", "c", "d"],
        ),
        var=pd.DataFrame(index=["G"]),
    )
    table = TableModel.parse(table, region="cells", region_key="region", instance_key="instance")
    sdata = sd.SpatialData(points={"cells": points}, tables={"table": table})
    actual = sca.terrain(sdata, "G", element="cells", resolution=16, support=0.001)
    expected = table.copy()
    expected.obsm["spatial"] = np.array([[120, 210], [130, 230], [100, 200], [110, 220]])
    reference = sca.terrain(expected, "G", resolution=16, support=0.001)
    np.testing.assert_allclose(actual.values, reference.values)
    np.testing.assert_allclose(actual.x, reference.x)
    assert "spatial" not in table.obsm


def test_table_obsm_and_multiple_tables():
    a = sca.demo(1000)
    sdata = sd.SpatialData(tables={"a": a, "b": a.copy()})
    with pytest.raises(ValueError, match="table="):
        sca.terrain(sdata)
    result = sca.terrain(sdata, table="a", resolution=30)
    reference = sca.terrain(a, resolution=30)
    np.testing.assert_allclose(result.values, reference.values)


def test_shapes_centroids():
    import geopandas as gpd
    from shapely.geometry import Point
    from spatialdata.models import ShapesModel

    shapes = ShapesModel.parse(
        gpd.GeoDataFrame(
            {"geometry": [Point(0, 0), Point(10, 10), Point(20, 0)], "radius": [1.0, 1.0, 1.0]},
            index=[2, 0, 1],
        ),
        transformations={"global": Identity()},
    )
    table = AnnData(
        np.ones((3, 1)),
        obs=pd.DataFrame(
            {"region": pd.Categorical(["cells"] * 3), "instance": [0, 1, 2]}, index=["a", "b", "c"]
        ),
        var=pd.DataFrame(index=["G"]),
    )
    table = TableModel.parse(table, region="cells", region_key="region", instance_key="instance")
    sdata = sd.SpatialData(shapes={"cells": shapes}, tables={"table": table})
    t = sca.terrain(sdata, "G", element="cells", resolution=16, support=0.001)
    np.testing.assert_allclose(t.values[t.mask], 1)
