import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy.sparse import csr_matrix

import scalps as sca


def tissue():
    x, y = np.meshgrid(np.arange(12), np.arange(8))
    xy = np.column_stack([x.ravel(), y.ravel()]).astype(float)
    a = AnnData(
        csr_matrix(np.full((len(xy), 2), 3.0)),
        obs=pd.DataFrame(
            {"score": np.full(len(xy), 4.0), "type": np.where(xy[:, 0] < 6, "A", "B")},
            index=[str(i) for i in range(len(xy))],
        ),
        var=pd.DataFrame(index=["G", "H"]),
    )
    a.obsm["spatial"] = xy
    a.layers["counts"] = a.X * 2
    return a


def test_constant_mean_and_sparse_layer():
    a = tissue()
    for value, layer, expected in [("score", None, 4), ("G", None, 3), ("G", "counts", 6)]:
        t = sca.terrain(a, value, layer=layer, resolution=30)
        np.testing.assert_allclose(t.values[t.mask], expected)
    assert isinstance(a.X, csr_matrix)
    assert set(a.obsm) == {"spatial"}


def test_cell_weighted_mean():
    a = tissue()
    xy = a.obsm["spatial"]
    # Put unequal numbers of cells at two coordinates that land in the same bin.
    a = AnnData(
        np.zeros((102, 1)),
        obs=pd.DataFrame({"v": [10.0] * 100 + [0.0, np.nan]}, index=[str(i) for i in range(102)]),
    )
    a.obsm["spatial"] = np.vstack([np.tile(xy[0], (101, 1)), [100, 100]])
    t = sca.terrain(a, "v", resolution=16, smooth=1)
    np.testing.assert_allclose(t.values[t.mask], 1000 / 101)


def test_grouped_mean_with_pandas_copy_on_write():
    # pandas 3 enables CoW by default; enable it on pandas 2 as well to catch
    # attempts to modify read-only arrays returned by Series.to_numpy().
    with pd.option_context("mode.copy_on_write", True):
        a = tissue()
        a.obs.loc["0", "score"] = np.nan
        original = a.obs.copy(deep=True)
        t = sca.terrain(a, "score", groupby="type", groups="A", color="G", resolution=30)
        np.testing.assert_allclose(t.values[t.mask], 4)
        np.testing.assert_allclose(t.color_values[t.color_mask], 3)
        assert t.metadata["n_selected"] == 47
        assert t.metadata["color"]["n_selected"] == 48
        pd.testing.assert_frame_equal(a.obs, original)


def test_density_units_and_fraction():
    a = tissue()
    total = sca.terrain(a, resolution=30, support=1e-8)
    selected = sca.terrain(a, groupby="type", groups="A", resolution=30, support=1e-8)
    frac = sca.terrain(
        a, groupby="type", groups="A", statistic="fraction", resolution=30, support=1e-8
    )
    np.testing.assert_array_equal(total.mask, selected.mask)
    np.testing.assert_allclose(
        frac.values[frac.mask], selected.values[frac.mask] / total.values[frac.mask]
    )
    a.obsm["spatial"] *= 2
    bigger = sca.terrain(a, resolution=30, support=1e-8)
    np.testing.assert_allclose(bigger.values, total.values / 4)


def test_signed_zero_nan_and_disconnected():
    a = tissue()
    a.obs["score"] = -2.0
    t = sca.terrain(a, "score", resolution=30)
    assert t.mesh().points[:, 2].max() < 0
    a.obs["score"] = 0.0
    a.obs.iloc[0, a.obs.columns.get_loc("score")] = np.nan
    t = sca.terrain(a, "score", resolution=30)
    assert np.isfinite(t.mesh().points).all()
    assert np.all(t.values[t.mask] == 0)
    a.obsm["spatial"][a.obsm["spatial"][:, 0] > 5, 0] += 100
    t = sca.terrain(a, resolution=100, smooth=1)
    middle = np.argmin(abs(t.x - 50))
    assert not t.mask[middle].any()


def test_orientation_and_bin_centers():
    a = tissue()
    a.obs["score"] = a.obsm["spatial"][:, 0]
    t = sca.terrain(a, "score", resolution=30, smooth=1)
    mesh = t.mesh()
    assert np.corrcoef(mesh.points[:, 0], mesh["value"])[0, 1] > 0.98
    assert np.isclose(np.diff(t.x)[0], np.diff(t.y)[0])


def test_backed_sparse(tmp_path):
    path = tmp_path / "a.h5ad"
    tissue().write_h5ad(path)
    a = read_h5ad(path, backed="r")
    try:
        t = sca.terrain(a, "G", color="H", resolution=30)
        np.testing.assert_allclose(t.values[t.mask], 3)
        np.testing.assert_allclose(t.color_values[t.color_mask], 3)
    finally:
        a.file.close()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"resolution": 0},
        {"resolution": 2.5},
        {"smooth": -1},
        {"support": 0},
        {"groupby": "type"},
        {"groups": "A"},
        {"value": "nope"},
        {"statistic": "fraction"},
        {"value": "G", "statistic": "density"},
        {"groupby": "type", "groups": "missing"},
    ],
)
def test_validation(kwargs):
    with pytest.raises((ValueError, KeyError)):
        sca.terrain(tissue(), **kwargs)


def test_ambiguous_names():
    a = tissue()
    a.obs["G"] = 9
    with pytest.raises(ValueError, match="Ambiguous"):
        sca.terrain(a, "G")
    t = sca.terrain(a, "obs:G", resolution=30)
    np.testing.assert_allclose(t.values[t.mask], 9)


def test_demo_reproducible():
    a, b = sca.demo(1000), sca.demo(1000)
    np.testing.assert_array_equal(a.obsm["spatial"], b.obsm["spatial"])
    assert (a.X != b.X).nnz == 0


def test_sparse_outskirts_trimmed_without_filtering_rare_types_or_zero_genes():
    rng = np.random.default_rng(31)
    xy = np.vstack([rng.normal(size=(3000, 2)), [[12.0, 12.0], [-12.0, -12.0]]])
    a = AnnData(np.zeros((len(xy), 1)), var=pd.DataFrame(index=["zero"]))
    a.obsm["spatial"] = xy
    a.obs["type"] = np.where(np.arange(len(xy)) == 0, "rare", "other")
    untrimmed = sca.terrain(a, resolution=100, density_percentile=0, support=0.001)
    trimmed = sca.terrain(a, resolution=100, support=0.001)
    rare = sca.terrain(a, groupby="type", groups="rare", resolution=100, support=0.001)
    zero = sca.terrain(a, "zero", resolution=100, support=0.001)
    edge = (np.argmin(abs(trimmed.x - 12)), np.argmin(abs(trimmed.y - 12)))
    assert untrimmed.mask[edge] and not trimmed.mask[edge]
    assert 0 < trimmed.metadata["n_low_density_cells"] <= 31
    np.testing.assert_array_equal(trimmed.mask, rare.mask)
    np.testing.assert_array_equal(trimmed.tissue_mask, zero.tissue_mask)
    assert np.max(rare.values[rare.mask]) > 0
    assert np.all(zero.values[zero.mask] == 0)
    # Trimming only changes display support, not the estimated field in retained tissue.
    np.testing.assert_allclose(trimmed.values[trimmed.mask], untrimmed.values[trimmed.mask])


def test_footprint_retains_all_tissue_for_selected_cell_mean():
    a = tissue()
    t = sca.terrain(a, "score", groupby="type", groups="A", resolution=40, smooth=2)
    assert t.tissue_mask.sum() > t.mask.sum()
    floor = t.footprint(z=-3)
    np.testing.assert_array_equal(floor.points[:, 2], -3)
    assert floor.bounds[1] > t.mesh().bounds[1]


@pytest.mark.parametrize("percentile", [-1, 100, np.nan])
def test_invalid_density_percentile(percentile):
    with pytest.raises(ValueError, match="density_percentile"):
        sca.terrain(tissue(), density_percentile=percentile)
