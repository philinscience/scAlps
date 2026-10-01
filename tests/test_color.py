"""Independent encodings must share coordinates, but never missing-value filters."""

import json

import numpy as np
import pytest
import pyvista as pv
from test_core import tissue

import scalps as sca


def test_color_does_not_change_density_or_geometry():
    a = tissue()
    a.obs["activation"] = np.where(a.obs["type"] == "A", a.obsm["spatial"][:, 1], 9999)
    a.obs.iloc[0, a.obs.columns.get_loc("activation")] = np.nan
    opts = dict(groupby="type", groups="A", resolution=40, smooth=2)
    height = sca.terrain(a, **opts)
    dual = sca.terrain(a, color="activation", **opts)
    reference = sca.terrain(a, "activation", **opts)
    np.testing.assert_array_equal(dual.mask, height.mask)
    np.testing.assert_allclose(dual.values, height.values)
    np.testing.assert_allclose(dual.mesh().points, height.mesh().points)
    np.testing.assert_array_equal(dual.color_mask, reference.mask)
    np.testing.assert_allclose(dual.color_values, reference.values)
    assert dual.metadata["n_selected"] == 48
    assert dual.metadata["color"]["n_selected"] == 47
    assert np.nanmax(dual.color_values) <= 7  # other cells never enter the mean
    assert np.any(dual.mask & ~dual.color_mask)


def test_missing_height_does_not_remove_color_measurements():
    a = tissue()
    a.obs.loc["0", "score"] = np.nan
    a.obs["color"] = 0.0
    a.obs.loc["0", "color"] = 100.0
    dual = sca.terrain(a, "score", color="color", resolution=30)
    ref = sca.terrain(a, "color", resolution=30)
    np.testing.assert_allclose(dual.color_values[dual.color_mask], ref.values[dual.color_mask])
    assert dual.metadata["color"]["n_selected"] == len(a)
    assert dual.metadata["n_selected"] == len(a) - 1


@pytest.mark.parametrize("statistic", ["mean", "sum"])
def test_sparse_color_layer_and_zero_values(statistic):
    a = tissue()
    t = sca.terrain(
        a, "score", color="gene:G", color_layer="counts", color_statistic=statistic, resolution=30
    )
    ref = sca.terrain(a, "G", layer="counts", statistic=statistic, resolution=30)
    np.testing.assert_allclose(t.color_values[t.color_mask], ref.values[t.color_mask])
    zero = sca.terrain(a, color=np.zeros(len(a)), resolution=30)
    assert np.all(zero.color_values[zero.color_mask] == 0)
    assert zero.color_mask.any()
    assert t.metadata["layer"] is None
    assert t.metadata["color"]["layer"] == "counts"


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"color": np.full(96, np.nan)}, "No finite color"),
        ({"color": "score", "color_statistic": "density"}, "color_statistic"),
        ({"color_layer": "counts"}, "require color"),
        ({"color_statistic": "sum"}, "require color"),
    ],
)
def test_invalid_color(kwargs, match):
    with pytest.raises(ValueError, match=match):
        sca.terrain(tissue(), **kwargs)


def test_both_fields_exported_without_pickling(tmp_path):
    t = sca.terrain(tissue(), color="G", color_layer="counts", resolution=30)
    path = t.save(tmp_path / "both.npz")
    with np.load(path, allow_pickle=False) as arrays:
        np.testing.assert_allclose(arrays["values"], t.values)
        np.testing.assert_allclose(arrays["color_values"], t.color_values)
        np.testing.assert_array_equal(arrays["color_mask"], t.color_mask)
    meta = json.loads(path.with_suffix(".npz.json").read_text())
    assert meta["color"]["statistic"] == "mean"
    assert meta["color"]["layer"] == "counts"
    np.testing.assert_allclose(meta["color_range"], [6, 6])
    mesh = pv.read(t.save(tmp_path / "both.vtp"))
    np.testing.assert_allclose(mesh["color"], 6)
    np.testing.assert_allclose(mesh["value"], t.mesh()["value"])


def test_renderer_uses_color_limits_and_keeps_height():
    a = tissue()
    a.obs["score"] = -2.0
    t = sca.terrain(a, color="score", resolution=30)
    p = t.plot(off_screen=True, show=False)
    try:
        actor = p.actors["terrain"]
        assert actor.mapper.array_name == "color"
        np.testing.assert_allclose(actor.mapper.scalar_range, [-2, 2])
        assert actor.mapper.dataset.points[:, 2].min() >= 0
    finally:
        p.close()
