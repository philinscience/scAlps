import json

import numpy as np
import pytest
from PIL import Image

import scalps as sca


def test_numerical_export(tmp_path):
    t = sca.terrain(sca.demo(1000), "MKI67", resolution=30)
    target = t.save(tmp_path / "terrain.npz")
    with np.load(target) as loaded:
        np.testing.assert_allclose(loaded["values"], t.values)
    meta = json.loads(target.with_suffix(".npz.json").read_text())
    assert meta["n_cells"] == 1000
    assert t.save(tmp_path / "terrain.vtp").stat().st_size > 1000


@pytest.mark.parametrize("preset", ["alpine", "ember", "glacier"])
def test_png(tmp_path, preset):
    t = sca.terrain(sca.demo(2000), "MKI67", resolution=40)
    path = t.save(tmp_path / f"{preset}.png", preset=preset, window_size=(400, 300))
    image = np.array(Image.open(path))
    assert image.shape[:2] == (300, 400)
    assert image.std() > 10
    np.testing.assert_array_equal(image[0, 0, :3], [255, 255, 255])


def test_transparent_png(tmp_path):
    t = sca.terrain(sca.demo(1000), "MKI67", resolution=30)
    path = t.save(tmp_path / "transparent.png", transparent_background=True, window_size=(400, 300))
    image = np.array(Image.open(path))
    assert image.shape == (300, 400, 4)
    assert image[0, 0, 3] == 0
    assert image[:, :, 3].max() == 255
    metadata = json.loads(path.with_suffix(".png.json").read_text())
    assert metadata["render"]["transparent_background"] is True


def test_html_gif(tmp_path):
    pytest.importorskip("trame")
    t = sca.terrain(sca.demo(1000), "MKI67", resolution=30)
    html = t.save(tmp_path / "terrain.html", window_size=(320, 240))
    assert "<html" in html.read_text().lower()
    gif = t.save(tmp_path / "terrain.gif", frames=4, fps=4, window_size=(320, 240))
    with Image.open(gif) as im:
        assert im.n_frames == 4
