"""Render the documentation gallery from synthetic data only."""

from pathlib import Path

import scalps as sca

out = Path(__file__).resolve().parents[1] / "docs" / "assets"
adata = sca.demo()
for preset, value, title in [
    ("ember", "MKI67", "MKI67 · Synthetic tissue"),
    ("alpine", None, "Cell density · Synthetic tissue"),
    ("glacier", "CD3D", "CD3D · Synthetic tissue"),
]:
    mountain = sca.terrain(adata, value, resolution=300, smooth=3)
    mountain.save(out / f"{preset}.png", preset=preset, title=title)

mountain = sca.terrain(adata, "MKI67", resolution=160, smooth=2)
mountain.save(out / "orbit.gif", preset="ember", frames=48, fps=16, window_size=(700, 500))
mountain.save("outputs/demo.html", preset="ember")
