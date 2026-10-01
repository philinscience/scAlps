"""Reproducible synthetic tissue; no downloads or patient data."""

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix


def demo(n_cells=30_000, seed=7):
    """Create a lobed synthetic tissue with cell types, gene counts and signed scores."""
    if n_cells < 100:
        raise ValueError("Use at least 100 cells for the demo.")
    rng = np.random.default_rng(seed)
    centers = np.array(
        [
            [950, 1100],
            [1700, 1900],
            [2650, 1300],
            [3300, 2300],
            [1800, 2900],
            [3000, 3400],
            [4000, 3100],
        ]
    )
    sizes = np.array(
        [[370, 500], [450, 400], [450, 350], [400, 440], [410, 400], [420, 400], [280, 310]]
    )
    component = rng.choice(len(centers), n_cells, p=[0.1, 0.18, 0.16, 0.22, 0.12, 0.16, 0.06])
    xy = centers[component] + rng.normal(size=(n_cells, 2)) * sizes[component]
    x, y = xy.T
    peak = lambda cx, cy, width: np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * width**2))  # noqa: E731
    proliferation = 3 * peak(3300, 2250, 300) + 2 * peak(1500, 1700, 250)
    immune = 2 * peak(2300, 2850, 350) + peak(3700, 3200, 220)
    expression = rng.poisson(
        np.column_stack(
            [0.1 + 5 * proliferation, 0.2 + 7 * immune, 0.2 + 4 * peak(2550, 1300, 400)]
        )
    ).astype(np.float32)
    cell_type = np.where(
        rng.random(n_cells) < np.clip(immune / 2, 0.08, 0.85),
        "T cell",
        np.where(component % 3 == 0, "Stromal", "Tumor"),
    )
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(cell_type),
            "prolif_score": proliferation + rng.normal(0, 0.08, n_cells),
            "immune_balance": immune - proliferation,
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    adata = AnnData(
        csr_matrix(expression), obs=obs, var=pd.DataFrame(index=["MKI67", "CD3D", "EPCAM"])
    )
    adata.obsm["spatial"] = xy
    adata.layers["counts"] = adata.X.copy()
    # Deliberately independent of the T-cell density field; an illustration only.
    activation = np.clip(0.1 + 0.8 * peak(3400, 2700, 700) + rng.normal(0, 0.06, n_cells), 0, 1)
    adata.obs["activation_score"] = np.where(cell_type == "T cell", activation, np.nan)
    adata.uns["scalps"] = {"synthetic": True, "seed": seed, "coordinate_unit": "µm"}
    return adata
