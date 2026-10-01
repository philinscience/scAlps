# Developing scAlps

```bash
uv venv --python 3.11
uv pip install -e '.[dev,notebook,spatial]'
.venv/bin/pytest -q
.venv/bin/ruff check .
.venv/bin/mkdocs build --strict
.venv/bin/python -m build
```

To work against the local PyVista checkout, install `-e ../pyvista` into the
same environment. The normal dependency uses released PyVista for portability.

The numerical core lives in `core.py`; input adapters in `_input.py`; rendering
and exports in `render.py`. Keep scientific aggregation separate from styling.
Do not normalize expression implicitly or mutate input objects. Test weighted
aggregation, missing/zero/signed values, sparse inputs and instance alignment
when changing those paths.

`examples/gallery.py` regenerates documentation images from seeded synthetic
data. Keep real datasets in ignored `data/`, and local exports in `outputs/`.
The CI workflow tests and builds docs but does not deploy or publish anything.
