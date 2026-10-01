# Developing scAlps

```bash
uv venv --python 3.11
uv pip install -e '.[dev,notebook,spatial]'
uv pip install -r docs/requirements.txt
.venv/bin/pytest -q
.venv/bin/ruff check .
.venv/bin/sphinx-build -b dirhtml -n -W --keep-going docs site
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
GitHub Actions tests Python 3.11 and 3.12, checks dependency consistency, and
builds the package. A separate documentation job uses the same requirements
as Read the Docs, without installing scAlps or its rendering dependencies.

The `spatial` extra currently constrains AnnData to `<0.13` for compatibility
with SpatialData `<0.8`. Upgrade and validate these dependencies together.

## Read the Docs

The live documentation is at [scalps.readthedocs.io](https://scalps.readthedocs.io/en/latest/).
The README badge reports the build status of the `latest` version.

The repository includes `.readthedocs.yaml` for Sphinx on Python 3.12, with
strict warning checks and pinned documentation tools in `docs/requirements.txt`.
It uses Sphinx Book Theme, matching cellpin, with MyST parsing the existing
Markdown pages. The `dirhtml` builder preserves URLs such as `/api/` and `/recipes/`.
The gallery is committed synthetic data; documentation builds do not load
private slides or regenerate plots.

To build the docs alone:

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -b dirhtml -n -W --keep-going docs site
python -m http.server --directory site
```

For the one-time hosting setup, [import the repository into Read the Docs](https://app.readthedocs.org/dashboard/import/).
Select `philinscience/scAlps`, use `main` as the default branch, and keep the
configuration path `.readthedocs.yaml`. Trigger the first `latest` build and
verify the GitHub integration/webhook so subsequent pushes rebuild the docs.
For a fork, update the documentation URL and badge project slug in the README
and the GitHub repository's website field. The site takes its canonical URL from
`READTHEDOCS_CANONICAL_URL`, so it also works if the assigned project slug differs.

See the [Read the Docs Sphinx guide](https://docs.readthedocs.com/platform/stable/intro/sphinx.html)
for the hosting integration. GitHub Actions only validates builds; Read the
Docs publishes the documentation after the project has been imported.
