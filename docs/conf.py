"""Documentation built with the same Sphinx Book Theme as cellpin."""

import os
import tomllib
from pathlib import Path

info = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())["project"]
project = "scAlps"
author = info["authors"][0]["name"]
copyright = f"2026, {author}"
release = version = info["version"]

# Read package metadata without importing scAlps or installing rendering libraries.
extensions = ["myst_parser", "sphinx_copybutton"]
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
nitpicky = True
myst_heading_anchors = 6

html_theme = "sphinx_book_theme"
html_title = project
html_logo = "assets/scalps-logo.png"
html_favicon = "assets/scalps-icon.png"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "https://scalps.readthedocs.io/en/latest/")
html_context = {"default_mode": "light"}
html_theme_options = {
    "repository_url": "https://github.com/philinscience/scAlps",
    "repository_branch": "main",
    "use_repository_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": False,
    "home_page_in_toc": True,
}
pygments_style = "default"
