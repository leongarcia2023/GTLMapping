from __future__ import annotations

from datetime import UTC, datetime
from importlib.metadata import version as package_version

project = "GTLMapping"
author = "León Garcia"
current_year = datetime.now(UTC).year
copyright = (
    f"2026, {author}"
    if current_year == 2026
    else f"2026-{current_year}, {author}"
)
release = package_version("GTLMapping")
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]
autosummary_generate = True
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "GTLMapping"
html_context = {
    "display_github": True,
    "github_user": "leongarcia2023",
    "github_repo": "GTLMapping",
    "github_version": "main",
    "conf_py_path": "/docs/",
}
