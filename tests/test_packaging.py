"""Release-facing consistency checks."""

from __future__ import annotations

from importlib.metadata import version as installed_version
from pathlib import Path
import tomllib

import gtlmapping


def test_version_is_consistent_across_package_metadata() -> None:
    project_root = Path(__file__).resolve().parents[1]
    with (project_root / "pyproject.toml").open("rb") as handle:
        configured_version = tomllib.load(handle)["project"]["version"]

    assert configured_version == gtlmapping.__version__
    assert configured_version == installed_version("GTLMapping")


def test_canonical_repository_url_is_published() -> None:
    project_root = Path(__file__).resolve().parents[1]
    with (project_root / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)["project"]

    assert project["urls"]["Repository"] == (
        "https://github.com/leongarcia2023/GTLMapping"
    )


def test_license_metadata_and_file_are_consistent() -> None:
    project_root = Path(__file__).resolve().parents[1]
    with (project_root / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)["project"]

    assert project["license"] == "BSD-3-Clause"
    assert project["license-files"] == ["LICENSE"]
    assert "BSD 3-Clause License" in (project_root / "LICENSE").read_text()
