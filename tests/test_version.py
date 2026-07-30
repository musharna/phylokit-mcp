"""The version the package reports must equal the version it ships as.

There was no test here at all, which is the weaker half of the same problem
plantcv-mcp had: it *did* have one, but it asserted only that `__version__` was
a non-empty string — and a wrong version is a non-empty string. That suite went
green while the published 0.2.0 reported "0.1.0" to anyone who asked it.

phylokit's two declarations agreed, so nothing was wrong here yet. A latent
version of a defect that has already shipped elsewhere is still worth a guard.
"""

from __future__ import annotations

import tomllib
from importlib.metadata import version
from pathlib import Path

from phylokit_mcp import __version__

PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


def test_reported_version_matches_the_one_the_project_declares():
    """Compare against pyproject, which is the only independent source.

    `__version__` now reads from installed metadata, so asserting that those two
    agree would compare a value to itself and pass no matter what. pyproject is
    what a reintroduced literal would disagree with, which is the failure this
    is here to catch.
    """
    assert PYPROJECT.is_file(), f"expected pyproject.toml at {PYPROJECT}"
    declared = tomllib.loads(PYPROJECT.read_text())["project"]["version"]
    assert __version__ == declared, (
        f"package reports {__version__!r} but pyproject declares {declared!r}"
    )
    # If installed metadata disagrees with the checkout, the assertion above was
    # comparing against a stale install and proved nothing.
    assert version("phylokit-mcp") == declared
