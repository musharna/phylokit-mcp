"""phylokit-mcp — phylogenetic inference over MCP, never a topology without support."""

from importlib.metadata import version

# Read from installed metadata rather than restated as a literal. The literal
# and pyproject.toml were two hand-maintained copies with nothing enforcing
# agreement — plantcv-mcp shipped 0.2.0 reporting "0.1.0" from exactly that
# setup. They agree here today; deriving removes the way they stop agreeing.
#
# No PackageNotFoundError fallback, deliberately: a sentinel like
# "0.0.0+unknown" answers a version question with a lie instead of an error.
__version__ = version("phylokit-mcp")

__all__ = ["__version__"]
