"""Reading a caller-supplied Newick string -- the tree tools' input boundary.

cogent3 is the parser, and it has two behaviours a caller cannot see through the
server. A malformed string raises its TreeParseError, which is not a ValueError,
so the server reported it as `Error executing tool` with the reason dropped. And
a repeated tip label is silently renamed (`a`, `a.2`), so two trees that name the
same taxon twice were compared as if they held a taxon called `a.2` -- a
distance over a taxon that does not exist, returned as a normal result.

Both are the caller's input, so both are refused here, as TreeInputError, with
the parser's own message.
"""

from __future__ import annotations

from collections import Counter
from typing import Any


class TreeInputError(ValueError):
    """A Newick string the tree tools cannot use, with the reason named."""


def _raw_tip_labels(text: str) -> list[str | None]:
    """Tip labels exactly as written, before cogent3 renames duplicates.

    Uses cogent3's own tokeniser and tree-builder callback, so what counts as a
    tip and a label is the parser's definition rather than a second one written
    here; the recording subclass only watches what it is handed.
    """
    from cogent3.core.tree import TreeBuilder
    from cogent3.parse.newick import parse_string

    labels: list[str | None] = []

    class _Recorder(TreeBuilder):
        def create_edge(self, children, name, *args: Any, **kwargs: Any):
            if not children:
                labels.append(name)
            return super().create_edge(children, name, *args, **kwargs)

    parse_string(text, _Recorder().create_edge)
    return labels


def parse_newick(text: str, label: str = "newick"):
    """Parse `text` into a cogent3 tree, or refuse it with the reason."""
    from cogent3 import make_tree
    from cogent3.parse.newick import TreeParseError

    text = text.strip()
    if not text:
        raise TreeInputError(f"{label} is empty.")
    try:
        tree = make_tree(text)
        tips = _raw_tip_labels(text)
    except TreeParseError as exc:
        raise TreeInputError(f"{label} is not valid Newick: {exc}") from exc

    unnamed = sum(1 for t in tips if not t)
    if unnamed:
        raise TreeInputError(
            f"{label} has {unnamed} unnamed tip(s). Every tip needs a label, "
            "or the trees cannot be compared taxon by taxon."
        )
    dupes = sorted(n for n, c in Counter(tips).items() if c > 1)  # type: ignore[type-var]
    if dupes:
        raise TreeInputError(
            f"{label} names the same tip more than once: {dupes[:8]}. A tree "
            "cannot hold one taxon at two places; the parser would otherwise "
            "rename the repeat (e.g. 'a' -> 'a.2') and the result would describe "
            "a taxon that is not in your data."
        )
    return tree
