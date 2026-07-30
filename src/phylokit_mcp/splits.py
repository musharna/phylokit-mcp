"""Bipartitions — the unit in which topology agreement is actually measured.

Two trees are compared by the SPLITS they contain, not by their Newick strings.
The same topology has many valid Newick representations (rotate any child order,
re-root anywhere on an unrooted tree) so string comparison answers a different
question from the one a caller means when they ask whether two trees agree.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Iterable

Split = frozenset[str]


def canonical_splits(tree, tip_names: Iterable[str]) -> set[Split]:
    """The non-trivial unrooted bipartitions of ``tree``.

    Each internal edge divides the taxa in two. Which half you name is arbitrary
    on an unrooted tree, so a split is canonicalised by anchoring on a fixed
    reference taxon — the half WITHOUT it is the representative. Without that
    anchoring the same biological split gets two different identities depending
    on where the tree happens to be rooted, and support counts silently split
    between them.

    Splits with fewer than two taxa on a side are dropped: every tree on the same
    taxa contains those, so they carry no information about topology and would
    inflate any agreement measure toward 1.
    """
    tips = frozenset(tip_names)
    if len(tips) < 4:
        # An unrooted tree on <4 taxa has no internal edge, so no split exists
        # to support. This is not an error; it is an empty set.
        return set()
    ref = min(tips)
    out: set[Split] = set()
    for edge in tree.get_edge_vector(include_root=False):
        if edge.is_tip():
            continue
        below = frozenset(t.name for t in edge.tips())
        side = below if ref not in below else tips - below
        if 2 <= len(side) <= len(tips) - 2:
            out.add(side)
    return out


def format_split(split: Split) -> str:
    """A stable, sorted, human-readable name for a clade."""
    return "|".join(sorted(split))


def robinson_foulds(splits_a: set[Split], splits_b: set[Split]) -> int:
    """Symmetric difference of split sets — the unnormalised RF distance."""
    return len(splits_a ^ splits_b)


def normalised_robinson_foulds(
    splits_a: set[Split], splits_b: set[Split]
) -> float | None:
    """RF as a fraction of the maximum possible for these trees.

    Returns None when the maximum is zero — i.e. neither tree has any internal
    structure to disagree about. Reporting 0.0 there would read as "identical",
    which is a stronger claim than "there was nothing to compare".
    """
    denom = len(splits_a) + len(splits_b)
    if denom == 0:
        return None
    return robinson_foulds(splits_a, splits_b) / denom
