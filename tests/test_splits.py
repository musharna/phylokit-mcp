"""Splits are the foundation: if canonicalisation is wrong, every support value
and every RF distance downstream is wrong in a way that still looks plausible.
"""

from __future__ import annotations

import pytest
from cogent3 import make_tree

from phylokit_mcp.splits import (
    canonical_splits,
    format_split,
    normalised_robinson_foulds,
    robinson_foulds,
)

from .conftest import TRUE_CLADES, TRUE_NEWICK


def _named(newick: str) -> set[str]:
    tree = make_tree(newick)
    tips = set(tree.get_tip_names())
    return {format_split(s) for s in canonical_splits(tree, tips)}


def test_extracts_exactly_the_expected_clades():
    assert _named(TRUE_NEWICK) == TRUE_CLADES


def test_an_unrooted_tree_on_n_taxa_has_n_minus_three_splits():
    # The count is a hard combinatorial fact for a fully resolved unrooted tree.
    # Getting more means trivial splits leaked in; fewer means real ones were
    # dropped. Either way support counts would be computed against the wrong set.
    assert len(_named(TRUE_NEWICK)) == 7 - 3


def test_rerooting_does_not_change_the_splits():
    """The same topology written from a different root must have one identity.

    This is the bug canonicalisation exists to prevent: without anchoring, a
    bootstrap replicate rooted elsewhere contributes its counts to a DIFFERENT
    key than the ML tree's clade, and support silently reads as zero.
    """
    rerooted = make_tree(TRUE_NEWICK).rooted_with_tip("G")
    tips = set(rerooted.get_tip_names())
    assert {format_split(s) for s in canonical_splits(rerooted, tips)} == TRUE_CLADES


def test_child_order_does_not_change_the_splits():
    swapped = (
        "(((B:0.05,A:0.05):0.05,(D:0.05,C:0.05):0.05):0.05,G:0.1,(F:0.05,E:0.05):0.05);"
    )
    assert _named(swapped) == TRUE_CLADES


def test_a_tree_rooted_at_a_tip_drops_the_trivial_split():
    """One-taxon-versus-the-rest is not information, and must not be counted.

    Found by a surviving mutant: removing the trivial-split guard changed
    nothing on the balanced 7-taxon fixture, because that shape never produces
    such a split. A LADDER ROOTED AT A TIP does — here the edge below A
    subtends {B,C,D,E}, i.e. every other taxon. Every tree on these taxa
    contains that split, so counting it would add a guaranteed-1.00 clade to
    every result and inflate `fraction_resolved` for free.
    """
    tree = make_tree("(A:0.1,(B:0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1):0.1);")
    tips = set(tree.get_tip_names())
    named = {format_split(s) for s in canonical_splits(tree, tips)}
    assert "B|C|D|E" not in named
    assert named == {"D|E", "C|D|E"}
    # 5 taxa unrooted -> exactly 2 internal splits, regardless of rooting.
    assert len(named) == 5 - 3


def test_the_same_topology_gives_one_split_set_across_three_rootings():
    """Rooted-balanced, rooted-at-a-tip and unrooted must agree exactly.

    A rooted tree yields the SAME bipartition from both root children; only
    returning a set keeps that from being counted twice, which would show up as
    a duplicated clade in every support table.
    """
    shapes = [
        "(A:0.1,(B:0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1):0.1)",
        "((A:0.1,B:0.1):0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1)",
        "(A:0.1,B:0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1)",
    ]
    seen = set()
    for shape in shapes:
        tree = make_tree(shape + ";")
        tips = set(tree.get_tip_names())
        splits = canonical_splits(tree, tips)
        assert len(splits) == 2, f"{shape} gave {len(splits)} splits, expected 2"
        seen.add(frozenset(format_split(s) for s in splits))
    assert len(seen) == 1, f"rootings disagreed: {seen}"


@pytest.mark.parametrize("n_taxa", [2, 3])
def test_too_few_taxa_yields_no_splits_rather_than_an_error(n_taxa: int):
    newick = "(" + ",".join(f"T{i}:0.1" for i in range(n_taxa)) + ");"
    tree = make_tree(newick)
    assert canonical_splits(tree, set(tree.get_tip_names())) == set()


def test_robinson_foulds_is_zero_for_the_same_topology_written_differently():
    a = canonical_splits(make_tree(TRUE_NEWICK), set("ABCDEFG"))
    swapped = (
        "(((B:0.05,A:0.05):0.05,(D:0.05,C:0.05):0.05):0.05,G:0.1,(F:0.05,E:0.05):0.05);"
    )
    b = canonical_splits(make_tree(swapped), set("ABCDEFG"))
    assert robinson_foulds(a, b) == 0
    assert normalised_robinson_foulds(a, b) == 0.0


def test_robinson_foulds_counts_both_directions_of_disagreement():
    a = canonical_splits(make_tree(TRUE_NEWICK), set("ABCDEFG"))
    # Move G next to the AB clade instead of the EF clade.
    other = (
        "(((A:0.05,B:0.05):0.05,G:0.1):0.05,(C:0.05,D:0.05):0.05,(E:0.05,F:0.05):0.05);"
    )
    b = canonical_splits(make_tree(other), set("ABCDEFG"))
    rf = robinson_foulds(a, b)
    assert rf > 0
    # Symmetric difference, so it must equal the two one-way differences summed.
    assert rf == len(a - b) + len(b - a)


def test_normalised_rf_returns_none_when_there_is_nothing_to_compare():
    """Zero would read as 'identical', which is a stronger claim than the truth."""
    assert normalised_robinson_foulds(set(), set()) is None
