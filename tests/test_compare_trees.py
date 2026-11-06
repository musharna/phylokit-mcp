"""compare_trees called directly, on trees small enough to work out by hand.

The other suites compare trees on the same taxa and read one or two fields.
These compare the whole result, including trees whose taxon sets differ, where
the distance is taken on the shared taxa after pruning the rest.
"""

from __future__ import annotations

import pytest

from phylokit_mcp.inputs import InputFileError
from phylokit_mcp.newick import TreeInputError
from phylokit_mcp.server import compare_trees

# On the shared taxa a-e, A is ((a,b),(c,d),e) and B is ((a,b),(c,e),d): both
# split {a,b} from {c,d,e}, and they disagree on {c,d} against {c,e}. x and y
# are each in one tree only; kept in, x would add the split {e,x} to A.
A = "((a,b),(c,d),(e,x));"
B = "((a,b),(c,e),(d,y));"

TAXA_DIFFER = (
    "The two trees do not have the same taxa, so the distance is computed on "
    "the shared subset only and is not comparable to an RF distance between "
    "trees on the full taxon set."
)


def test_trees_on_different_taxa_are_compared_on_the_shared_ones():
    assert compare_trees(A, B) == {
        "robinson_foulds": 2,
        "normalised_robinson_foulds": 0.5,
        "identical_topology": False,
        "n_shared_clades": 1,
        "only_in_a": ["c|d"],
        "only_in_b": ["c|e"],
        "shared_taxa": 5,
        "warnings": [
            {
                "code": "taxon_sets_differ",
                "message": TAXA_DIFFER,
                "only_in_a": ["x"],
                "only_in_b": ["y"],
            }
        ],
    }


def test_trees_on_the_same_taxa_get_no_warning():
    same = "((a,b),(c,d),e);"
    other = "((a,b),(c,e),d);"
    assert compare_trees(same, other) == {
        "robinson_foulds": 2,
        "normalised_robinson_foulds": 0.5,
        "identical_topology": False,
        "n_shared_clades": 1,
        "only_in_a": ["c|d"],
        "only_in_b": ["c|e"],
        "shared_taxa": 5,
        "warnings": [],
    }
    # A tree pruned to the shared taxa is the same tree as one written on them.
    assert compare_trees(A, same)["identical_topology"] is True


def test_four_shared_taxa_are_enough_and_three_are_not():
    four = compare_trees("((a,b),(c,d));", "((a,c),(b,d));")
    assert four["robinson_foulds"] == 2
    assert four["normalised_robinson_foulds"] == 1.0
    assert four["shared_taxa"] == 4
    with pytest.raises(ValueError) as exc:
        compare_trees("((a,b),(c,d));", "((a,b),(c,e));")
    assert str(exc.value) == (
        "Only 3 shared taxa; at least 4 are needed for an unrooted tree to have "
        "any internal edge to compare."
    )


def test_the_inputs_are_named_in_every_refusal(tmp_path):
    with pytest.raises(ValueError) as exc:
        compare_trees(A)
    assert str(exc.value) == "Give both newick_a and newick_b, or trees_path."
    with pytest.raises(ValueError) as exc:
        compare_trees(newick_b=B)
    assert str(exc.value) == "Give both newick_a and newick_b, or trees_path."

    path = tmp_path / "trees.nwk"
    path.write_text(f"{A}\n{B}\n")
    assert compare_trees(trees_path=str(path)) == compare_trees(A, B)
    with pytest.raises(ValueError) as exc:
        compare_trees(newick_b=B, trees_path=str(path))
    assert str(exc.value) == "Give newick_a and newick_b, or trees_path, not both."

    missing = tmp_path / "nope.nwk"
    with pytest.raises(InputFileError) as exc:
        compare_trees(trees_path=str(missing))
    assert str(exc.value) == f"trees_path: no such file: {missing}."

    with pytest.raises(TreeInputError) as exc:
        compare_trees(" ", B)
    assert str(exc.value) == "newick_a is empty."
    with pytest.raises(TreeInputError) as exc:
        compare_trees(A, "")
    assert str(exc.value) == "newick_b is empty."
