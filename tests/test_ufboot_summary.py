"""The UFBoot2 summary and support reader, on inputs whose answer is known.

Through `ufboot_tree` these two only ever see real IQ-TREE support, which moves
between runs and never sits exactly on a threshold, and `ufboot_tree` passes
the same 0.95 as both lines. So 25 of their mutants survived the 2026-11-07
nightly: `>=` read as `>`, a count of 2 per clade, the strong line dropped,
the empty-tree defaults, a trivial split let in, the message's engine version.
Here `_summarise` gets support values chosen to sit on its lines, and
`support_by_split` gets Newick trees parsed by cogent3, the same parser that
reads IQ-TREE's output.

Four of `support_by_split`'s survivors cannot change its result, so no test
here is aimed at them. `include_root=None` and `include_root=False` are the
same to cogent3's `postorder` (`if include_self or node is not self`). With the
root included, the root's side of the split is empty, so the size check skips
it. A cogent3 `PhyloNode` always has `.support` (None when unlabelled), so
dropping `getattr`'s default changes nothing. `break` for `continue` stops the
scan at the first edge the size check skips; the only such edge is a child of
a two-way root covering all tips but one, and `get_edge_vector` is postorder,
so nothing but a tip comes after it.
"""

from __future__ import annotations

import pytest
from cogent3 import make_tree

from phylokit_mcp import bootstrap
from phylokit_mcp.bootstrap import UFBOOT, CladeSupport, _summarise, support_by_split
from phylokit_mcp.engine import piqtree
from phylokit_mcp.splits import format_split

TIPS = {"a", "b", "c", "d", "e"}


def test_a_clade_exactly_on_a_line_counts_and_both_lines_are_kept():
    clades = [
        CladeSupport(["a", "b"], 0.95, True),
        CladeSupport(["a", "b", "c"], 0.9, True),
        CladeSupport(["d", "e"], 0.8, True),
        CladeSupport(["c", "d"], 0.5, True),
    ]
    conflicting = [CladeSupport(["b", "c"], 0.75, False)]
    result = _summarise(clades, conflicting, 1000, 7, UFBOOT, 0.8, 0.9)

    assert result.n_ml_clades == 4
    # 0.8 is on the supported line and 0.9 on the strong one: both count.
    assert result.n_supported == 3
    assert result.n_strongly_supported == 2
    assert result.mean_support == pytest.approx(0.7875)
    assert result.min_support == 0.5  # the last clade, not the first
    assert (result.supported_at, result.strong_at) == (0.8, 0.9)
    assert (result.replicates, result.seed, result.method) == (1000, 7, UFBOOT)
    assert result.clades == clades
    assert result.conflicting_clades == conflicting


def test_a_tree_with_no_clades_summarises_to_zeros_and_no_minimum():
    result = _summarise([], [], 1000, 1, UFBOOT, 0.95, 0.95)

    assert (result.n_ml_clades, result.n_supported, result.n_strongly_supported) == (
        0,
        0,
        0,
    )
    assert result.mean_support == 0.0
    assert result.min_support is None


def test_each_clade_is_read_once_and_the_all_but_one_split_is_left_out():
    # Rooted on tip a: the edge above (b,c,d,e) splits off a alone, a trivial
    # split, so its label 80 is not a clade's support.
    tree = make_tree("(a,(((b,c)60,d)70,e)80);")

    assert support_by_split(tree, TIPS, 1000) == {
        frozenset({"b", "c"}): 0.6,
        frozenset({"b", "c", "d"}): 0.7,
    }


def _no_support_message(engine_version: str) -> str:
    return (
        f"IQ-TREE returned no UFBoot support for clade {format_split(frozenset({'c', 'd'}))} "
        f"after 1000 replicates (piqtree {engine_version}). "
        "The engine's output format has changed; the support cannot be reported."
    )


def test_an_unlabelled_clade_is_refused_naming_the_engine_version():
    labelled = make_tree("(a,b,((c,d)90,e)70);")
    assert support_by_split(labelled, TIPS, 1000) == {
        frozenset({"c", "d"}): 0.9,
        frozenset({"c", "d", "e"}): 0.7,
    }

    version = piqtree().__version__
    assert version[0].isdigit()  # the real engine reports one
    with pytest.raises(RuntimeError) as caught:
        support_by_split(make_tree("(a,b,((c,d),e)70);"), TIPS, 1000)
    assert str(caught.value) == _no_support_message(version)


def test_an_engine_without_a_version_is_named_with_a_question_mark(monkeypatch):
    monkeypatch.setattr(bootstrap, "piqtree", lambda: object())
    with pytest.raises(RuntimeError) as caught:
        support_by_split(make_tree("(a,b,((c,d),e)70);"), TIPS, 1000)
    assert str(caught.value) == _no_support_message("?")
