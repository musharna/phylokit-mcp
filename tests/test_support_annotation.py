"""`newick_with_support` and the `conflicting_support` advisory, against values
worked out by hand.

Until 2026-11-10 no test named either. Every `_annotate` mutant and every
`conflicting_support` mutant could survive, because the tests that ran them
checked the clade list and the warning codes, never the Newick labels or the
advisory's text.

The hand-derived cases come first. `_annotate` labels a node with the support
of the clade on its far side from the alphabetically first taxon, so on
TRUE_NEWICK the four internal nodes carry C|D|E|F|G, C|D, E|F|G and E|F, in
the order the string writes them. Then the real engine: on the hard fixture,
where the ML tree is wrong, every label must equal the clade list's value.
Last, the conflicting-clade path through `infer_tree`: the ML tree is the real
engine's, the replicate trees are a script, so each count is known in advance.
"""

from __future__ import annotations

from cogent3 import make_tree

from phylokit_mcp import bootstrap
from phylokit_mcp.bootstrap import BootstrapResult, CladeSupport
from phylokit_mcp.diagnostics import collect, conflicting_support
from phylokit_mcp.server import _annotate, infer_tree
from phylokit_mcp.splits import canonical_splits, format_split

from .conftest import TRUE_CLADES, TRUE_NEWICK
from .test_bootstrap_counts import _ScriptedEngine

TIPS = set("ABCDEFG")


def test_each_internal_node_is_labeled_with_its_own_clades_support():
    support = {"C|D": 0.57, "E|F": 0.8, "E|F|G": 0.914, "C|D|E|F|G": 1.0}
    assert set(support) == TRUE_CLADES
    out = _annotate(make_tree(TRUE_NEWICK), support, TIPS)
    # (A,B) is labeled with C|D|E|F|G: its far side from A. Two decimals,
    # rounded (0.914 -> 0.91), and branch lengths kept.
    assert out == (
        "(((A:0.05,B:0.05)1.00:0.05,(C:0.05,D:0.05)0.57:0.05)0.91:0.05,"
        "(E:0.05,F:0.05)0.80:0.05,G:0.1);"
    )


def test_labels_follow_the_tree_not_the_order_of_the_support_dict():
    """Same clades, rotated tree and reversed dict: each label moves with its
    clade. A label assigned by position would put 1.00 on the wrong node."""
    support = {"C|D|E|F|G": 0.6, "E|F|G": 0.7, "E|F": 0.85, "C|D": 0.95}
    rotated = (
        "(G:0.1,(F:0.05,E:0.05):0.05,((D:0.05,C:0.05):0.05,(B:0.05,A:0.05):0.05):0.05);"
    )
    assert _annotate(make_tree(rotated), support, TIPS) == (
        "(G:0.1,(F:0.05,E:0.05)0.85:0.05,((D:0.05,C:0.05)0.95:0.05,"
        "(B:0.05,A:0.05)0.60:0.05)0.70:0.05);"
    )


def test_the_labels_on_a_wrong_tree_are_the_clade_lists_values(fasta_hard):
    """Real engine, on the fixture where the ML tree is wrong (it contains
    C|D|G, which the true tree does not). Each label is read back from the
    returned string and matched to the clade it sits on."""
    result = infer_tree(fasta=fasta_hard, model="JC", replicates=20, seed=1)
    listed = {c["clade"]: f"{c['support']:.2f}" for c in result["support"]["clades"]}
    assert "C|D|G" in listed  # the wrong clade is there to be labeled

    labeled = {}
    tree = make_tree(result["newick_with_support"])
    tips = set(tree.get_tip_names())
    for edge in tree.get_edge_vector(include_root=False):
        if edge.is_tip():
            continue
        below = frozenset(edge.get_tip_names())
        side = below if "A" not in below else frozenset(tips) - below
        # cogent3 parses a numeric node label into `.support`, not the name.
        labeled[format_split(side)] = f"{edge.support:.2f}"
    assert labeled == listed
    # The annotated string is the same tree as the unannotated one.
    assert canonical_splits(tree, tips) == canonical_splits(
        make_tree(result["newick_topology"]), tips
    )


def _result(conflicting, supported_at=0.70):
    return BootstrapResult(
        replicates=20,
        clades=[CladeSupport(["C", "D"], 1.0, True)],
        conflicting_clades=conflicting,
        n_ml_clades=1,
        n_supported=1,
        n_strongly_supported=1,
        mean_support=1.0,
        min_support=1.0,
        seed=1,
        supported_at=supported_at,
    )


def test_the_conflicting_support_advisory_word_for_word():
    conflicting = [
        CladeSupport(["A", "C"], 0.75, False),
        CladeSupport(["B", "D"], 0.7, False),
    ]
    expected = {
        "code": "conflicting_support",
        "message": (
            "2 clade(s) reached at least 0.70 support across replicates but are "
            "ABSENT from the maximum-likelihood tree. The data support a grouping "
            "the reported topology does not show, so the single best tree is not "
            "the whole story."
        ),
        "clades": [
            {"taxa": ["A", "C"], "clade": "A|C", "support": 0.75, "in_ml_tree": False},
            {"taxa": ["B", "D"], "clade": "B|D", "support": 0.7, "in_ml_tree": False},
        ],
    }
    assert conflicting_support(_result(conflicting)) == expected
    assert expected in collect(result=_result(conflicting))
    # The line in the text is the method's own, not a constant.
    at95 = conflicting_support(_result(conflicting[:1], supported_at=0.95))
    assert at95 is not None
    assert at95["message"].startswith("1 clade(s) reached at least 0.95 support ")
    # Positive control: with nothing conflicting there is no advisory at all.
    assert conflicting_support(_result([])) is None
    codes = {w["code"] for w in collect(result=_result([]))}
    assert "conflicting_support" not in codes


def test_infer_tree_reports_a_conflicting_clade_it_was_scripted_to_find(
    monkeypatch, fasta_easy
):
    """The ML tree is the real engine's (the true tree, on the easy fixture).
    The 20 replicate trees are a script: 15 group A with C and B with D, 5 are
    the true tree. Anchored on A, the 15 contribute B|D and B|D|E|F|G, neither
    in the ML tree, each at 15/20 = 0.75; E|F and E|F|G are in all 20."""
    other = "(((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1):0.1,(E:0.1,F:0.1):0.1,G:0.1);"
    engine = _ScriptedEngine([other] * 15 + [TRUE_NEWICK] * 5)
    monkeypatch.setattr(bootstrap, "piqtree", lambda: engine)
    result = infer_tree(fasta=fasta_easy, model="JC", replicates=20, seed=1)
    assert len(engine.calls) == 20

    clades = {c["clade"]: c["support"] for c in result["support"]["clades"]}
    assert clades == {"C|D": 0.25, "C|D|E|F|G": 0.25, "E|F": 1.0, "E|F|G": 1.0}
    conflicting = [
        {"taxa": ["B", "D"], "clade": "B|D", "support": 0.75, "in_ml_tree": False},
        {
            "taxa": ["B", "D", "E", "F", "G"],
            "clade": "B|D|E|F|G",
            "support": 0.75,
            "in_ml_tree": False,
        },
    ]
    assert result["support"]["conflicting_clades"] == conflicting
    warned = [w for w in result["warnings"] if w["code"] == "conflicting_support"]
    assert len(warned) == 1
    assert warned[0]["clades"] == conflicting
    assert warned[0]["message"].startswith("2 clade(s) reached at least 0.70 ")
    # The Newick labels are the ML tree's own clades, never a conflicting one.
    assert result["newick_with_support"] == _annotate(
        make_tree(result["newick"]), clades, TIPS
    )
    assert "0.25" in result["newick_with_support"]
    assert "0.75" not in result["newick_with_support"]
