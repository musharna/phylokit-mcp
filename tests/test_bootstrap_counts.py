"""`bootstrap_support`'s counting, asserted exactly, on scripted replicate trees.

The real engine's support values drift between calls (`IN_PROCESS_DRIFT`), so
the tests that drive it can only assert tolerances and orderings. That left the
arithmetic itself unasserted: the nightly mutation job (2026-11-08) found 31
mutants of `bootstrap_support` alive -- the count's starting value, both sort
orders, every field of a conflicting clade, the `>=` at 0.70 and 0.95, and the
empty-tree fallbacks.

Here the engine is replaced by a script of replicate trees, so every count is
known in advance: with 20 replicates, 14 is exactly 0.70 and 19 exactly 0.95.
Everything else is real: the column resampling, cogent3's alignment (which is
what refuses a wrong moltype) and the split canonicalisation. The last test
drives the real engine for what a script cannot stand in for: the defaults.
"""

from __future__ import annotations

import numpy as np
import pytest
from cogent3 import make_tree

from phylokit_mcp import bootstrap
from phylokit_mcp.alignment import AlignmentError
from phylokit_mcp.bootstrap import CladeSupport, bootstrap_support
from phylokit_mcp.inference import build_ml_tree

ML = "((A,B),(C,D),(E,F));"
SEQS = {t: "ACGTACGTAC" for t in "ABCDEF"}


class _ScriptedEngine:
    """Returns the scripted trees in order and records every call."""

    def __init__(self, newicks: list[str]):
        self.trees = [make_tree(n) for n in newicks]
        self.calls: list[tuple[str, dict]] = []

    def build_tree(self, alignment, model, **kwargs):
        self.calls.append((model, kwargs))
        return self.trees[len(self.calls) - 1]


def _run(monkeypatch, newicks, seqs=SEQS, ml=ML, **kwargs):
    engine = _ScriptedEngine(newicks)
    monkeypatch.setattr(bootstrap, "piqtree", lambda: engine)
    result = bootstrap_support(seqs, make_tree(ml), "JC", replicates=20, **kwargs)
    return result, engine


def test_ml_clades_are_counted_ranked_and_judged_at_the_lines_inclusively(
    monkeypatch,
):
    """C|D|E|F in 20 of 20, E|F in 19 (0.95) and C|D in 14 (0.70): both lines
    are met exactly, and both are inclusive. The 7 trees that disagree add
    D|E|F (0.25), C|E and D|F, none of which is a conflicting clade."""
    newicks = [ML] * 14 + ["((A,B),C,(D,(E,F)));"] * 5 + ["((A,B),(C,E),(D,F));"]
    result, engine = _run(monkeypatch, newicks, seed=5)

    assert result.clades == [
        CladeSupport(["C", "D", "E", "F"], 1.0, True),
        CladeSupport(["E", "F"], 0.95, True),
        CladeSupport(["C", "D"], 0.70, True),
    ]
    assert result.conflicting_clades == []
    assert result.n_ml_clades == 3
    assert result.n_supported == 3
    assert result.n_strongly_supported == 2
    assert result.mean_support == float(np.mean([1.0, 0.95, 0.70]))
    assert result.min_support == 0.70
    assert result.seed == 5
    # One search per replicate, every one on the caller's seed: the
    # randomness measured is the resampling, not the search's start.
    assert engine.calls == [("JC", {"rand_seed": 5})] * 20


def test_conflicting_clades_are_counted_from_zero_and_ranked_by_support(
    monkeypatch,
):
    """The ML tree has one clade, D|E|F. The replicates add C|D|E|F (16 of 20,
    0.80) and E|F (14, exactly 0.70), neither in the ML tree.

    Every split counted is nested in the next (E|F < D|E|F < C|D|E|F), so set
    inclusion orders them totally: a sort that falls back to comparing the
    clades themselves puts E|F first, every time. Only the descending sort by
    support puts C|D|E|F first. A first draft used two disjoint clades, and
    its order under that fallback moved with the hash seed (2026-11-08).
    """
    newicks = (
        ["(A,B,(C,(D,(E,F))));"] * 10
        + ["(A,B,(C,(D,E,F)));"] * 6
        + ["(A,B,C,(D,(E,F)));"] * 4
    )
    result, _ = _run(monkeypatch, newicks, ml="(A,B,C,(D,E,F));", seed=5)

    assert result.conflicting_clades == [
        CladeSupport(["C", "D", "E", "F"], 0.80, False),
        CladeSupport(["E", "F"], 0.70, False),
    ]
    assert result.clades == [CladeSupport(["D", "E", "F"], 1.0, True)]
    assert (result.n_supported, result.n_strongly_supported) == (1, 1)


def test_a_tree_with_no_clades_reports_none_not_a_number(monkeypatch):
    """Three taxa have no internal edge: nothing to support, so the mean is
    0.0 and the minimum is None, not nan, not 1.0, not an exception."""
    result, _ = _run(
        monkeypatch, ["(A,B,C);"] * 20, seqs={t: "ACGT" for t in "ABC"}, ml="(A,B,C);"
    )
    assert result.clades == []
    assert result.conflicting_clades == []
    assert result.n_ml_clades == 0
    assert result.mean_support == 0.0
    assert result.min_support is None
    assert result.fraction_resolved == 0.0


def test_the_defaults_are_seed_1_and_dna_on_the_real_engine():
    """No seed, no moltype: the real engine runs on DNA with seed 1.

    The negative half shows the moltype default is load-bearing: the same call
    with `"DNA"` is refused by the alignment, so a default other than `"dna"`
    would refuse every DNA caller who left it out.
    """
    seqs = {
        "a": "ACGTACGTACGTACGTACGTACGTACGTAC",
        "b": "ACGTACGTACGTACGTACGTACGTACGTAA",
        "c": "TCGAACGTTCGTACCTACGAACGTTCGTAC",
        "d": "TCGAACGTTCGTACCTACGAACGTTCGTAA",
    }
    tree = build_ml_tree(seqs, "JC", seed=1)
    result = bootstrap_support(seqs, tree, "JC", replicates=20)
    assert result.seed == 1
    assert len(result.clades) == 1
    with pytest.raises(AlignmentError) as exc:
        bootstrap_support(seqs, tree, "JC", replicates=20, moltype="DNA")
    assert str(exc.value) == "Unknown sequence_type 'DNA'. Valid: ['dna', 'protein']."
