"""`ufboot_tree` called directly, not through `infer_tree`.

`infer_tree` validates its arguments and always passes `seed` and `moltype`, so
through it this function's own defaults, its own replicate check and its own
seeding were never exercised: 18 of its mutants survived the 2026-11-02 nightly
(seed default 2, moltype default "XXdnaXX", the ufboot check run as
"bootstrap", `rand_seed=None`, clade order reversed, `in_ml_tree=False`, the
summary's seed None, the coverage check's message). Every test here runs real
IQ-TREE. Two wrap the real engine: one to see the seed it is handed, one to
corrupt the tree it returns.

The hard fixture is used throughout: a UFBoot run on it takes 1-2 s, against
7-9 s on the easy one (measured 2026-11-03), and nothing here needs a correct
topology.
"""

from __future__ import annotations

import pytest

from phylokit_mcp import bootstrap
from phylokit_mcp.alignment import AlignmentError, parse_fasta
from phylokit_mcp.bootstrap import UFBOOT, UFBOOT_MIN_REPLICATES, ufboot_tree
from phylokit_mcp.splits import format_split

from .test_protein import PROTEIN


@pytest.fixture(scope="module")
def hard(fasta_hard) -> dict[str, str]:
    return parse_fasta(fasta_hard)


def _supports(result) -> list[tuple[str, float]]:
    return [(format_split(frozenset(c.taxa)), c.support) for c in result.clades]


def test_the_defaults_are_seed_1_dna_and_the_floor_and_refuse_below_it(hard):
    """Called with only an alignment and a model, the run is the documented
    one: seed 1, DNA, 1000 replicates, every clade marked as in the ML tree,
    and the same clades as the call that spells those defaults out.

    The refusal is in the same test as the run it guards: a count below
    IQ-TREE's floor is refused as a ufboot count, and the floor runs.
    """
    _, result = ufboot_tree(hard, "JC")
    assert result.seed == 1
    assert result.replicates == UFBOOT_MIN_REPLICATES
    assert result.method == UFBOOT
    assert result.clades and all(c.in_ml_tree is True for c in result.clades)
    assert result.conflicting_clades == []
    _, explicit = ufboot_tree(
        hard, "JC", replicates=UFBOOT_MIN_REPLICATES, seed=1, moltype="dna"
    )
    # Same clades. Not the same values: IQ-TREE is not bit-exact at one seed
    # (engine.reproducibility("ufboot")), and on this fixture one clade moved 0.07
    # between two seed-1 runs (C|D|G 0.30, 0.23 on 2026-11-03).
    assert dict(_supports(result)).keys() == dict(_supports(explicit)).keys()
    # The "bootstrap" range (20-1000) would accept 100 and hand IQ-TREE a count
    # it refuses; the refusal must name the ufboot floor.
    with pytest.raises(ValueError, match="ufboot replicates must be between 1000"):
        ufboot_tree(hard, "JC", replicates=100)


def test_the_seed_reaches_iqtree_and_is_the_one_reported(hard, monkeypatch):
    """The engine is handed the caller's seed, and the result reports it.

    Watched at the call, not inferred from the output: IQ-TREE reads the wall
    clock, so two runs at one seed already differ (seed 1 gave E|F 0.80, then
    0.75, in one process on 2026-11-03), and that cannot be told apart from
    piqtree drawing a fresh random seed because it was handed None.
    """
    real_engine = bootstrap.piqtree
    handed: list[object] = []

    class Watched:
        def build_tree(self, *args, **kwargs):
            handed.append(kwargs.get("rand_seed"))
            return real_engine().build_tree(*args, **kwargs)

    monkeypatch.setattr(bootstrap, "piqtree", Watched)
    _, result = ufboot_tree(hard, "JC", seed=7)
    assert handed == [7]
    assert result.seed == 7
    assert result.clades


def test_clades_are_listed_strongest_first_then_by_name(hard):
    """The order a caller reads: highest support first, ties by clade name.

    Not vacuous: the hard case's clades sit far apart (0.98, 0.91, 0.80, 0.26
    on 2026-11-03), so an order by anything but support breaks it.
    """
    _, result = ufboot_tree(hard, "JC", seed=1)
    got = _supports(result)
    assert got == sorted(got, key=lambda cs: (-cs[1], cs[0]))
    assert len({v for _, v in got}) > 1, got


def test_a_protein_alignment_is_read_as_protein():
    """`moltype` reaches the alignment: protein residues under the DNA default
    are refused naming DNA, and the same rows declared protein give a tree."""
    with pytest.raises(AlignmentError, match="Unrecognised characters for dna"):
        ufboot_tree(PROTEIN, "LG")
    tree, result = ufboot_tree(PROTEIN, "LG", moltype="protein")
    assert set(tree.get_tip_names()) == set(PROTEIN)
    assert [format_split(frozenset(c.taxa)) for c in result.clades] == ["c|d"]


def test_engine_output_that_does_not_cover_the_tree_is_refused_by_name(
    hard, monkeypatch
):
    """Two ways IQ-TREE's annotation can stop matching its own tree, each made
    from a real run: a clade with no support label, and a support map missing a
    clade. Each is refused with its own message; the first names the replicate
    count the caller asked for. The unpatched run is the control.
    """
    real_engine, real_support = bootstrap.piqtree, bootstrap.support_by_split
    _, ok = ufboot_tree(hard, "JC")
    assert ok.clades

    class Unlabelled:
        """The real engine, with one internal edge's support label removed."""

        __version__ = "test"

        def build_tree(self, *args, **kwargs):
            tree = real_engine().build_tree(*args, **kwargs)
            labelled = [
                e
                for e in tree.get_edge_vector(include_root=False)
                if not e.is_tip() and getattr(e, "support", None) is not None
            ]
            labelled[0].support = None
            return tree

    monkeypatch.setattr(bootstrap, "piqtree", Unlabelled)
    with pytest.raises(
        RuntimeError, match="no UFBoot support .* after 1000 replicates"
    ):
        ufboot_tree(hard, "JC")
    monkeypatch.setattr(bootstrap, "piqtree", real_engine)

    def drop_one(tree, tips, replicates):
        support = real_support(tree, tips, replicates)
        support.pop(next(iter(support)))
        return support

    monkeypatch.setattr(bootstrap, "support_by_split", drop_one)
    with pytest.raises(
        RuntimeError,
        match=r"^UFBoot support did not cover every clade of the ML tree\.$",
    ):
        ufboot_tree(hard, "JC")
