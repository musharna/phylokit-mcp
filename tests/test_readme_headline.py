"""The README's headline table, re-run from the calls the README states.

Until 2026-11-10 the README said "seed 1" over a table that only simulation
seed 42 reproduces: seed 1 gives 56 and 7 informative sites and a correct
60-site tree, so the table's point (a wrong tree that looks as confident as a
right one) did not follow from the stated run. This test reads the tree, the
two calls and the table out of README.md and runs the calls, so the text and
the numbers cannot drift apart again.
"""

from __future__ import annotations

import re
from pathlib import Path

from phylokit_mcp.engine import IN_PROCESS_DRIFT
from phylokit_mcp.server import compare_trees, infer_tree, simulate_alignment

README = (Path(__file__).resolve().parent.parent / "README.md").read_text()


def _headline():
    section = README.split("## Why the rule", 1)[1].split("\n## ", 1)[0]
    newick = re.search(r"`(\(\(\(A:[^`]+;)`", section).group(1)
    sim = re.search(
        r'`simulate_alignment\(newick, model="(\w+)", length=sites, seed=(\d+)\)`',
        section,
    )
    inf = re.search(
        r'`infer_tree\(fasta, model="(\w+)", replicates=(\d+), seed=(\d+)\)`', section
    )
    rows = re.findall(
        r"^\|\s*(\d+) \|\s*(\d+) \| (yes, exactly|\*\*no, RF (\d+)\*\*)\s*\|\s*([\d.]+) \|$",
        section,
        re.MULTILINE,
    )
    return newick, sim, inf, rows, section


def test_the_readme_table_is_what_its_stated_calls_return():
    newick, sim, inf, rows, section = _headline()
    assert sim is not None and inf is not None, "the README no longer states its calls"
    assert len(rows) == 2, rows
    for sites, informative, _, rf, min_support in rows:
        aln = simulate_alignment(
            newick=newick, model=sim.group(1), length=int(sites), seed=int(sim.group(2))
        )
        assert aln["alignment"]["n_parsimony_informative"] == int(informative), sites
        tree = infer_tree(
            fasta=aln["fasta"],
            model=inf.group(1),
            replicates=int(inf.group(2)),
            seed=int(inf.group(3)),
        )
        cmp = compare_trees(tree["newick"], newick)
        assert cmp["robinson_foulds"] == (int(rf) if rf else 0), sites
        # Support drifts between calls by up to IN_PROCESS_DRIFT (engine.py).
        assert (
            abs(tree["support"]["min_support"] - float(min_support)) <= IN_PROCESS_DRIFT
        )
        if rf:
            # The prose under the table names the two clades.
            assert cmp["only_in_a"] == ["C|D|G"] and cmp["only_in_b"] == ["E|F|G"]
            assert "(`C,D,G`) that does not exist" in section
            assert "one that does (`E,F,G`)" in section
