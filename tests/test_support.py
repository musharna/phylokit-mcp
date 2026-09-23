"""The tests that carry this server's argument.

The claim is that a maximum-likelihood topology, on its own, gives a caller no
way to tell a correct tree from an incorrect one — and that bootstrap support
does. Testing that requires a case where the tree is ACTUALLY WRONG. A suite
built only on easy alignments would pass against a `support()` that returned 1.0
unconditionally, which is the mutant these tests exist to kill.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

from phylokit_mcp.bootstrap import (
    MAX_REPLICATES,
    MIN_REPLICATES,
    WEAK_SUPPORT,
    resample_columns,
    validate_replicates,
)
from phylokit_mcp.engine import IN_PROCESS_DRIFT
from phylokit_mcp.server import compare_trees, infer_tree

from .conftest import TRUE_CLADES


def test_a_long_alignment_recovers_the_true_tree(fasta_easy, easy_alignment):
    result = infer_tree(fasta=fasta_easy, model="JC", replicates=50, seed=1)
    comparison = compare_trees(
        newick_a=result["newick"], newick_b=easy_alignment["true_newick"]
    )
    assert comparison["identical_topology"], "the easy case must recover the truth"
    assert comparison["robinson_foulds"] == 0
    assert result["support"]["fraction_resolved"] == 1.0


def test_a_short_alignment_produces_a_WRONG_tree(fasta_hard, hard_alignment):
    """The positive control for the whole premise.

    If this ever starts passing the RF==0 check, the 'hard' fixture has stopped
    being hard and every test below it is only measuring the easy regime.
    """
    result = infer_tree(fasta=fasta_hard, model="JC", replicates=100, seed=1)
    comparison = compare_trees(
        newick_a=result["newick"], newick_b=hard_alignment["true_newick"]
    )
    assert not comparison["identical_topology"]
    assert comparison["robinson_foulds"] > 0


def test_support_separates_the_false_clade_from_the_true_ones(fasta_hard):
    """The product claim, stated as an assertion.

    On the wrong tree, the clade that does not exist in the true topology must
    carry LOWER support than every clade that does. This is what makes the
    support values load-bearing rather than decorative.
    """
    result = infer_tree(fasta=fasta_hard, model="JC", replicates=100, seed=1)
    clades = {c["clade"]: c["support"] for c in result["support"]["clades"]}
    false_clades = {k: v for k, v in clades.items() if k not in TRUE_CLADES}
    true_clades = {k: v for k, v in clades.items() if k in TRUE_CLADES}

    assert false_clades, "fixture must yield at least one clade absent from the truth"
    assert true_clades, "fixture must retain at least one true clade"
    assert max(false_clades.values()) < min(true_clades.values()), (
        f"a false clade was supported at least as well as a true one: "
        f"false={false_clades} true={true_clades}"
    )


def test_the_two_cases_are_indistinguishable_without_support(fasta_easy, fasta_hard):
    """Both runs return a fully resolved Newick string of the same shape.

    This is the negative half of the argument: the tree ALONE cannot tell the
    caller which run to trust, so returning one without support would be
    returning a result that cannot be evaluated.
    """
    easy = infer_tree(fasta=fasta_easy, model="JC", replicates=MIN_REPLICATES, seed=1)
    hard = infer_tree(fasta=fasta_hard, model="JC", replicates=MIN_REPLICATES, seed=1)
    assert easy["newick"].count(",") == hard["newick"].count(",")
    assert len(easy["support"]["clades"]) == len(hard["support"]["clades"])
    # ...and yet the support tells them apart.
    assert hard["support"]["min_support"] < easy["support"]["min_support"]


def test_unsupported_clades_warning_names_the_weak_clade(fasta_hard):
    result = infer_tree(fasta=fasta_hard, model="JC", replicates=100, seed=1)
    warned = [w for w in result["warnings"] if w["code"] == "unsupported_clades"]
    assert warned, "a tree with a sub-threshold clade must warn"
    flagged = {c["clade"] for c in warned[0]["clades"]}
    weak = {
        c["clade"] for c in result["support"]["clades"] if c["support"] < WEAK_SUPPORT
    }
    assert flagged == weak


def test_a_well_resolved_tree_raises_no_support_warnings(fasta_easy):
    """The negative control: the advisories must be capable of staying silent.

    Without this, a `collect()` that returned every warning unconditionally
    would pass every other test in this file.
    """
    result = infer_tree(fasta=fasta_easy, model="JC", replicates=50, seed=1)
    codes = {w["code"] for w in result["warnings"]}
    assert not codes & {"unsupported_clades", "poorly_resolved", "conflicting_support"}


def test_repeated_calls_agree_to_within_the_measured_engine_drift(fasta_hard):
    """Same seed, same input, repeated IN ONE PROCESS.

    This asserts a tolerance rather than equality because bit-exact equality is
    a property IQ-TREE does not have: it reads the wall clock during its search
    (see the frozen-clock control below), so building the same tree three times
    in a row can give call 1 == call 2 but call 3 different. Measured over six repeated
    50-replicate calls, three of four clades were bit-identical and one moved
    0.02 — one replicate flipping, well inside the bootstrap's own sampling
    error. The engine module documents the same number, and the server reports
    it to callers rather than claiming a determinism it does not have.
    """
    a = infer_tree(fasta=fasta_hard, model="JC", replicates=50, seed=7)
    b = infer_tree(fasta=fasta_hard, model="JC", replicates=50, seed=7)
    assert [c["clade"] for c in a["support"]["clades"]] == [
        c["clade"] for c in b["support"]["clades"]
    ], "the TOPOLOGY must be identical; only support values may drift"
    for x, y in zip(a["support"]["clades"], b["support"]["clades"], strict=True):
        assert abs(x["support"] - y["support"]) <= IN_PROCESS_DRIFT


# One fresh interpreter, one MCP session, one infer_tree call. Reports what a
# caller would compare: branch lengths (in the newick), the log-likelihood, the
# support values, and what the server claims about reproducibility.
_FRESH_RUN = """
import json, sys
import anyio
from mcp.client.client import Client
from phylokit_mcp.server import build_server

async def go():
    async with Client(build_server()) as c:
        r = await c.call_tool("infer_tree", {
            "fasta": sys.stdin.read(), "model": "JC", "replicates": 20, "seed": 7,
        })
        assert not r.is_error, r.content[0].text
        return r.structured_content

out = anyio.run(go)
print(json.dumps({
    "newick": out["newick"],
    "lnL": out["log_likelihood"],
    "support": [[c["clade"], c["support"]] for c in out["support"]["clades"]],
    "claim": out["reproducibility"],
}))
"""

# gettimeofday() pinned to a constant. Freezing it -- and only it; time(),
# clock() and getrusage() were each tried alone and changed nothing -- made
# IQ-TREE bit-exact across and within processes, which is how the source of the
# drift was identified.
_FROZEN_CLOCK_C = """
#include <sys/time.h>
int gettimeofday(struct timeval *tv, void *tz) {
    (void)tz; tv->tv_sec = 1000000000; tv->tv_usec = 0; return 0;
}
"""


def _fresh_runs(n: int, fasta: str, env: dict | None = None) -> list[dict]:
    """`n` fresh server processes, at most 3 at a time."""
    import json
    import os
    from concurrent.futures import ThreadPoolExecutor

    def one(_):
        proc = subprocess.run(
            [sys.executable, "-c", _FRESH_RUN],
            input=fasta,
            capture_output=True,
            text=True,
            timeout=600,
            env={**os.environ, **(env or {})},
            check=False,  # the returncode is asserted below, with stderr
        )
        assert proc.returncode == 0, proc.stderr[-2000:]
        return json.loads(proc.stdout.strip().splitlines()[-1])

    with ThreadPoolExecutor(max_workers=3) as pool:
        return list(pool.map(one, range(n)))


def _numbers(run: dict) -> tuple:
    return (run["newick"], run["lnL"], tuple(map(tuple, run["support"])))


def test_with_the_clock_frozen_fresh_processes_are_bit_identical(fasta_hard, tmp_path):
    """The positive control for the test below, and the proof of the cause.

    If the comparison could not see equality, "the runs differ" would be a
    property of the harness. Here the same harness, with only the wall clock
    taken away from IQ-TREE, gets byte-identical branch lengths, likelihood and
    support from every fresh process.
    """
    import shutil

    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        pytest.skip("no C compiler to build the frozen-clock shim")
    src = tmp_path / "frozen_clock.c"
    lib = tmp_path / "frozen_clock.so"
    src.write_text(_FROZEN_CLOCK_C)
    subprocess.run([cc, "-shared", "-fPIC", "-o", str(lib), str(src)], check=True)

    runs = _fresh_runs(4, fasta_hard, env={"LD_PRELOAD": str(lib)})
    assert len({_numbers(r) for r in runs}) == 1, [r["lnL"] for r in runs]


def test_fresh_processes_are_not_bit_exact_and_the_server_says_so(fasta_hard):
    """M1 of the 2026-09-22 audit. The claim must match the measurement.

    The server reported `deterministic_across_processes: True`, backed by a test
    that compared SUPPORT VALUES only -- the one output coarse enough to hide
    the drift. Branch lengths and the log-likelihood differ between fresh
    processes given the same seed (IQ-TREE reads the wall clock during its
    search; see the frozen-clock control above). This compares those, not just
    support.

    Measured on this fixture: 8 fresh processes gave 3 distinct results, so 6
    runs all agreeing by chance is about 2%; a second batch of 6 takes it to
    about 1 in 2000. If 12 runs ever agree, the engine has plausibly become
    deterministic -- re-measure and revisit the claim, do not loosen this.
    """
    runs = _fresh_runs(6, fasta_hard)
    if len({_numbers(r) for r in runs}) == 1:
        runs += _fresh_runs(6, fasta_hard)
    distinct = {_numbers(r) for r in runs}
    assert len(distinct) > 1, (
        "12 fresh processes returned bit-identical trees. The engine may now be "
        "deterministic across processes; re-measure before changing the claim."
    )

    for r in runs:
        claim = r["claim"]
        assert claim["deterministic_across_processes"] is False, (
            "the runs above differ, so the server must not claim cross-process "
            f"determinism; it said {claim}"
        )
        assert claim["bit_exact_on_repeat_within_process"] is False

    # What IS stable, stated as the tolerance the server reports: the topology
    # and the clade list, with support within the advertised drift.
    first = runs[0]
    for r in runs[1:]:
        assert [c for c, _ in r["support"]] == [c for c, _ in first["support"]]
        for (_, x), (_, y) in zip(r["support"], first["support"], strict=True):
            assert abs(x - y) <= first["claim"]["max_support_drift_across_processes"]


def test_a_different_seed_changes_support_but_not_wildly(fasta_hard):
    """Support is a sampling estimate, so it must MOVE with the seed.

    Identical values across seeds would mean the resampling is not actually
    varying — the failure mode where a bug makes every replicate the original
    alignment, which yields support 1.0 everywhere and looks like success.
    """
    a = infer_tree(fasta=fasta_hard, model="JC", replicates=50, seed=7)
    b = infer_tree(fasta=fasta_hard, model="JC", replicates=50, seed=99)
    sa = [c["support"] for c in a["support"]["clades"]]
    sb = [c["support"] for c in b["support"]["clades"]]
    assert sa != sb, "support did not move with the seed; is resampling happening?"
    assert max(abs(x - y) for x, y in zip(sa, sb, strict=False)) < 0.5


def test_resampling_draws_columns_with_replacement():
    """Columns must be kept intact and drawn WITH replacement.

    Resampling characters independently within a column would destroy the
    cross-taxon correlation that carries phylogenetic signal, while still
    producing a plausible-looking alignment.
    """
    seqs = {"a": "AAAA", "b": "CCCC", "c": "GGGG", "d": "TTTT"}
    out = resample_columns(seqs, np.random.default_rng(0))
    assert {len(v) for v in out.values()} == {4}
    # Each taxon is constant in this construction, so any column permutation
    # preserves it exactly. A character-level shuffle across taxa would not.
    assert out == seqs

    varied = {"a": "ACGT", "b": "ACGT", "c": "ACGT", "d": "ACGT"}
    draws = {
        tuple(resample_columns(varied, np.random.default_rng(s))["a"])
        for s in range(20)
    }
    assert len(draws) > 1, "with replacement, 20 seeds must not give one alignment"
    assert any(len(set(d)) < 4 for d in draws), (
        "with replacement, some draw must repeat a column"
    )


@pytest.mark.parametrize("bad", [0, 1, MIN_REPLICATES - 1, MAX_REPLICATES + 1])
def test_replicate_count_is_bounded(bad: int):
    with pytest.raises(ValueError, match="replicates must be between"):
        validate_replicates(bad)


def test_replicate_bounds_are_inclusive():
    validate_replicates(MIN_REPLICATES)
    validate_replicates(MAX_REPLICATES)
