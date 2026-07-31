"""EVAL: given data generated on a KNOWN tree, is that tree recovered?

Every other test here checks that the code does what the code intends. This one
checks inference against a ground truth that existed before it ran.

**The simulator is written here, deliberately.** Using cogent3 to evolve the
sequences would share a library with the inference path; forty lines of
Jukes-Cantor share nothing. The generating process is also fully visible: each
branch of length `t` mutates a site with probability

    p = 3/4 · (1 − e^(−4t/3))

to a uniformly-chosen different base. That is JC69, the model IQ-TREE is then
asked to fit — so this measures whether inference recovers a topology from data
that genuinely came from one, not whether two implementations agree.

The score is Robinson-Foulds distance to the true topology. RF counts SPLITS, so
it is invariant to the many Newick spellings of one tree; string comparison
would answer a different question.
"""

import math
import random

import pytest
from cogent3 import make_tree

from phylokit_mcp.inference import build_ml_tree
from phylokit_mcp.splits import canonical_splits, robinson_foulds

BASES = "ACGT"
TIPS = ["a", "b", "c", "d"]

# Two regimes, both MEASURED before being written down (10 seeds each):
#
#   internal  tip   1000 sites  100 sites  30 sites
#   0.25      0.05  10/10       10/10      10/10
#   0.02      0.30   9/10        5/10       3/10
#
# The discriminating variable is the INTERNAL EDGE, not the alignment length —
# which is this server's own doctrine ("alignment length is not evidence")
# arrived at from the other direction.
EASY = (0.25, 0.05)
HARD = (0.02, 0.30)


def _evolve(seq: str, branch_length: float, rng: random.Random) -> str:
    """JC69 along one branch."""
    p = 0.75 * (1.0 - math.exp(-4.0 * branch_length / 3.0))
    return "".join(
        rng.choice([b for b in BASES if b != ch]) if rng.random() < p else ch
        for ch in seq
    )


def simulate(n_sites: int, internal: float, tip: float, seed: int) -> dict[str, str]:
    """Evolve `n_sites` down ((a,b),(c,d)) and return the tip alignment."""
    rng = random.Random(seed)
    root = "".join(rng.choice(BASES) for _ in range(n_sites))
    left = _evolve(root, internal, rng)  # ancestor of (a, b)
    right = _evolve(root, internal, rng)  # ancestor of (c, d)
    return {
        "a": _evolve(left, tip, rng),
        "b": _evolve(left, tip, rng),
        "c": _evolve(right, tip, rng),
        "d": _evolve(right, tip, rng),
    }


def _rf_to_truth(seqs: dict[str, str], internal: float, tip: float, seed: int) -> int:
    tree = build_ml_tree(seqs, "JC", seed=seed)
    truth = make_tree(f"((a:{tip},b:{tip}):{internal},(c:{tip},d:{tip}):{internal});")
    return robinson_foulds(canonical_splits(tree, TIPS), canonical_splits(truth, TIPS))


def _recovery_rate(n_sites: int, regime: tuple[float, float], seeds) -> int:
    internal, tip = regime
    return sum(
        _rf_to_truth(simulate(n_sites, internal, tip, s), internal, tip, s) == 0
        for s in seeds
    )


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_the_true_topology_is_recovered_from_data_generated_on_it(seed):
    """RF distance to truth must be zero on well-conditioned data."""
    internal, tip = EASY
    got = _rf_to_truth(simulate(1000, internal, tip, seed), internal, tip, seed)
    assert got == 0, (
        "inference did not recover the topology the data was generated on, at "
        "1000 sites with a 0.25 internal edge — conditions under which the "
        f"answer is not in doubt (RF={got})"
    )


def test_the_eval_can_fail_when_the_signal_is_destroyed():
    """THE control, and the reason the test above is evidence at all.

    Shuffling each SITE independently across taxa preserves every base
    composition and every alignment dimension while destroying shared ancestry.
    If recovery still succeeded, the test above would be measuring something
    other than phylogenetic signal — a default topology, or tip order.

    Inference still returns a tree, confidently. That is the point: an
    unsupported topology looks exactly like a supported one from the outside.
    """
    internal, tip = EASY
    rng = random.Random(99)
    seqs = simulate(1000, internal, tip, 1)
    columns = [[seqs[t][i] for t in TIPS] for i in range(len(seqs["a"]))]
    for col in columns:
        rng.shuffle(col)
    shuffled = {t: "".join(col[k] for col in columns) for k, t in enumerate(TIPS)}

    # Sanity: the shuffle preserved shape and base composition exactly.
    assert sorted("".join(shuffled.values())) == sorted("".join(seqs.values()))

    distances = {_rf_to_truth(shuffled, internal, tip, s) for s in (1, 2, 3)}
    assert distances != {0}, (
        "the true topology was recovered from an alignment whose phylogenetic "
        "signal was destroyed — the recovery test is not measuring ancestry"
    )


def test_the_eval_has_resolution_on_the_internal_edge():
    """A benchmark that always scores full marks cannot detect a regression.

    Measured: the easy regime recovers 10/10 at 1000 sites AND at 30, so
    alignment length is not what carries the topology here. Shortening the
    INTERNAL EDGE while lengthening the tips is what makes the problem hard —
    exactly the case this server warns about, reached from the other side.
    """
    seeds = range(1, 11)
    easy_long = _recovery_rate(1000, EASY, seeds)
    easy_short_alignment = _recovery_rate(30, EASY, seeds)
    hard = _recovery_rate(100, HARD, seeds)

    assert easy_long == 10, f"well-conditioned data should always recover: {easy_long}"
    assert easy_short_alignment == 10, (
        f"a 30-site alignment on a 0.25 internal edge recovered "
        f"{easy_short_alignment}/10 — LENGTH is not the signal, and this eval "
        "should show that rather than conflating the two"
    )
    assert hard < easy_long, (
        f"a 0.02 internal edge with 0.30 tips recovered {hard}/10, the same as "
        "well-conditioned data. This eval has no resolution and would not notice "
        "inference getting worse."
    )
