"""Nonparametric bootstrap support (Felsenstein 1985, Evolution 39:783-791).

This is computed here rather than read back from IQ-TREE because piqtree 0.8.3
runs `bootstrap_replicates` but does not expose the resulting support values:
internal nodes come back carrying only `mprobs`, node names are `edge.N`, and
nothing is written to disk. Measured on 0.8.3 — if a later release exposes them,
this module becomes a fallback rather than the only route.

Computing it here has a side benefit worth keeping either way: the full split
frequency distribution is available, including well-supported clades that are
NOT in the maximum-likelihood tree. A conflicting clade at 0.70 is the single
most informative thing about a tree and no support-annotated Newick string can
express it, because such a clade has nowhere to attach.

Why this matters, measured on a simulated 7-taxon tree with a known topology:

    sites   ML topology     lowest true clade   false clade in ML tree
    300     correct         1.00                --
     60     WRONG           0.86                CDG at 0.57

At 60 sites the maximum-likelihood tree contains a clade that does not exist and
omits one that does. The two runs produce topologies that look equally
authoritative. The support values are the only thing distinguishing them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .engine import piqtree
from .splits import Split, canonical_splits, format_split

# Below this a support value is mostly sampling noise in its own right: the
# standard error of a proportion at n=100 is about 0.05, so ranking clades that
# differ by less than that across 20 replicates is reading tea leaves.
MIN_REPLICATES = 20
MAX_REPLICATES = 1000
DEFAULT_REPLICATES = 100

# The conventional reading of a bootstrap proportion. 0.70 is the customary
# "supported" line and 0.95 the strong one; both are conventions, not theorems,
# and the server says so rather than implying a test was performed.
WEAK_SUPPORT = 0.70
STRONG_SUPPORT = 0.95


@dataclass(frozen=True)
class CladeSupport:
    taxa: list[str]
    support: float
    in_ml_tree: bool

    def as_dict(self) -> dict:
        return {
            "taxa": self.taxa,
            "clade": format_split(frozenset(self.taxa)),
            "support": round(self.support, 3),
            "in_ml_tree": self.in_ml_tree,
        }


@dataclass(frozen=True)
class BootstrapResult:
    replicates: int
    clades: list[CladeSupport]
    conflicting_clades: list[CladeSupport]
    n_ml_clades: int
    n_supported: int
    n_strongly_supported: int
    mean_support: float
    min_support: float | None
    seed: int

    @property
    def fraction_resolved(self) -> float:
        """Share of the ML tree's clades that clear the conventional 0.70 line.

        This is the headline number: a tree whose clades are mostly unsupported
        is a drawing, not a result.
        """
        if self.n_ml_clades == 0:
            return 0.0
        return self.n_supported / self.n_ml_clades

    def as_dict(self) -> dict:
        return {
            "replicates": self.replicates,
            "method": "nonparametric bootstrap over alignment columns",
            "seed": self.seed,
            "n_clades": self.n_ml_clades,
            "n_supported": self.n_supported,
            "n_strongly_supported": self.n_strongly_supported,
            "fraction_resolved": round(self.fraction_resolved, 3),
            "mean_support": round(self.mean_support, 3),
            "min_support": None
            if self.min_support is None
            else round(self.min_support, 3),
            "support_thresholds": {"supported": WEAK_SUPPORT, "strong": STRONG_SUPPORT},
            "clades": [c.as_dict() for c in self.clades],
            "conflicting_clades": [c.as_dict() for c in self.conflicting_clades],
        }


def resample_columns(seqs: dict[str, str], rng: np.random.Generator) -> dict[str, str]:
    """Draw ``n_sites`` columns with replacement, keeping columns intact.

    Columns, not characters: the whole point is that a site is the unit of
    independent evidence. Resampling characters within a column would destroy
    the very correlations across taxa that carry the phylogenetic signal, and
    would produce uniformly high support for whatever the starting tree was.
    """
    names = list(seqs)
    n_sites = len(seqs[names[0]])
    idx = rng.integers(0, n_sites, size=n_sites)
    return {nm: "".join([seqs[nm][i] for i in idx]) for nm in names}


def validate_replicates(replicates: int) -> None:
    if not MIN_REPLICATES <= replicates <= MAX_REPLICATES:
        raise ValueError(
            f"replicates must be between {MIN_REPLICATES} and {MAX_REPLICATES}, "
            f"got {replicates}. Fewer than {MIN_REPLICATES} gives a support value "
            "whose own sampling error is larger than the differences it is being "
            "used to judge."
        )


def bootstrap_support(
    seqs: dict[str, str],
    ml_tree,
    model: str,
    replicates: int = DEFAULT_REPLICATES,
    seed: int = 1,
) -> BootstrapResult:
    """Resample, rebuild, and count how often each clade reappears."""
    validate_replicates(replicates)
    from .alignment import to_cogent3

    tips = set(seqs)
    ml_splits = canonical_splits(ml_tree, tips)
    counts: dict[Split, int] = dict.fromkeys(ml_splits, 0)

    rng = np.random.default_rng(seed)
    pt = piqtree()
    for _ in range(replicates):
        replicate = to_cogent3(resample_columns(seqs, rng))
        # The engine seed is held FIXED across replicates on purpose. The
        # randomness being measured is the resampling of sites, not the search
        # heuristic's starting point; varying both would fold search noise into
        # a number that is reported as sampling uncertainty.
        tree = pt.build_tree(replicate, model, rand_seed=seed)
        for split in canonical_splits(tree, tips):
            counts[split] = counts.get(split, 0) + 1

    ml_clades = [
        CladeSupport(sorted(s), counts[s] / replicates, True)
        for s in sorted(ml_splits, key=lambda s: (-counts[s], format_split(s)))
    ]
    conflicting = [
        CladeSupport(sorted(s), c / replicates, False)
        for s, c in sorted(counts.items(), key=lambda kv: -kv[1])
        if s not in ml_splits and c / replicates >= WEAK_SUPPORT
    ]

    values = [c.support for c in ml_clades]
    return BootstrapResult(
        replicates=replicates,
        clades=ml_clades,
        conflicting_clades=conflicting,
        n_ml_clades=len(ml_clades),
        n_supported=sum(1 for v in values if v >= WEAK_SUPPORT),
        n_strongly_supported=sum(1 for v in values if v >= STRONG_SUPPORT),
        mean_support=float(np.mean(values)) if values else 0.0,
        min_support=min(values) if values else None,
        seed=seed,
    )
