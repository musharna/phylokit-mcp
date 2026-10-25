"""Nonparametric bootstrap support (Felsenstein 1985, Evolution 39:783-791).

Two methods. The nonparametric bootstrap is computed here, one full ML search
per replicate. UFBoot2 (`ufboot_support`) is IQ-TREE's own, read back from the
tree piqtree returns.

An earlier version of this docstring said piqtree 0.8.3 did not expose UFBoot
support. It does, as each internal node's `.support` (cogent3 parses a numeric
internal label into that attribute, not into the name). The measurement that
said otherwise looked at node names and `params` only. Re-measured 2026-10-25.

Why both. The nonparametric loop costs one ML search per replicate: on the
field benchmark's 89-taxon, 435-site alignment one replicate took 25-93 s on
one thread, so 100 replicates outran a 30-minute client budget every time.
UFBoot2's 1000 replicates took 83 s on the same alignment. The nonparametric
method keeps the side benefit below, which UFBoot read back this way cannot
offer: the full split frequency distribution, including well-supported clades
that are NOT in the maximum-likelihood tree. A conflicting clade at 0.70 is the single
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
from typing import Any

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

# UFBoot2 (Hoang et al. 2018, Mol Biol Evol 35:518-522, doi:10.1093/molbev/msx281).
# Its values are not on the nonparametric scale: the authors' own reading is
# that a clade is supported at 95 or above, so 0.95 is the one line for both
# "supported" and "strong". IQ-TREE refuses fewer than 1000 replicates.
UFBOOT_SUPPORT = 0.95
UFBOOT_MIN_REPLICATES = 1000
UFBOOT_MAX_REPLICATES = 10000

NONPARAMETRIC = "nonparametric bootstrap over alignment columns"
UFBOOT = "ultrafast bootstrap (UFBoot2, IQ-TREE), read from the ML tree"
SUPPORT_METHODS = ("bootstrap", "ufboot")


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
    method: str = NONPARAMETRIC
    supported_at: float = WEAK_SUPPORT
    strong_at: float = STRONG_SUPPORT

    @property
    def fraction_resolved(self) -> float:
        """Share of the ML tree's clades that clear the conventional 0.70 line.

        This is the headline number: a tree whose clades are mostly unsupported
        is a drawing, not a result. "Clear" is the method's own line
        (`supported_at`), 0.70 nonparametric or 0.95 UFBoot2.
        """
        if self.n_ml_clades == 0:
            return 0.0
        return self.n_supported / self.n_ml_clades

    def as_dict(self) -> dict:
        return {
            "replicates": self.replicates,
            "method": self.method,
            "seed": self.seed,
            "n_clades": self.n_ml_clades,
            "n_supported": self.n_supported,
            "n_strongly_supported": self.n_strongly_supported,
            "fraction_resolved": round(self.fraction_resolved, 3),
            "mean_support": round(self.mean_support, 3),
            "min_support": None
            if self.min_support is None
            else round(self.min_support, 3),
            "support_thresholds": {
                "supported": self.supported_at,
                "strong": self.strong_at,
            },
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


def validate_support_method(method: str) -> None:
    if method not in SUPPORT_METHODS:
        raise ValueError(
            f"support_method must be one of {list(SUPPORT_METHODS)}, got {method!r}."
        )


def default_replicates(method: str) -> int:
    """The method's default count. Not a check: `validate_replicates` is."""
    return UFBOOT_MIN_REPLICATES if method == "ufboot" else DEFAULT_REPLICATES


def validate_replicates(replicates: int, method: str = "bootstrap") -> None:
    validate_support_method(method)
    if method == "ufboot":
        if not UFBOOT_MIN_REPLICATES <= replicates <= UFBOOT_MAX_REPLICATES:
            raise ValueError(
                f"ufboot replicates must be between {UFBOOT_MIN_REPLICATES} and "
                f"{UFBOOT_MAX_REPLICATES}, got {replicates}. IQ-TREE refuses UFBoot2 "
                f"below {UFBOOT_MIN_REPLICATES}; it is cheap enough that the floor "
                "costs little."
            )
        return
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
    moltype: str = "dna",
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
        replicate = to_cogent3(resample_columns(seqs, rng), moltype)
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


def _summarise(
    clades: list[CladeSupport],
    conflicting: list[CladeSupport],
    replicates: int,
    seed: int,
    method: str,
    supported_at: float,
    strong_at: float,
) -> BootstrapResult:
    values = [c.support for c in clades]
    return BootstrapResult(
        replicates=replicates,
        clades=clades,
        conflicting_clades=conflicting,
        n_ml_clades=len(clades),
        n_supported=sum(1 for v in values if v >= supported_at),
        n_strongly_supported=sum(1 for v in values if v >= strong_at),
        mean_support=float(np.mean(values)) if values else 0.0,
        min_support=min(values) if values else None,
        seed=seed,
        method=method,
        supported_at=supported_at,
        strong_at=strong_at,
    )


def ufboot_tree(
    seqs: dict[str, str],
    model: str,
    replicates: int = UFBOOT_MIN_REPLICATES,
    seed: int = 1,
    moltype: str = "dna",
) -> tuple[Any, BootstrapResult]:
    """The ML tree and its UFBoot2 support, from ONE IQ-TREE run.

    UFBoot2 is not a separate pass over a finished tree: IQ-TREE samples the
    replicates during its own ML search, so the tree returned here IS the ML
    tree, and it replaces `build_ml_tree` for this method. No conflicting-clade
    list: IQ-TREE annotates only the clades of its tree, so `conflicting_clades`
    is empty by construction, and the `method` string says so.
    """
    validate_replicates(replicates, "ufboot")
    from .alignment import to_cogent3

    tree = piqtree().build_tree(
        to_cogent3(seqs, moltype),
        model,
        rand_seed=seed,
        bootstrap_replicates=replicates,
    )
    tips = set(seqs)
    ref = min(tips)
    support: dict[Split, float] = {}
    for edge in tree.get_edge_vector(include_root=False):
        if edge.is_tip():
            continue
        below = frozenset(t.name for t in edge.tips())
        side = below if ref not in below else frozenset(tips) - below
        if not 2 <= len(side) <= len(tips) - 2:
            continue
        value = getattr(edge, "support", None)
        if value is None:
            # Fail loud: a silent 0 or a skipped clade would read as "unsupported".
            raise RuntimeError(
                f"IQ-TREE returned no UFBoot support for clade {format_split(side)} "
                f"after {replicates} replicates (piqtree {getattr(piqtree(), '__version__', '?')}). "
                "The engine's output format has changed; the support cannot be reported."
            )
        support[side] = float(value) / 100.0
    if set(support) != canonical_splits(tree, tips):
        raise RuntimeError("UFBoot support did not cover every clade of the ML tree.")
    clades = [
        CladeSupport(sorted(s), v, True)
        for s, v in sorted(
            support.items(), key=lambda kv: (-kv[1], format_split(kv[0]))
        )
    ]
    result = _summarise(
        clades, [], replicates, seed, UFBOOT, UFBOOT_SUPPORT, UFBOOT_SUPPORT
    )
    return tree, result
