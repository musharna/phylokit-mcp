"""Tree building and model selection.

Model selection returns the runners-up, not just the winner. IQ-TREE's
`model_finder` reports a single best model per criterion, and a caller handed one
name will state it as a finding. Measured on a 400-site alignment simulated under
JC, the AIC winner is F81 — a model the data were not generated under — and the
margin over JC is small enough that the ranking would not survive a different
seed. A winner without its margin is a claim the numbers do not support.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .engine import piqtree

# Conventional readings of an information-criterion gap. Burnham & Anderson
# (2002) treat models within ~2 of the best as essentially indistinguishable.
DELTA_INDISTINGUISHABLE = 2.0
CRITERIA = ("AIC", "AICc", "BIC")


@dataclass(frozen=True)
class ModelScore:
    name: str
    log_likelihood: float
    n_free_parameters: int
    aic: float
    aicc: float | None
    bic: float
    delta: float

    def as_dict(self) -> dict:
        return {
            "model": self.name,
            "log_likelihood": round(self.log_likelihood, 4),
            "n_free_parameters": self.n_free_parameters,
            "aic": round(self.aic, 3),
            "aicc": None if self.aicc is None else round(self.aicc, 3),
            "bic": round(self.bic, 3),
            "delta": round(self.delta, 3),
        }


def _information_criteria(
    lnl: float, k: int, n_sites: int
) -> tuple[float, float | None, float]:
    aic = 2 * k - 2 * lnl
    bic = k * math.log(n_sites) - 2 * lnl
    # AICc is undefined once the parameter count reaches the sample size; a
    # negative or infinite correction term is worse than declining to report it.
    aicc = aic + (2 * k * (k + 1)) / (n_sites - k - 1) if n_sites - k - 1 > 0 else None
    return aic, aicc, bic


def rank_models(
    model_stats: dict, n_sites: int, criterion: str = "AIC"
) -> list[ModelScore]:
    """Every candidate model scored and ranked, best first, with delta to best."""
    if criterion not in CRITERIA:
        raise ValueError(f"criterion must be one of {CRITERIA}, got {criterion!r}")

    scored: list[ModelScore] = []
    for name, stat in model_stats.items():
        lnl = getattr(stat, "lnL", None)
        k = getattr(stat, "nfp", None)
        if lnl is None or k is None:
            continue
        aic, aicc, bic = _information_criteria(float(lnl), int(k), n_sites)
        scored.append(
            ModelScore(str(name), float(lnl), int(k), aic, aicc, bic, delta=0.0)
        )

    def key(m: ModelScore) -> float:
        value = {
            "AIC": m.aic,
            "AICc": m.aicc if m.aicc is not None else m.aic,
            "BIC": m.bic,
        }
        return value[criterion]

    scored.sort(key=key)
    if not scored:
        return []
    best = key(scored[0])
    return [
        ModelScore(
            m.name,
            m.log_likelihood,
            m.n_free_parameters,
            m.aic,
            m.aicc,
            m.bic,
            key(m) - best,
        )
        for m in scored
    ]


def select_model(
    seqs: dict[str, str], criterion: str = "AIC", seed: int = 1, top_n: int = 5
) -> dict:
    """Rank substitution models, and say how much the winner actually won by."""
    from .alignment import to_cogent3

    n_sites = len(next(iter(seqs.values())))
    result = piqtree().model_finder(to_cogent3(seqs), rand_seed=seed)
    ranked = rank_models(result.model_stats, n_sites, criterion)
    if not ranked:
        raise RuntimeError("model_finder returned no scoreable models.")

    ties = [m for m in ranked[1:] if m.delta <= DELTA_INDISTINGUISHABLE]
    best_by = {c: rank_models(result.model_stats, n_sites, c)[0].name for c in CRITERIA}
    return {
        "criterion": criterion,
        "best_model": ranked[0].name,
        "n_models_compared": len(ranked),
        "ranking": [m.as_dict() for m in ranked[: max(1, top_n)]],
        "indistinguishable_from_best": [m.name for m in ties],
        "best_by_criterion": best_by,
        "criteria_agree": len(set(best_by.values())) == 1,
        "seed": seed,
    }


def build_ml_tree(seqs: dict[str, str], model: str, seed: int = 1):
    from .alignment import to_cogent3

    return piqtree().build_tree(to_cogent3(seqs), model, rand_seed=seed)


def tree_log_likelihood(tree) -> float | None:
    params = getattr(tree, "params", None) or {}
    lnl = params.get("lnL")
    return None if lnl is None else float(lnl)
