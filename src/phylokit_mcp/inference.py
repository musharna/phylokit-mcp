"""Tree building and model selection.

Model selection returns the runners-up, not just the winner. IQ-TREE's
`model_finder` reports a single best model per criterion, and a caller handed one
name will state it as a finding. Measured on the test suite's 300-site alignment
simulated under JC (piqtree 0.8.3, seed 1), the AIC winner is TPM2u (tied with
TPM2) — a model the data were not generated under — and JC trails it by only 0.7 AIC units. A
winner without its margin is a claim the numbers do not support.
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


def _stats_by_name(model_stats: dict) -> dict[str, object]:
    """One entry per model name.

    piqtree 0.8.3's `ModelFinderResult.model_stats` holds every model under its
    string name, then adds the best AIC, AICc and BIC models again under `Model`
    object keys, which do not compare equal to those strings. Iterated as is,
    the best model appeared up to three times and was listed as
    indistinguishable from itself. The extra keys carry the same statistics as
    the string key; two entries under one name that DISAGREE would be a
    different, unexplained state, so that is refused rather than resolved.
    """
    by_name: dict[str, object] = {}
    for key, stat in model_stats.items():
        name = str(key)
        if name not in by_name:
            by_name[name] = stat
            continue
        seen = by_name[name]
        fields = ("lnL", "nfp")
        if [getattr(seen, f, None) for f in fields] != [
            getattr(stat, f, None) for f in fields
        ]:
            raise RuntimeError(
                f"model_finder reported {name!r} twice with different "
                f"statistics: {seen!r} and {stat!r}."
            )
    return by_name


def rank_models(
    model_stats: dict, n_sites: int, criterion: str = "AIC"
) -> list[ModelScore]:
    """Every candidate model scored and ranked, best first, with delta to best."""
    if criterion not in CRITERIA:
        raise ValueError(f"criterion must be one of {CRITERIA}, got {criterion!r}")

    scored: list[ModelScore] = []
    for name, stat in _stats_by_name(model_stats).items():
        lnl = getattr(stat, "lnL", None)
        k = getattr(stat, "nfp", None)
        if lnl is None or k is None:
            continue
        aic, aicc, bic = _information_criteria(float(lnl), int(k), n_sites)
        scored.append(ModelScore(name, float(lnl), int(k), aic, aicc, bic, delta=0.0))

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


# Textbook names IQ-TREE does not use, and the IQ-TREE model each one means.
# HKY85 and F81 estimate base frequencies; IQ-TREE's bare "HKY" and "F81" are
# scored by model_finder with EQUAL frequencies (measured: bare F81 has the
# log-likelihood and parameter count of JC), so HKY85 is "HKY+F". Bare "F81" is
# an IQ-TREE name in its own right and is not rewritten; the same-likelihood
# warning in `restrict_models` is what flags it.
TEXTBOOK_NAMES = {"JC69": "JC", "K80": "K2P", "HKY85": "HKY+F"}


def restrict_models(
    model_stats: dict, models: list[str]
) -> tuple[dict, dict[str, str], list[dict]]:
    """`model_stats` cut down to the named models.

    Returns the subset, the textbook names that were resolved to IQ-TREE names,
    and warnings for requested models that are the same model under two names
    (identical log-likelihood and parameter count). An unknown name is refused
    with the names that were compared, never dropped: a silently shorter set
    would rank the survivors as if the missing model had lost.
    """
    if not models:
        raise ValueError("models is empty; omit it to rank every model.")
    by_name = _stats_by_name(model_stats)
    resolved = {
        m: TEXTBOOK_NAMES[m] for m in models if m not in by_name and m in TEXTBOOK_NAMES
    }
    wanted = [resolved.get(m, m) for m in models]
    unknown = [m for m, w in zip(models, wanted) if w not in by_name]
    if unknown:
        bases = sorted({n.split("+")[0] for n in by_name})
        raise ValueError(
            f"model_finder did not score {unknown}. Use IQ-TREE names with their "
            f"suffixes, e.g. 'HKY+F', 'F81+F', 'K2P', 'JC', 'GTR+F+G4'. Base names "
            f"it scored: {bases}."
        )
    subset = {w: by_name[w] for w in dict.fromkeys(wanted)}
    warnings: list[dict] = []
    seen: dict[tuple, str] = {}
    for name, stat in subset.items():
        sig = (getattr(stat, "lnL", None), getattr(stat, "nfp", None))
        if sig in seen:
            warnings.append(
                {
                    "code": "same_model_two_names",
                    "message": (
                        f"{name!r} and {seen[sig]!r} have the same log-likelihood and "
                        "parameter count: they are one model here. In IQ-TREE a bare "
                        "'F81' or 'HKY' has equal base frequencies (F81 = JC, HKY = "
                        "K2P); the textbook F81 and HKY85 are 'F81+F' and 'HKY+F'."
                    ),
                    "models": [seen[sig], name],
                }
            )
        else:
            seen[sig] = name
    return subset, resolved, warnings


def model_finder_stats(aln, seed: int) -> dict:
    """ModelFinder's per-model statistics, keyed by IQ-TREE's model names.

    piqtree 0.8.3's `model_finder` builds a ModelFinderResult whose constructor
    parses the three best-model names with `make_model`, and that parser does
    not know every name ModelFinder reports. On the agent benchmark's q23 and
    q26 alignments the AIC winner is bare "K3Pu", and the whole call failed
    with "Unknown substitution model: 'K3Pu'", every statistic lost with it.
    The statistics are all in IQ-TREE's raw output, so they are read from it
    here, with piqtree's own call and parser for each value. This uses
    piqtree's private `_model_finder` module; `tests/test_models.py` checks
    that it returns exactly what the public call returns where that call works.
    """
    piqtree()  # the engine module, imported under engine.py's thread pin
    # PyYAML is piqtree's dependency, and this parses piqtree's output exactly
    # as its own model_finder does; it goes if that private module changes.
    import yaml  # type: ignore[import-untyped]
    from piqtree.iqtree import _model_finder as mf
    from piqtree.util import process_rand_seed_nonzero

    names = aln.names
    raw = yaml.safe_load(
        mf.iq_model_finder(
            names,
            [str(s) for s in aln.iter_seqs(names)],
            process_rand_seed_nonzero(seed),
            "",
            "",
            "",
            1,
            "",
        )
    )
    return {
        key: mf.ModelResultValue.from_string(val)
        for key, val in raw.items()
        if not key.startswith(("best_", "initTree")) and isinstance(val, str)
    }


def select_model(
    seqs: dict[str, str],
    criterion: str = "AIC",
    seed: int = 1,
    top_n: int = 5,
    moltype: str = "dna",
    models: list[str] | None = None,
) -> dict:
    """Rank substitution models, and say how much the winner actually won by.

    With `models`, the ranking, deltas, ties and per-criterion winners are all
    computed within that set, from the same model_finder run.
    """
    from .alignment import to_cogent3

    if models is not None and not models:
        raise ValueError("models is empty; omit it to rank every model.")
    n_sites = len(next(iter(seqs.values())))
    stats = model_finder_stats(to_cogent3(seqs, moltype), seed)
    resolved: dict[str, str] = {}
    warnings: list[dict] = []
    if models is not None:
        stats, resolved, warnings = restrict_models(stats, models)
    ranked = rank_models(stats, n_sites, criterion)
    if not ranked:
        raise RuntimeError("model_finder returned no scoreable models.")

    ties = [m for m in ranked[1:] if m.delta <= DELTA_INDISTINGUISHABLE]
    best_by = {c: rank_models(stats, n_sites, c)[0].name for c in CRITERIA}
    out = {
        "criterion": criterion,
        "best_model": ranked[0].name,
        "n_models_compared": len(ranked),
        # A named set is returned whole: it was asked for model by model.
        "ranking": [
            m.as_dict()
            for m in (ranked if models is not None else ranked[: max(1, top_n)])
        ],
        "indistinguishable_from_best": [m.name for m in ties],
        "best_by_criterion": best_by,
        "criteria_agree": len(set(best_by.values())) == 1,
        "seed": seed,
    }
    if models is not None:
        out["restricted_to"] = list(stats)
        out["resolved_names"] = resolved
        out["model_warnings"] = warnings
    return out


# piqtree's own classification of a substitution model, mapped onto this
# server's sequence_type vocabulary.
_MODEL_KIND_TO_MOLTYPE = {"nucleotide": "dna", "protein": "protein"}


def model_moltype(model: str) -> str:
    """The molecule type a substitution model is defined over, per piqtree.

    Parsed with piqtree's `make_model`, the same parser the engine call uses, so
    a name it does not know is refused here with its message rather than later.
    """
    from piqtree.model import make_model

    kind = make_model(model).submod_type.model_type()
    try:
        return _MODEL_KIND_TO_MOLTYPE[kind]
    except KeyError:
        raise ValueError(
            f"Model {model!r} is a {kind} model; this server supports "
            f"{sorted(_MODEL_KIND_TO_MOLTYPE.values())} only."
        ) from None


def require_model_for(model: str, moltype: str) -> None:
    """Refuse a model defined over the other molecule type, before the engine.

    IQ-TREE's own answer to `LG` on DNA, or `GTR` on protein, is "File not found
    LG" -- a RuntimeError the server could only report as a crash, for what is a
    mismatch in the caller's arguments.
    """
    wanted = model_moltype(model)
    if wanted != moltype:
        raise ValueError(
            f"Model {model!r} is a {wanted} model but sequence_type is "
            f"{moltype!r}. Use a {moltype} model"
            + (" such as 'GTR+G' or 'HKY'" if moltype == "dna" else " such as 'LG+G'")
            + f", or pass sequence_type={wanted!r} if the data are {wanted}."
        )


def build_ml_tree(
    seqs: dict[str, str], model: str, seed: int = 1, moltype: str = "dna"
):
    from .alignment import to_cogent3

    return piqtree().build_tree(to_cogent3(seqs, moltype), model, rand_seed=seed)


def tree_log_likelihood(tree) -> float | None:
    params = getattr(tree, "params", None) or {}
    lnl = params.get("lnL")
    return None if lnl is None else float(lnl)
