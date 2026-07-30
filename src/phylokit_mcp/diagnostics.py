"""Advisories — the things a caller would otherwise state as findings.

Each keys on a MEASUREMENT, never on a name or a setting. The distinction
matters: "the model was chosen by AIC" is a setting and says nothing about
whether the choice was close, whereas the delta to the runner-up does.

Nothing here withholds a result. Every value returned by this server is already
accompanied by its own uncertainty, so a caller can see when an answer is thin.
These explain rather than censor.
"""

from __future__ import annotations

from .alignment import AlignmentStats
from .bootstrap import STRONG_SUPPORT, WEAK_SUPPORT, BootstrapResult
from .inference import DELTA_INDISTINGUISHABLE

# Below this share of resolved clades the topology is mostly unresolved and
# should not be read as a phylogeny.
POORLY_RESOLVED = 0.7
# Parsimony-informative sites per taxon, below which support will be low for
# reasons of information content rather than genuine conflict.
THIN_INFORMATION = 10


def unsupported_clades(result: BootstrapResult) -> dict | None:
    weak = [c.as_dict() for c in result.clades if c.support < WEAK_SUPPORT]
    if not weak:
        return None
    return {
        "code": "unsupported_clades",
        "message": (
            f"{len(weak)} of {result.n_ml_clades} clades in this tree fall below "
            f"{WEAK_SUPPORT:.2f} bootstrap support. Those groupings are not "
            "distinguishable from alternatives given these data — the tree draws "
            "them with the same confidence as the rest, but they are not results."
        ),
        "clades": weak,
    }


def poorly_resolved(result: BootstrapResult) -> dict | None:
    if result.fraction_resolved >= POORLY_RESOLVED:
        return None
    return {
        "code": "poorly_resolved",
        "message": (
            f"Only {result.fraction_resolved:.0%} of clades clear {WEAK_SUPPORT:.2f} "
            "support. Treat this as a partially resolved tree; collapsing the "
            "unsupported nodes to a polytomy is the honest representation."
        ),
        "fraction_resolved": round(result.fraction_resolved, 3),
    }


def conflicting_support(result: BootstrapResult) -> dict | None:
    """Well-supported clades that the ML tree does not contain.

    This is the advisory with no equivalent in a support-annotated Newick
    string: a conflicting clade has nowhere to attach in the winning topology,
    so the standard output format cannot express it and it goes unreported.
    """
    if not result.conflicting_clades:
        return None
    return {
        "code": "conflicting_support",
        "message": (
            f"{len(result.conflicting_clades)} clade(s) reached at least "
            f"{WEAK_SUPPORT:.2f} support across replicates but are ABSENT from the "
            "maximum-likelihood tree. The data support a grouping the reported "
            "topology does not show, so the single best tree is not the whole story."
        ),
        "clades": [c.as_dict() for c in result.conflicting_clades],
    }


def thin_information(
    stats: AlignmentStats, result: BootstrapResult | None = None
) -> dict | None:
    """Flag an alignment thin enough that low support is an information problem.

    This EXPLAINS a measured shortfall; it does not predict one. When a bootstrap
    result is in hand, support has been measured, and a proxy that disagrees with
    the measurement is simply wrong — an early version fired here on a tree whose
    every clade sat at 1.00 support, because 7.3 informative sites per taxon was
    below a threshold picked a priori. Suppressing the forecast whenever the
    outcome is known removes that whole class of contradiction; no choice of
    threshold does, since the proxy and the measurement are different quantities.

    Without a result — `simulate_alignment`, where no tree has been built — it is
    the only signal available, so it fires as a forecast and says so.
    """
    per_taxon = stats.n_parsimony_informative / max(1, stats.n_taxa)
    if per_taxon >= THIN_INFORMATION:
        return None
    if result is not None and (result.min_support or 0.0) >= WEAK_SUPPORT:
        # Every clade cleared the line, so the information was evidently
        # sufficient. Keyed on the weakest clade rather than on the RESOLVED
        # FRACTION deliberately: what this advisory explains is an unsupported
        # clade, and a tree can be 75% resolved while still having one.
        return None
    tail = (
        "At least one clade is unsupported, and this is the likely reason: too "
        "little signal, rather than biological conflict."
        if result is not None
        else "Expect low support for reasons of information content."
    )
    return {
        "code": "thin_information",
        "message": (
            f"{stats.n_parsimony_informative} parsimony-informative sites across "
            f"{stats.n_taxa} taxa ({per_taxon:.1f} per taxon). Alignment LENGTH is "
            f"{stats.n_sites}, but length is not evidence — only variable sites "
            f"shared by at least two taxa carry topological signal. {tail}"
        ),
        "n_parsimony_informative": stats.n_parsimony_informative,
        "n_sites": stats.n_sites,
    }


def duplicate_sequences(stats: AlignmentStats) -> dict | None:
    if not stats.duplicate_sequences:
        return None
    return {
        "code": "duplicate_sequences",
        "message": (
            "Some sequences are identical. Their branching order relative to one "
            "another is arbitrary — the tree will show one, drawn no differently "
            "from a resolved node — and their branch lengths will be zero."
        ),
        "groups": stats.duplicate_sequences,
    }


def gappy_alignment(stats: AlignmentStats) -> dict | None:
    if stats.fraction_gaps < 0.5:
        return None
    return {
        "code": "gappy_alignment",
        "message": (
            f"{stats.fraction_gaps:.0%} of the alignment is gaps or missing data. "
            "Sites are still resampled whole by the bootstrap, so support values "
            "remain interpretable, but the effective information is far below the "
            f"nominal {stats.n_sites} sites."
        ),
        "fraction_gaps": stats.fraction_gaps,
    }


def model_choice_is_close(selection: dict) -> dict | None:
    ties = selection.get("indistinguishable_from_best") or []
    if not ties:
        return None
    return {
        "code": "model_choice_is_close",
        "message": (
            f"{len(ties)} model(s) score within {DELTA_INDISTINGUISHABLE} of the "
            f"winner ({selection['best_model']}), which is the conventional line "
            "for 'indistinguishable'. Reporting the best model as THE model "
            "overstates what the comparison established."
        ),
        "models": ties,
    }


def criteria_disagree(selection: dict) -> dict | None:
    if selection.get("criteria_agree", True):
        return None
    return {
        "code": "criteria_disagree",
        "message": (
            "AIC, AICc and BIC do not choose the same model: "
            f"{selection['best_by_criterion']}. BIC penalises parameters more "
            "heavily, so this pattern usually means the extra parameters buy a "
            "real but small likelihood gain. Which to prefer is a judgement, not "
            "a result the data settle."
        ),
        "best_by_criterion": selection["best_by_criterion"],
    }


def threads_not_pinned(pinned: bool) -> dict | None:
    if pinned:
        return None
    return {
        "code": "threads_not_pinned",
        "message": (
            "OMP_NUM_THREADS is not 1. Likelihood sums are accumulated in "
            "thread-completion order and floating-point addition is not "
            "associative, so results may not reproduce even with a fixed seed."
        ),
    }


def collect(
    stats: AlignmentStats | None = None,
    result: BootstrapResult | None = None,
    selection: dict | None = None,
    pinned: bool = True,
) -> list[dict]:
    out: list[dict | None] = []
    if stats is not None:
        out += [
            thin_information(stats, result),
            duplicate_sequences(stats),
            gappy_alignment(stats),
        ]
    if result is not None:
        out += [
            poorly_resolved(result),
            unsupported_clades(result),
            conflicting_support(result),
        ]
    if selection is not None:
        out += [model_choice_is_close(selection), criteria_disagree(selection)]
    out.append(threads_not_pinned(pinned))
    return [w for w in out if w]


__all__ = ["POORLY_RESOLVED", "STRONG_SUPPORT", "THIN_INFORMATION", "collect"]
