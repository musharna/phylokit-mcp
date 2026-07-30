"""Model selection, and the reason it returns a ranking rather than a name.

A caller handed one model name will state it as a finding. These tests pin the
behaviour that makes that statement checkable: the margin over the runner-up,
and whether the three criteria even agree.
"""

from __future__ import annotations

import math

import pytest

from phylokit_mcp.inference import CRITERIA, DELTA_INDISTINGUISHABLE, rank_models
from phylokit_mcp.server import (
    compare_trees,
    select_substitution_model,
    simulate_alignment,
)

from .conftest import TRUE_NEWICK


class FakeStat:
    """Constructed rather than simulated, so likelihood and parameter count can
    be varied INDEPENDENTLY — real model fits never let you hold one fixed."""

    def __init__(self, lnL: float, nfp: int) -> None:
        self.lnL = lnL
        self.nfp = nfp


def test_aic_penalises_parameters():
    """Equal fit, more parameters -> worse AIC. The whole point of the criterion."""
    ranked = rank_models(
        {"simple": FakeStat(-100.0, 5), "complex": FakeStat(-100.0, 9)}, 500
    )
    assert ranked[0].name == "simple"
    assert ranked[1].delta == pytest.approx(2 * (9 - 5))


def test_bic_penalises_parameters_more_heavily_than_aic():
    """At n=1000, ln(n) ~ 6.9 per parameter against AIC's 2.

    This is why the criteria disagree, and why the server reports all three.
    """
    # The 4 extra parameters must buy MORE than 4 log-likelihood units to win
    # under AIC (2 per parameter against 2 per lnL unit) and far more than that
    # under BIC. 5 units sits deliberately between the two thresholds; at 3, as
    # an earlier version of this test had it, AIC also prefers the simple model
    # and the test asserts nothing about the difference between the criteria.
    stats = {"simple": FakeStat(-100.0, 5), "complex": FakeStat(-95.0, 9)}
    by_aic = rank_models(stats, 1000, "AIC")[0].name
    by_bic = rank_models(stats, 1000, "BIC")[0].name
    assert by_aic == "complex"
    assert by_bic == "simple"


def test_aic_matches_its_definition():
    ranked = rank_models({"m": FakeStat(-250.0, 7)}, 400)
    assert ranked[0].aic == pytest.approx(2 * 7 - 2 * -250.0)
    assert ranked[0].bic == pytest.approx(7 * math.log(400) - 2 * -250.0)


def test_aicc_is_withheld_when_undefined():
    """n - k - 1 <= 0 makes the correction term negative or infinite.

    Returning a number there would be worse than returning nothing, because it
    would be compared against a well-defined one.
    """
    assert rank_models({"m": FakeStat(-10.0, 20)}, 21)[0].aicc is None
    assert rank_models({"m": FakeStat(-10.0, 5)}, 500)[0].aicc is not None


def test_delta_is_zero_for_the_winner_and_positive_after():
    ranked = rank_models(
        {"a": FakeStat(-100.0, 5), "b": FakeStat(-101.0, 5), "c": FakeStat(-110.0, 5)},
        500,
    )
    assert ranked[0].delta == 0.0
    assert ranked[1].delta > 0
    assert ranked[1].delta < ranked[2].delta


def test_entries_without_a_likelihood_are_skipped_not_scored_as_zero():
    ranked = rank_models({"good": FakeStat(-100.0, 5), "broken": object()}, 500)
    assert [m.name for m in ranked] == ["good"]


@pytest.mark.parametrize("criterion", CRITERIA)
def test_every_advertised_criterion_works(criterion: str):
    ranked = rank_models(
        {"a": FakeStat(-100.0, 5), "b": FakeStat(-140.0, 9)}, 500, criterion
    )
    assert ranked[0].name == "a"


def test_an_unknown_criterion_is_refused():
    with pytest.raises(ValueError, match="criterion must be one of"):
        rank_models({"a": FakeStat(-100.0, 5)}, 500, "AICx")


def test_model_selection_returns_the_runners_up_not_just_a_winner(fasta_easy):
    result = select_substitution_model(
        fasta=fasta_easy, criterion="AIC", seed=1, top_n=5
    )
    assert result["n_models_compared"] > 50
    assert len(result["ranking"]) == 5
    assert result["ranking"][0]["model"] == result["best_model"]
    assert result["ranking"][0]["delta"] == 0.0


def test_a_close_model_choice_is_flagged(fasta_easy):
    """Data simulated under JC do not single out one model.

    Measured: the AIC winner on this fixture is not JC, and several models fall
    within the conventional indistinguishability margin — so reporting the
    winner alone would overstate the comparison.
    """
    result = select_substitution_model(fasta=fasta_easy, criterion="AIC", seed=1)
    assert result["indistinguishable_from_best"], (
        "expected near-ties on JC-simulated data"
    )
    assert "model_choice_is_close" in {w["code"] for w in result["warnings"]}
    ties = [m for m in result["ranking"][1:] if m["delta"] <= DELTA_INDISTINGUISHABLE]
    assert ties


def test_criteria_agreement_is_reported(fasta_easy):
    result = select_substitution_model(fasta=fasta_easy, seed=1)
    assert set(result["best_by_criterion"]) == set(CRITERIA)
    expected = len(set(result["best_by_criterion"].values())) == 1
    assert result["criteria_agree"] is expected


def test_simulation_recovers_its_own_tree(easy_alignment):
    """The end-to-end positive control: truth in, truth out.

    If this fails, no result elsewhere in the suite can be trusted, because the
    engine cannot recover a topology it was handed.
    """
    from phylokit_mcp.server import infer_tree

    inferred = infer_tree(
        fasta=easy_alignment["fasta"], model="JC", replicates=20, seed=1
    )
    assert compare_trees(
        newick_a=inferred["newick"], newick_b=easy_alignment["true_newick"]
    )["identical_topology"]


def test_simulation_is_reproducible_under_a_fixed_seed():
    a = simulate_alignment(newick=TRUE_NEWICK, model="JC", length=200, seed=5)
    b = simulate_alignment(newick=TRUE_NEWICK, model="JC", length=200, seed=5)
    assert a["fasta"] == b["fasta"]


def test_simulation_varies_with_the_seed():
    a = simulate_alignment(newick=TRUE_NEWICK, model="JC", length=200, seed=5)
    b = simulate_alignment(newick=TRUE_NEWICK, model="JC", length=200, seed=6)
    assert a["fasta"] != b["fasta"]
