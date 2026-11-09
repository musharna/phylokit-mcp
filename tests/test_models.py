"""Model selection, and the reason it returns a ranking rather than a name.

A caller handed one model name will state it as a finding. These tests pin the
behaviour that makes that statement checkable: the margin over the runner-up,
and whether the three criteria even agree.
"""

from __future__ import annotations

import math

import pytest

from phylokit_mcp.diagnostics import criteria_disagree, model_choice_is_close
from phylokit_mcp.inference import (
    CRITERIA,
    DELTA_INDISTINGUISHABLE,
    ModelScore,
    _stats_by_name,
    rank_models,
)
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


def test_each_model_is_ranked_once_and_never_its_own_runner_up(fasta_easy):
    """piqtree adds the best AIC/AICc/BIC models to `model_stats` a second time,
    under `Model` keys that do not equal their string names. Ranked as is, the
    winner appeared up to three times and was listed as indistinguishable from
    itself (measured on this fixture: TPM2u x3, F81 x2).

    Real engine output, not a constructed dict: the duplication is piqtree's,
    and a constructed dict would only encode a belief about what piqtree returns.
    """
    result = select_substitution_model(
        fasta=fasta_easy, criterion="AIC", seed=1, top_n=1000
    )
    names = [m["model"] for m in result["ranking"]]
    assert len(names) == len(set(names)), sorted(
        n for n in set(names) if names.count(n) > 1
    )
    assert result["n_models_compared"] == len(names)
    assert result["best_model"] not in result["indistinguishable_from_best"]
    # Positive control: deduplicating must not drop the real near-ties. JC
    # generated these data and must still be ranked, and still be close.
    assert "JC" in names
    assert "JC" in result["indistinguishable_from_best"]
    for criterion in CRITERIA:
        assert result["best_by_criterion"][criterion] in names


def test_a_name_reported_twice_with_different_statistics_is_refused():
    """Duplicate keys are collapsed only when they carry the same statistics.

    Two entries under one name that disagree are not something piqtree is
    known to produce; picking one of them would hide it.
    """

    class Key:
        def __init__(self, name: str) -> None:
            self.name = name

        def __str__(self) -> str:
            return self.name

    agreeing = {"m": FakeStat(-100.0, 5), Key("m"): FakeStat(-100.0, 5)}
    assert [m.name for m in rank_models(agreeing, 500)] == ["m"]
    clashing = {"m": FakeStat(-100.0, 5), Key("m"): FakeStat(-90.0, 5)}
    with pytest.raises(RuntimeError, match="twice with different statistics"):
        rank_models(clashing, 500)


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


def _benchmark_alignment(qid: str) -> str:
    import json
    from pathlib import Path

    key = Path(__file__).parent.parent / "benchmarks" / "agent_qa" / "questions.json"
    (q,) = [q for q in json.loads(key.read_text())["questions"] if q["id"] == qid]
    return q["files"]["alignment.fasta"]


def test_raw_model_stats_equal_piqtrees_own_where_its_call_works(fasta_easy):
    from phylokit_mcp.alignment import parse_fasta, to_cogent3
    from phylokit_mcp.engine import piqtree
    from phylokit_mcp.inference import model_finder_stats

    aln = to_cogent3(parse_fasta(fasta_easy))
    public = piqtree().model_finder(aln, rand_seed=1).model_stats
    ours = model_finder_stats(aln, 1)
    assert ours == {k: v for k, v in public.items() if isinstance(k, str)}
    assert len(ours) > 100


def test_a_best_model_piqtree_cannot_parse_does_not_lose_the_ranking():
    """q23's AIC winner is a name piqtree 0.8.3's make_model rejects; its own
    model_finder fails on it, and ours must still rank every model."""
    from phylokit_mcp.alignment import parse_fasta, to_cogent3
    from phylokit_mcp.engine import piqtree

    fasta = _benchmark_alignment("q23")
    with pytest.raises(ValueError, match="Unknown substitution model"):
        piqtree().model_finder(to_cogent3(parse_fasta(fasta)), rand_seed=1)
    out = select_substitution_model(
        fasta, criterion="BIC", models=["JC", "K2P", "F81+F", "HKY+F"]
    )
    assert out["best_model"] == "F81+F"  # the generating model (questions.json)
    assert len(out["ranking"]) == 4


# Known answers. The tests above check orderings and use `approx` on one model,
# so the nightly (2026-11-08) left 54 mutants of the arithmetic alive: the AICc
# correction term and its guard, the AICc fallback, the skip rules, the
# duplicate check's fields, and every key and rounding in `as_dict`.


def test_aicc_is_its_definition_and_is_reported_down_to_one_spare_site():
    """AICc = AIC + 2k(k+1)/(n-k-1), exactly. It is defined while n-k-1 >= 1:
    at n = k+2 it is a number, at n = k+1 it is withheld."""
    [m] = rank_models({"m": FakeStat(-100.0, 5)}, 500)
    assert (m.aic, m.aicc, m.bic) == (210.0, 210.0 + 60 / 494, 5 * math.log(500) + 200)
    assert m.n_free_parameters == 5
    assert rank_models({"m": FakeStat(-100.0, 5)}, 7)[0].aicc == 210.0 + 60 / 1
    assert rank_models({"m": FakeStat(-100.0, 5)}, 6)[0].aicc is None


def test_aicc_ranks_by_aicc_and_falls_back_to_aic_only_where_undefined():
    """At n=20, b wins on AIC (203 < 204) and loses on AICc (205.67 > 204.71):
    the small-sample correction is what AICc is for. Where AICc is undefined
    (k = n-1) the model is ranked by its AIC instead."""
    stats = {"a": FakeStat(-100.0, 2), "b": FakeStat(-97.5, 4)}
    assert [(m.name, m.delta) for m in rank_models(stats, 20, "AIC")] == [
        ("b", 0.0),
        ("a", 1.0),
    ]
    by_aicc = rank_models(stats, 20, "AICc")
    assert [m.name for m in by_aicc] == ["a", "b"]
    assert by_aicc[1].delta == (203.0 + 40 / 15) - (204.0 + 12 / 17)

    fallback = rank_models(
        {"undefined": FakeStat(-100.0, 9), "defined": FakeStat(-108.0, 2)}, 10, "AICc"
    )
    assert [(m.name, m.aicc) for m in fallback] == [
        ("undefined", None),
        ("defined", 220.0 + 12 / 7),
    ]
    assert fallback[1].delta == (220.0 + 12 / 7) - 218.0


def test_an_entry_missing_either_statistic_is_skipped_wherever_it_sits():
    """Missing lnL or missing parameter count: either one makes the entry
    unscorable. A broken entry listed first must not end the scan."""

    class OnlyLikelihood:
        lnL = -90.0

    ranked = rank_models(
        {"broken": object(), "no_k": OnlyLikelihood(), "good": FakeStat(-100.0, 5)},
        500,
    )
    assert [m.name for m in ranked] == ["good"]


def test_duplicates_must_agree_on_the_parameter_count_too():
    """Same name, same lnL, different parameter count is a disagreement. Two
    entries with no statistics at all agree (nothing differs), and collapse."""

    class Key:
        def __init__(self, name: str) -> None:
            self.name = name

        def __str__(self) -> str:
            return self.name

    blank, blank_again = object(), object()
    assert _stats_by_name({"m": blank, Key("m"): blank_again}) == {"m": blank}
    with pytest.raises(RuntimeError) as exc:
        _stats_by_name({"m": FakeStat(-100.0, 5), Key("m"): FakeStat(-100.0, 6)})
    assert str(exc.value).startswith(
        "model_finder reported 'm' twice with different statistics: "
    )


def test_a_model_score_reports_every_field_under_its_name_at_its_precision():
    """Every value carries more digits than its rounding keeps, so one more or
    one fewer decimal place, or none, changes the output."""
    score = ModelScore(
        "m", -100.123456789, 5, 210.2469135, 210.368371, 231.3194567, 1.2345678
    )
    assert score.as_dict() == {
        "model": "m",
        "log_likelihood": -100.1235,
        "n_free_parameters": 5,
        "aic": 210.247,
        "aicc": 210.368,
        "bic": 231.319,
        "delta": 1.235,
    }
    undefined = ModelScore("m", -100.0, 20, 240.0, None, 300.0, 0.0)
    assert undefined.as_dict()["aicc"] is None


def test_criteria_disagree_speaks_only_on_a_recorded_disagreement():
    """The advisory fires on `criteria_agree: False` and on nothing else.

    A selection that does not say whether the criteria agree is not evidence
    that they disagree: the advisory stays silent rather than inventing a
    disagreement and then failing to name it.
    """
    best_by = {"AIC": "GTR+G", "AICc": "HKY+G", "BIC": "HKY+G"}
    assert criteria_disagree({"criteria_agree": True, "best_by_criterion": {}}) is None
    assert criteria_disagree({}) is None
    assert criteria_disagree(
        {"criteria_agree": False, "best_by_criterion": best_by}
    ) == {
        "code": "criteria_disagree",
        "message": (
            "AIC, AICc and BIC do not choose the same model: "
            "{'AIC': 'GTR+G', 'AICc': 'HKY+G', 'BIC': 'HKY+G'}. BIC penalises "
            "parameters more heavily, so this pattern usually means the extra "
            "parameters buy a real but small likelihood gain. Which to prefer is "
            "a judgement, not a result the data settle."
        ),
        "best_by_criterion": best_by,
    }


def test_model_choice_is_close_names_the_winner_and_every_near_tie():
    """No near-ties, however the absence is spelled, is silence; any near-tie
    is reported with its count, the margin, the winner and the tied models."""
    for quiet in (
        {},
        {"indistinguishable_from_best": None},
        {"indistinguishable_from_best": []},
    ):
        assert model_choice_is_close({"best_model": "GTR", **quiet}) is None
    ties = ["HKY", "JC"]
    assert model_choice_is_close(
        {"best_model": "GTR", "indistinguishable_from_best": ties}
    ) == {
        "code": "model_choice_is_close",
        "message": (
            "2 model(s) score within 2.0 of the winner (GTR), which is the "
            "conventional line for 'indistinguishable'. Reporting the best "
            "model as THE model overstates what the comparison established."
        ),
        "models": ties,
    }
