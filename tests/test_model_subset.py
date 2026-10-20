"""select_substitution_model(models=[...]): ranking within a named set."""

from __future__ import annotations

import pytest

from phylokit_mcp.server import select_substitution_model


@pytest.fixture(scope="module")
def hky_fasta():
    from phylokit_mcp.server import simulate_alignment

    tree = "((a:0.1,b:0.1):0.05,(c:0.1,d:0.1):0.05,(e:0.1,f:0.1):0.05);"
    return simulate_alignment(
        tree, model="HKY{4}+F{0.35,0.15,0.15,0.35}", length=1500, seed=3
    )["fasta"]


def test_the_subset_ranking_is_the_full_ranking_cut_down(hky_fasta):
    full = select_substitution_model(hky_fasta, criterion="BIC", top_n=1000)
    want = ["JC", "K2P", "F81+F", "HKY+F"]
    # top_n=1 must not cut a named set: each model in it was asked for.
    sub = select_substitution_model(hky_fasta, criterion="BIC", models=want, top_n=1)
    assert sub["restricted_to"] == want
    assert sorted(m["model"] for m in sub["ranking"]) == sorted(want)
    by_full = {m["model"]: m for m in full["ranking"]}
    order = sorted(want, key=lambda n: by_full[n]["bic"])
    assert [m["model"] for m in sub["ranking"]] == order
    assert sub["best_model"] == order[0] == "HKY+F"
    assert sub["ranking"][0]["delta"] == 0
    for m in sub["ranking"]:
        assert m["bic"] == by_full[m["model"]]["bic"]
    # Every per-criterion winner is from inside the set.
    assert set(sub["best_by_criterion"].values()) <= set(want)
    assert sub["n_models_compared"] == 4


def test_textbook_names_resolve_and_bare_f81_is_flagged(hky_fasta):
    sub = select_substitution_model(
        hky_fasta, criterion="BIC", models=["JC69", "K80", "F81", "HKY85"]
    )
    assert sub["resolved_names"] == {"JC69": "JC", "K80": "K2P", "HKY85": "HKY+F"}
    assert sub["restricted_to"] == ["JC", "K2P", "F81", "HKY+F"]
    same = [w for w in sub["warnings"] if w["code"] == "same_model_two_names"]
    assert len(same) == 1 and set(same[0]["models"]) == {"JC", "F81"}
    # Positive control: a set with no duplicate model raises no such warning.
    clean = select_substitution_model(
        hky_fasta, criterion="BIC", models=["JC", "F81+F"]
    )
    assert not [w for w in clean["warnings"] if w["code"] == "same_model_two_names"]


def test_an_unknown_name_is_refused_not_dropped(hky_fasta):
    assert (
        select_substitution_model(hky_fasta, models=["GTR+F"])["best_model"] == "GTR+F"
    )
    with pytest.raises(ValueError, match=r"did not score \['HKY99'\]"):
        select_substitution_model(hky_fasta, models=["HKY+F", "HKY99"])
    with pytest.raises(ValueError, match="models is empty"):
        select_substitution_model(hky_fasta, models=[])
