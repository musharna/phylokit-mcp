"""restrict_models called directly, on constructed statistics.

Through the server it is only reached after `select_model` has refused an
empty list and run model_finder, so its own empty check never runs there, and
the real fits never hold a likelihood fixed while the parameter count changes.
Constructed statistics vary the two independently.
"""

from __future__ import annotations

import pytest

from phylokit_mcp.inference import restrict_models


class Stat:
    def __init__(self, lnL: float, nfp: int) -> None:
        self.lnL = lnL
        self.nfp = nfp


SAME_MODEL_MESSAGE = (
    "'F81' and 'JC' have the same log-likelihood and parameter count: they are "
    "one model here. In IQ-TREE a bare 'F81' or 'HKY' has equal base frequencies "
    "(F81 = JC, HKY = K2P); the textbook F81 and HKY85 are 'F81+F' and 'HKY+F'."
)


def test_an_empty_list_is_refused_with_its_remedy():
    stats = {"JC": Stat(-100.0, 9)}
    assert restrict_models(stats, ["JC"])[0] == {"JC": stats["JC"]}
    with pytest.raises(ValueError) as exc:
        restrict_models(stats, [])
    assert str(exc.value) == "models is empty; omit it to rank every model."


def test_an_unknown_name_lists_the_base_names_that_were_scored():
    stats = {
        "JC": Stat(-110.0, 9),
        "HKY+F": Stat(-100.0, 13),
        "GTR+F+G4": Stat(-95.0, 18),
    }
    assert list(restrict_models(stats, ["HKY+F"])[0]) == ["HKY+F"]
    with pytest.raises(ValueError) as exc:
        restrict_models(stats, ["HKY+F", "TN93"])
    assert str(exc.value) == (
        "model_finder did not score ['TN93']. Use IQ-TREE names with their "
        "suffixes, e.g. 'HKY+F', 'F81+F', 'K2P', 'JC', 'GTR+F+G4'. Base names "
        "it scored: ['GTR', 'HKY', 'JC']."
    )


def test_same_model_needs_both_the_likelihood_and_the_parameter_count():
    stats = {
        "JC": Stat(-100.0, 9),
        "F81": Stat(-100.0, 9),  # JC under another name
        "K2P": Stat(-100.0, 10),  # same likelihood, one more parameter
        "TNe": Stat(-99.0, 9),  # same parameter count, better likelihood
    }
    _, _, warnings = restrict_models(stats, ["JC", "F81", "K2P", "TNe"])
    assert warnings == [
        {
            "code": "same_model_two_names",
            "message": SAME_MODEL_MESSAGE,
            "models": ["JC", "F81"],
        }
    ]
    # Neither half alone is enough.
    assert restrict_models(stats, ["JC", "K2P"])[2] == []
    assert restrict_models(stats, ["JC", "TNe"])[2] == []


def test_an_entry_without_statistics_is_kept_not_a_crash():
    # rank_models skips such an entry (test_models.py); restricting to it first
    # must not raise before the ranking gets the chance.
    stats = {"JC": Stat(-100.0, 9), "broken": object()}
    subset, resolved, warnings = restrict_models(stats, ["JC", "broken"])
    assert list(subset) == ["JC", "broken"]
    assert resolved == {}
    assert warnings == []
