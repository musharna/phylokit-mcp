"""The held-out P2/P3 key for the confirmatory run (PREREG amendment
2026-10-25): built by the pilot's builder from draws 2000+, checked before any
arm runs. A failed control here stops the run.

The pilot's own invariants (bands, tolerances, truths recomputed from the
generator) are re-run on this key by pointing ``test_agent_qa_field``'s key
reader at it. Counts the pilot only observed (no +G, same relative position)
are reported, not asserted.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from . import test_agent_qa_field as F

CONFIRM = F.BENCH / "field_2026-10-25_confirm"
FAMILIES = {"supported_clades": "P2", "align_tree_homology": "P3"}


def _key() -> dict[str, Any]:
    return json.loads((CONFIRM / "questions.json").read_text())


@pytest.fixture
def confirm_key(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """``test_agent_qa_field``'s helpers, reading this key instead of the pilot's."""
    key = _key()
    monkeypatch.setattr(F, "_key", lambda: key)
    return key


def test_the_held_out_key_has_ten_per_family_and_the_oracle_agrees(
    confirm_key: dict[str, Any],
) -> None:
    qs = confirm_key["questions"]
    counts: dict[str, int] = {}
    for q in qs:
        counts[q["category"]] = counts.get(q["category"], 0) + 1
    assert counts == dict.fromkeys(FAMILIES, 10)
    assert len({q["id"] for q in qs}) == 20
    assert [F._v(q, F.controls.key_value(q)) for q in qs] == ["right"] * 20
    assert [F._v(q, F.controls.wrong_value(q)) for q in qs] == ["wrong"] * 20
    assert F.controls.oracle_check(qs) == []
    assert confirm_key["alt_web"] is False


def test_the_held_out_draws_start_at_2000_and_share_nothing_with_the_pilot(
    confirm_key: dict[str, Any],
) -> None:
    draws = confirm_key["draws"]
    assert set(draws) == set(FAMILIES)
    for fam, d in draws.items():
        assert d["first_draw"] == 2000, fam
        assert d["kept"] == 10, fam
        assert d["drawn"] == d["kept"] + d["discarded"], fam
        assert sum(d["discard_reasons"].values()) == d["discarded"], fam
    pilot = json.loads((F.FIELD / "questions.json").read_text())["questions"]
    held = confirm_key["questions"]
    for field in ("alisim_seed", "alignment_sha256"):
        seen = {q["generator"].get(field) for q in pilot} - {None}
        new = {q["generator"].get(field) for q in held} - {None}
        assert new, field  # positive control: the field is there to compare
        assert not seen & new, field
    assert not {q["question"] for q in pilot} & {q["question"] for q in held}


def test_the_reference_passes_every_held_out_question_and_the_shortcuts_fail(
    confirm_key: dict[str, Any],
) -> None:
    v = F._control_verdicts()
    for cat in FAMILIES:
        assert v[cat]["reference pipeline (seed 2)"] == ["right"] * 10, cat
    assert v["supported_clades"]["ML tree, every clade reported"] == ["wrong"] * 10
    assert v["align_tree_homology"]["no alignment (pad at end) + ML"] == ["wrong"] * 10
    rel = v["align_tree_homology"]["same relative position + reference tree"]
    print(f"same relative position + reference tree: {rel.count('wrong')}/10 wrong")


@pytest.mark.parametrize(
    "check",
    [
        F.test_no_field_question_prescribes_a_substitution_model,
        F.test_each_field_truth_is_recomputed_from_the_generator_not_trusted,
        F.test_the_p2_band_is_the_reference_minus_three_binomial_se,
        F.test_the_p3_tolerance_is_mean_plus_max_of_three_sd_and_one_split,
        F.test_p3_admits_one_split_beyond_the_reference_and_not_two,
        F.test_the_builders_tolerance_step_reproduces_the_stored_bands,
    ],
    ids=lambda f: f.__name__.removeprefix("test_"),
)
def test_the_pilots_key_invariants_hold_on_the_held_out_key(
    check: Any, confirm_key: dict[str, Any]
) -> None:
    check()


def test_the_key_reader_really_is_redirected(confirm_key: dict[str, Any]) -> None:
    """Without this, every check above could be re-testing the pilot's key."""
    assert F._questions() == confirm_key["questions"]
    pilot = json.loads((F.FIELD / "questions.json").read_text())
    assert F._questions() != pilot["questions"]
    assert Path(CONFIRM / "questions.json").exists()
