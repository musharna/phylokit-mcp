"""The agent benchmark's scorer and answer key.

Scorer tests pair each wrong or declined reply with the right one for the same
question, so a scorer that called everything right (or wrong) fails here. Key
tests recompute answers from the text the question shows, with the benchmark's
own independent code, so a hand edit to questions.json or a drift between
writer and key cannot pass.
"""

from __future__ import annotations

import ast
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

BENCH = Path(__file__).resolve().parent.parent / "benchmarks" / "agent_qa"
sys.path.insert(0, str(BENCH))

import truth as T
from score import extract, score

MODEL_ALIASES = {
    "JC69": ["JC69", "JC"],
    "K80": ["K80", "K2P"],
    "F81": ["F81", "F81+F", "JC+F"],
    "HKY85": ["HKY85", "HKY", "HKY+F", "K2P+F"],
}
SISTER = {
    "category": "sister_taxon",
    "truth": {"taxon": "t4", "taxa": ["t1", "t2", "t3", "t4"]},
}
RF = {"category": "rf_distance", "truth": {"rf": 6, "max_rf": 14}}
HOMOLOG = {
    "category": "homologous_position",
    "truth": {"position": 57, "moltype": "dna"},
}
MODEL = {
    "category": "model_selection",
    "truth": {"model": "HKY85", "aliases": MODEL_ALIASES},
}
DIST = {
    "category": "pairwise_distance",
    "truth": {"distance": 0.40, "tolerance": 0.075},
}
INFORMATIVE = {"category": "informative_sites", "truth": {"count": 42}}


def _v(q: dict, reply: str) -> str:
    return score(q, reply)["verdict"]


def test_the_last_answer_line_decides_and_markdown_around_it_is_ignored() -> None:
    reply = "ANSWER: 4\nChecking again...\n**ANSWER:** `6`"
    assert extract(reply) == "6"
    assert _v(RF, reply) == "right"
    assert _v(RF, "**ANSWER:** 4") == "wrong"


def test_unknown_is_declined_and_a_reply_without_the_line_is_neither() -> None:
    assert _v(RF, "Not sure.\nANSWER: UNKNOWN") == "declined"
    assert _v(RF, "It is 6.") == "no_answer"
    assert _v(RF, "It is 6.\nANSWER: 6") == "right"


def test_an_integer_is_read_whole_and_never_out_of_a_decimal() -> None:
    assert _v(RF, "ANSWER: RF = 6") == "right"
    assert _v(RF, "ANSWER: 6 (normalised 0.43)") == "right"
    assert _v(RF, "ANSWER: 6.") == "right"
    # "0.6" holds a 6 but is not the integer 6; "16" is not 6.
    assert _v(RF, "ANSWER: 0.6") == "wrong"
    assert _v(RF, "ANSWER: 16") == "wrong"
    assert _v(INFORMATIVE, "ANSWER: 42 sites") == "right"
    assert _v(INFORMATIVE, "ANSWER: 41") == "wrong"
    assert _v(HOMOLOG, "ANSWER: 57") == "right"
    assert _v(HOMOLOG, "ANSWER: 58") == "wrong"
    assert _v(HOMOLOG, "ANSWER: gap") == "wrong"


def test_a_taxon_is_matched_by_name_and_a_wrong_one_is_flagged() -> None:
    assert _v(SISTER, "ANSWER: t4") == "right"
    assert _v(SISTER, "ANSWER: T4") == "right"
    assert _v(SISTER, "ANSWER: taxon t4") == "right"
    assert _v(SISTER, "ANSWER: t3") == "wrong"
    assert score(SISTER, "ANSWER: t3")["wrong_id"] is True  # a real taxon, not the one
    assert score(SISTER, "ANSWER: t9")["wrong_id"] is False  # not in the alignment
    assert score(SISTER, "ANSWER: t4")["wrong_id"] is False  # right, so not a wrong id
    # "t44" is not t4: the match is on whole words.
    assert _v(SISTER, "ANSWER: t44") == "wrong"


def test_every_name_of_the_right_model_counts_and_a_wrong_candidate_is_flagged() -> (
    None
):
    for name in ("HKY85", "HKY", "hky+F", "K2P+F", "HKY85 (lowest BIC)"):
        assert _v(MODEL, f"ANSWER: {name}") == "right", name
    assert _v(MODEL, "ANSWER: K80") == "wrong"
    assert score(MODEL, "ANSWER: F81")["wrong_id"] is True
    assert score(MODEL, "ANSWER: GTR")["wrong_id"] is False  # not one of the four
    # IQ-TREE's "+F" changes the model: K2P+F is HKY, so it is not K80.
    k80 = {
        "category": "model_selection",
        "truth": {"model": "K80", "aliases": MODEL_ALIASES},
    }
    assert _v(k80, "ANSWER: K2P") == "right"
    assert _v(k80, "ANSWER: K2P+F") == "wrong"


def test_a_distance_is_right_inside_the_tolerance_and_wrong_outside_it() -> None:
    assert _v(DIST, "ANSWER: 0.40") == "right"
    assert _v(DIST, "ANSWER: 0.474") == "right"
    assert _v(DIST, "ANSWER: 0.326") == "right"
    assert _v(DIST, "ANSWER: ~0.33 substitutions/site") == "right"
    assert _v(DIST, "ANSWER: 0.476") == "wrong"
    assert _v(DIST, "ANSWER: 0.324") == "wrong"
    assert _v(DIST, "ANSWER: 40") == "wrong"


def test_an_unknown_category_is_refused_rather_than_scored() -> None:
    with pytest.raises(ValueError, match="unknown category"):
        score({"category": "nope", "truth": {}}, "ANSWER: 1")
    assert _v(RF, "ANSWER: 6") == "right"


# --- the key ---------------------------------------------------------------------


def test_the_key_builder_never_imports_the_server() -> None:
    """The point of an independent key. Read from the AST, so a mention in a
    docstring is not an import and an import cannot hide in one."""
    for name in ("truth.py", "build_questions.py", "score.py"):
        tree = ast.parse((BENCH / name).read_text())
        mods: set[str] = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                mods.update(a.name for a in n.names)
            elif isinstance(n, ast.ImportFrom) and n.module:
                mods.add(n.module)
        assert not {m for m in mods if m.startswith("phylokit_mcp")}, name
        # Positive control: the walk does see this file's real imports.
        assert "re" in mods if name == "score.py" else "numpy" in mods, name


def test_jc_transition_probabilities_match_the_closed_form() -> None:
    pi = np.full(4, 0.25)
    p = T.transition_matrix(T.rate_matrix(1.0, pi), pi, 0.3)
    same = 0.25 + 0.75 * math.exp(-4 * 0.3 / 3)
    assert p[0, 0] == pytest.approx(same, abs=1e-12)
    assert p[0, 1] == pytest.approx((1 - same) / 3, abs=1e-12)
    # kappa reaches the matrix: transitions (A<->G) beat transversions.
    k = T.transition_matrix(T.rate_matrix(6.0, pi), pi, 0.3)
    assert k[0, 2] > 2 * k[0, 1]


def test_every_newick_spelling_of_a_tree_has_its_splits() -> None:
    rng = np.random.default_rng(3)
    adj = T.random_tree([f"t{i}" for i in range(1, 10)], rng, (0.05, 0.1), (0.05, 0.2))
    want = T.splits(adj)
    assert len(want) == 9 - 3
    for root_on_edge in (False, True):
        for lengths in (False, True):
            nw = T.to_newick(adj, rng, lengths=lengths, root_on_edge=root_on_edge)
            assert T.splits(T.parse_newick(nw)) == want, nw
    other = T.random_tree(
        [f"t{i}" for i in range(1, 10)], rng, (0.05, 0.1), (0.05, 0.2)
    )
    assert T.robinson_foulds(T.splits(other), want) > 0


def test_neighbour_joining_recovers_a_clear_tree_and_not_a_shuffled_one() -> None:
    rng = np.random.default_rng(5)
    taxa = [f"t{i}" for i in range(1, 9)]
    adj = T.random_tree(taxa, rng, (0.1, 0.2), (0.05, 0.1))
    seqs = T.simulate(adj, 1000, rng)
    names, d = T.jc_distance_matrix(seqs)
    assert T.robinson_foulds(T.splits(T.neighbor_joining(names, d)), T.splits(adj)) == 0
    cols = [list(c) for c in zip(*seqs.values(), strict=True)]
    for c in cols:
        rng.shuffle(c)
    shuffled = {t: "".join(c[i] for c in cols) for i, t in enumerate(seqs)}
    names, d = T.jc_distance_matrix(shuffled)
    assert T.robinson_foulds(T.splits(T.neighbor_joining(names, d)), T.splits(adj)) > 0


def test_the_indel_simulator_keeps_every_residue_in_its_column() -> None:
    rng = np.random.default_rng(11)
    adj = T.random_tree([f"t{i}" for i in range(1, 7)], rng, (0.05, 0.1), (0.05, 0.1))
    seqs, aln = T.simulate_with_indels(adj, 150, rng, alphabet=T.AMINO, indel_rate=0.1)
    assert len({len(s) for s in seqs.values()}) > 1  # indels happened
    for t, cols in aln.items():
        assert [c for c in cols if c is not None] == list(range(1, len(seqs[t]) + 1))
    assert all(any(aln[t][c] is not None for t in aln) for c in range(len(aln["t1"])))


def test_parsimony_informative_counts_two_states_twice_each() -> None:
    seqs = {"a": "AAAC", "b": "AACC", "c": "GACT", "d": "GATT"}
    # site 1: A,A,G,G -> yes; 2: all A -> no; 3: A,C,C,T (only C twice) -> no;
    # 4: C,C,T,T -> yes
    assert T.parsimony_informative(seqs) == 2


def _questions() -> list[dict]:
    return json.loads((BENCH / "questions.json").read_text())["questions"]


def test_the_frozen_key_matches_the_data_shown_in_each_question() -> None:
    """Recompute what can be recomputed from the question text alone."""
    checked = {"rf_distance": 0, "informative_sites": 0, "homologous_position": 0}
    for q in _questions():
        if q["category"] == "rf_distance":
            nw_a, nw_b = q["files"]["trees.nwk"].split()
            assert f"Tree A: {nw_a}\nTree B: {nw_b}" in q["question"]
            got = T.robinson_foulds(
                T.splits(T.parse_newick(nw_a)), T.splits(T.parse_newick(nw_b))
            )
            assert got == q["truth"]["rf"], q["id"]
        elif q["category"] == "informative_sites":
            text = q["files"]["alignment.fasta"]
            assert text in q["question"]
            seqs = {
                block.split("\n", 1)[0]: block.split("\n", 1)[1].replace("\n", "")
                for block in text.strip(">\n").split("\n>")
            }
            assert T.parsimony_informative(seqs) == q["truth"]["count"], q["id"]
        elif q["category"] == "homologous_position":
            assert q["files"]["sequences.fasta"] in q["question"]
        else:
            continue
        checked[q["category"]] += 1
    assert checked == {
        "rf_distance": 8,
        "informative_sites": 6,
        "homologous_position": 8,
    }


def _right_and_wrong(category: str, t: dict) -> tuple[str, str]:
    if category == "sister_taxon":
        return t["taxon"], next(x for x in t["taxa"] if x != t["taxon"])
    if category == "rf_distance":
        return str(t["rf"]), str(t["rf"] + 2)
    if category == "homologous_position":
        return str(t["position"]), str(t["position"] + 1)
    if category == "model_selection":
        return t["model"], "JC69" if t["model"] != "JC69" else "K80"
    if category == "pairwise_distance":
        return str(t["distance"]), f"{t['distance'] + 1.01 * t['tolerance']:.4f}"
    return str(t["count"]), str(t["count"] - 1)


def test_the_key_has_every_category_and_unique_ids() -> None:
    qs = _questions()
    counts: dict[str, int] = {}
    for q in qs:
        counts[q["category"]] = counts.get(q["category"], 0) + 1
    assert len({q["id"] for q in qs}) == len(qs) == 46
    assert counts == {
        "sister_taxon": 8,
        "rf_distance": 8,
        "homologous_position": 8,
        "model_selection": 8,
        "pairwise_distance": 8,
        "informative_sites": 6,
    }
    # Every stored truth scores right against itself, and a nearby wrong value
    # scores wrong, so no question is right whatever is answered.
    for q in qs:
        right, wrong = _right_and_wrong(q["category"], q["truth"])
        assert _v(q, f"ANSWER: {right}") == "right", q["id"]
        assert _v(q, f"ANSWER: {wrong}") == "wrong", q["id"]
