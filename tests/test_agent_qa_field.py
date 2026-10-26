"""The 2026-10-24 field families (P1 long-branch attraction, P2 reportable
clades at 95%, P3 unaligned indel-rich DNA): scorer, key, controls, and the
harness's web-off alt arm.

Every refusal sits beside the legitimate case in the same test: a key answer
that scores right beside a wrong one that scores wrong, a shortcut that fails
beside the reference that passes, an alt arm without the web beside the
default arm that has it.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest

BENCH = Path(__file__).resolve().parent.parent / "benchmarks" / "agent_qa"
FIELD = BENCH / "field_2026-10-24"
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(FIELD))

import build_field as B
import controls
import run
import truth as T
from score import score

CATEGORIES = {"lba_sister": "P1", "supported_clades": "P2", "align_tree_homology": "P3"}


def _key() -> dict[str, Any]:
    return json.loads((FIELD / "questions.json").read_text())


def _questions(category: str | None = None) -> list[dict[str, Any]]:
    qs = _key()["questions"]
    return [q for q in qs if category is None or q["category"] == category]


def _v(q: dict[str, Any], value: str) -> str:
    return score(q, controls.reply(value))["verdict"]


# --- scorer --------------------------------------------------------------------

# A binary 8-taxon tree: splits ab, cd, abcd, ef, efg (n - 3 = 5).
TRUE8 = "(((a,b),(c,d)),((e,f),g),h);"
P2Q = {
    "category": "supported_clades",
    "truth": {
        "newick": TRUE8,
        "taxa": list("abcdefgh"),
        "n_true_splits": 5,
        "precision_min": 0.9,
        "recall_min": 0.5,
    },
}
P3Q = {
    "category": "align_tree_homology",
    "truth": {
        "position": 57,
        "newick": TRUE8,
        "taxa": list("abcdefgh"),
        "nrf_max": 0.2,
    },
}


def test_supported_clades_is_scored_on_precision_and_recall_of_the_reported_splits() -> (
    None
):
    # Every true split, in another spelling and rooting: right.
    assert _v(P2Q, "(h,((g,(f,e)),((d,c),(b,a))));") == "right"
    # Three of five true splits and nothing false: recall 0.6, precision 1.
    got = score(P2Q, controls.reply("(((a,b),(c,d)),e,f,g,h);"))
    assert got["verdict"] == "right"
    assert got["detail"] == {"parsed": True, "n_reported": 3, "precision": 1.0,
                             "recall": 0.6}  # fmt: skip
    # Two of five: recall 0.4 < 0.5, wrong.
    assert _v(P2Q, "((a,b),(c,d),e,f,g,h);") == "wrong"
    # One false split (gh) among five reported: precision 0.8 < 0.9, wrong.
    got = score(P2Q, controls.reply("(((a,b),(c,d)),((e,f),(g,h)));"))
    assert got["verdict"] == "wrong" and got["detail"]["precision"] == 0.8
    # The star tree reports nothing: recall 0.
    assert _v(P2Q, "(a,b,c,d,e,f,g,h);") == "wrong"
    # A tree missing a taxon, or with an extra one, or no tree at all: wrong.
    assert _v(P2Q, "(((a,b),(c,d)),((e,f),g));") == "wrong"
    assert _v(P2Q, "(((a,b),(c,d)),((e,f),g),h,x);") == "wrong"
    assert _v(P2Q, "a b c") == "wrong"
    assert score(P2Q, "ANSWER: UNKNOWN")["verdict"] == "declined"


def test_align_tree_homology_needs_both_the_exact_residue_and_a_tree_in_tolerance() -> (
    None
):
    one_off = "(((a,c),(b,d)),((e,f),g),h);"  # ab, cd -> ac, bd: RF 4 of 10
    assert _v(P3Q, f"57 | {TRUE8}") == "right"
    # Tree first, residue second, is read the same way.
    assert _v(P3Q, f"{TRUE8} | 57") == "right"
    got = score(P3Q, controls.reply(f"58 | {TRUE8}"))
    assert got["verdict"] == "wrong" and got["detail"]["tree_right"] is True
    assert got["detail"]["position_right"] is False
    got = score(P3Q, controls.reply(f"57 | {one_off}"))
    assert got["verdict"] == "wrong" and got["detail"]["nrf"] == 0.4
    assert got["detail"]["position_right"] is True
    assert _v(P3Q, "57") == "wrong"  # no tree
    assert _v(P3Q, f"57 {TRUE8}") == "wrong"  # no separator


def test_lba_sister_reads_a_taxon_name_like_sister_taxon() -> None:
    q = {"category": "lba_sister", "truth": {"taxon": "t4", "taxa": ["t1", "t4", "t9"]}}
    assert _v(q, "taxon t4") == "right"
    got = score(q, controls.reply("t9"))
    assert got["verdict"] == "wrong" and got["wrong_id"] is True


# --- the frozen key -----------------------------------------------------------------


def test_the_field_key_has_ten_questions_per_family_and_the_oracle_agrees() -> None:
    qs = _questions()
    counts: dict[str, int] = {}
    for q in qs:
        counts[q["category"]] = counts.get(q["category"], 0) + 1
    assert counts == dict.fromkeys(CATEGORIES, 10)
    assert len({q["id"] for q in qs}) == 30
    right = [_v(q, controls.key_value(q)) for q in qs]
    wrong = [_v(q, controls.wrong_value(q)) for q in qs]
    assert right == ["right"] * 30
    assert wrong == ["wrong"] * 30
    assert controls.oracle_check(qs) == []


def test_the_oracle_check_names_a_field_question_whose_key_was_broken() -> None:
    qs = json.loads(json.dumps(_questions()))
    p2 = next(q for q in qs if q["category"] == "supported_clades")
    p2["truth"]["recall_min"] = 1.01  # no reply can reach it, the key's own included
    p3 = next(q for q in qs if q["category"] == "align_tree_homology")
    p3["truth"]["nrf_max"] = -0.01
    bad = controls.oracle_check(qs)
    assert [b.split(":")[0] for b in bad] == [p2["id"], p3["id"]]
    assert controls.oracle_check(_questions()) == []


def test_no_field_question_prescribes_a_substitution_model() -> None:
    names = re.compile(
        r"\b(JC\d*|K2P|K80|F81|HKY\d*|TN93|TrN|GTR|LG|WAG|gamma|invariant|"
        r"substitution model|model)\b|\+G|\+I",
        re.IGNORECASE,
    )
    for q in _questions():
        assert not names.search(q["question"] + q["answer_format"]), q["id"]
    # Positive control: the same rule finds the model the old set prescribes.
    old = json.loads((BENCH / "questions.json").read_text())["questions"]
    assert names.search(old[0]["question"])


def test_each_field_truth_is_recomputed_from_the_generator_not_trusted() -> None:
    """P1's answer is the asked taxon's cherry partner in the generating tree;
    P2's and P3's truth tree is the generating tree; every data file is the
    one whose sha256 the builder recorded."""
    for q in _questions():
        gen = T.parse_newick(q["generator"]["newick"])
        if q["category"] == "lba_sister":
            (asked,) = re.findall(r"cherry\) with (t\d+)\?", q["question"])
            assert B.cherry_partner(gen, asked) == q["truth"]["taxon"], q["id"]
            assert asked in q["generator"]["long_tips"], q["id"]
            assert B.cherry_partner(gen, q["truth"]["taxon"]) == asked
            long_a, long_b = q["generator"]["long_tips"]
            assert B.cherry_partner(gen, long_a) != long_b, q["id"]  # not sisters
        else:
            truth_tree = T.parse_newick(q["truth"]["newick"])
            assert T.splits(truth_tree) == T.splits(gen), q["id"]
        for name, text in q["files"].items():
            assert q["generator"]["artifacts_sha256"][name] == B.sha256(text), q["id"]
            assert f"`{name}`" in q["question"], q["id"]
            assert text.splitlines()[1] not in q["question"], q["id"]
        n = len(q["truth"]["taxa"])
        assert 3 < n <= 200
        seqs = B.read_fasta(next(iter(q["files"].values())))
        assert sorted(seqs) == q["truth"]["taxa"]
        assert max(len(s) for s in seqs.values()) <= 100_000  # the server's caps


def test_the_draw_bookkeeping_adds_up() -> None:
    draws = _key()["draws"]
    assert set(draws) == set(CATEGORIES)
    for fam, d in draws.items():
        assert d["kept"] == 10, fam
        assert d["drawn"] == d["kept"] + d["discarded"], fam
        assert sum(d["discard_reasons"].values()) == d["discarded"], fam


# --- controls: reference pipeline and shortcuts ----------------------------------------


def _control_verdicts() -> dict[str, dict[str, list[str]]]:
    out: dict[str, dict[str, list[str]]] = {}
    for r in controls.field_controls(_questions()):
        out.setdefault(r["category"], {}).setdefault(r["control"], []).append(
            r["verdict"]
        )
    return out


def test_the_reference_pipeline_passes_every_question_and_each_shortcut_fails() -> None:
    """Positive control: the reference pipeline re-run at another seed (never
    the server) is right on all 30. Shortcut controls: the JC tree under
    long-branch attraction, the unbootstrapped ML tree reporting every
    clade, and skipping the alignment, each fail where the family says."""
    v = _control_verdicts()
    for cat in CATEGORIES:
        assert v[cat]["reference pipeline (seed 2)"] == ["right"] * 10, cat
    assert v["lba_sister"]["JC ML tree"] == ["wrong"] * 10
    # Dropping +G alone (GTR, no rate variation) still recovers the cherry on
    # 2 of 10 (p1q02, p1q06): JC's equal rates and frequencies add to it.
    assert v["lba_sister"]["GTR ML tree, no +G"].count("wrong") == 8
    assert v["supported_clades"]["ML tree, every clade reported"] == ["wrong"] * 10
    assert v["align_tree_homology"]["no alignment (pad at end) + ML"] == ["wrong"] * 10
    # Not the named shortcut, a cruder one, reported beside it: the same
    # relative position lands on the right residue once (p3q08, offset -18).
    rel = v["align_tree_homology"]["same relative position + reference tree"]
    assert rel.count("wrong") == 9


def test_the_p1_shortcut_fails_by_attraction_not_by_chance() -> None:
    """The JC tree puts the asked long branch with the OTHER long branch on
    every kept draw: the failure is long-branch attraction."""
    for q in _questions("lba_sister"):
        assert q["check"]["jc_pairs_the_long_branches"] is True, q["id"]
        assert q["check"]["partner_of_asked"]["jc"] == q["generator"]["long_tips"][1]


def test_the_p2_band_is_the_reference_minus_three_binomial_se() -> None:
    for q in _questions("supported_clades"):
        ref = q["check"]["reference"]
        n_rep, n_true = ref["n_reported"], q["truth"]["n_true_splits"]
        hits = ref["hits"]
        assert hits == round(ref["precision"] * n_rep) == round(ref["recall"] * n_true)
        p_min = hits / n_rep - 3 * B.binomial_se(hits, n_rep)
        r_min = hits / n_true - 3 * B.binomial_se(hits, n_true)
        assert q["truth"]["precision_min"] == pytest.approx(p_min, abs=1e-12), q["id"]
        assert q["truth"]["recall_min"] == pytest.approx(r_min, abs=1e-12), q["id"]
        assert 0 < q["truth"]["recall_min"] < ref["recall"]
        # The shortcut's failure is on precision: it reports more, not less.
        ml = q["check"]["ml_all_clades"]
        assert ml["recall"] >= ref["recall"] and ml["precision"] < p_min, q["id"]


def test_the_p3_tolerance_is_mean_plus_max_of_three_sd_and_one_split() -> None:
    import statistics

    for q in _questions("align_tree_homology"):
        rf = q["check"]["reference_rf"]
        assert len(rf) == 20
        nrf = [k / (2 * (len(q["truth"]["taxa"]) - 3)) for k in rf]
        one_split = 2 / (2 * len(q["truth"]["taxa"]) - 6)
        tol = statistics.mean(nrf) + max(3 * statistics.stdev(nrf), one_split)
        assert q["truth"]["nrf_max"] == pytest.approx(tol, abs=1e-12), q["id"]
        assert q["check"]["nrf_one_split"] == pytest.approx(one_split, abs=1e-15)
        # The reference's own nRF is inside its tolerance: the band edge was
        # once rounded below it, and the seed-2 reference failed its own key.
        assert max(nrf) <= q["truth"]["nrf_max"] + 1e-9, q["id"]
        assert q["check"]["padded_nrf"] > q["truth"]["nrf_max"], q["id"]
        assert abs(q["check"]["offset"]) >= B.HOMOLOGY_MIN_OFFSET


def _tree_at_rf(newick: str, target: int, seed: int) -> str:
    """A tree whose RF to ``newick`` is exactly ``target``, by seeded NNIs."""
    import numpy as np
    from build_questions import nni

    rng = np.random.default_rng(seed)
    true = T.splits(T.parse_newick(newick))
    for _ in range(200):
        adj = T.parse_newick(newick)
        for _ in range(target // 2):
            nni(adj, rng)
        if T.robinson_foulds(T.splits(adj), true) == target:
            return T.to_newick(adj, rng, lengths=False)
    raise AssertionError(f"no tree at RF {target}")


def test_p3_admits_one_split_beyond_the_reference_and_not_two() -> None:
    """The floor, scored through the real scorer on every P3 question: a tree
    one split worse than the reference pipeline's is right, two splits worse
    is wrong, both with the true residue (so only the tree part decides)."""
    for q in _questions("align_tree_homology"):
        (k,) = set(q["check"]["reference_rf"])  # all 20 runs agree
        pos, truth_tree = q["truth"]["position"], q["truth"]["newick"]
        one_worse = _tree_at_rf(truth_tree, k + 2, seed=1)
        two_worse = _tree_at_rf(truth_tree, k + 4, seed=2)
        assert _v(q, f"{pos} | {one_worse}") == "right", q["id"]
        assert _v(q, f"{pos} | {two_worse}") == "wrong", q["id"]
        assert q["check"]["nrf_floor_binds"] is True, q["id"]


# --- harness: question set flag, web-off alt, dry run -------------------------------


def test_the_field_set_runs_the_alt_arm_without_the_web(tmp_path: Path) -> None:
    key = run.load_key(FIELD / "questions.json")
    assert key["alt_web"] is False
    assert run.load_key(BENCH / "questions.json").get("alt_web", True) is True
    tmp = tmp_path / "cwd"
    tmp.mkdir()
    flags, _ = run.arm_setup("alt", tmp, web=False)
    assert flags[flags.index("--tools") + 1] == "Bash,Read"
    assert flags[flags.index("--allowedTools") + 1] == "Bash,Read"
    # Positive control: the default sets keep the web.
    flags, _ = run.arm_setup("alt", tmp, web=True)
    assert flags[flags.index("--tools") + 1] == "Bash,Read,WebSearch,WebFetch"


def _init(tools: list[str]) -> list[dict[str, Any]]:
    own = "/x/tmp/phyloqa-alt-aaaa1111/cwd"
    return [
        {"type": "system", "subtype": "init", "cwd": own, "tools": tools,
         "mcp_servers": []},
        {"type": "assistant", "message": {"content": [
            {"type": "tool_use", "name": "Bash", "input": {"command": "ls"}}]}},
        {"type": "result", "subtype": "success", "is_error": False},
    ]  # fmt: skip


def test_a_web_tool_in_a_web_off_alt_run_makes_it_invalid() -> None:
    assert run.validity("alt", _init(["Bash", "Read"]), 0, web=False)["valid"]
    web = run.validity(
        "alt", _init(["Bash", "Read", "WebSearch", "WebFetch"]), 0, False
    )
    assert not web["valid"] and web["invalid_reasons"][0].startswith("tools")
    # The same transcript is valid for a web-on alt arm, and a web-off one is not.
    assert run.validity("alt", _init(["Bash", "Read", "WebSearch", "WebFetch"]), 0)[
        "valid"
    ]
    assert not run.validity("alt", _init(["Bash", "Read"]), 0)["valid"]


def test_the_dry_run_prints_the_alt_tools_and_prompts_and_calls_no_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    if not shutil.which("mafft") and not os.environ.get("MAFFT"):
        pytest.fail("MAFFT must be on PATH for the benchmark tests")
    monkeypatch.setattr(run, "BIN", tmp_path / "bin")
    real = run.subprocess.run

    def no_claude(cmd: Any, *a: Any, **k: Any) -> Any:
        if "claude" in str(cmd[0]):
            raise AssertionError(f"dry run called {cmd[0]}")
        return real(cmd, *a, **k)

    monkeypatch.setattr(run.subprocess, "run", no_claude)
    monkeypatch.setattr(
        sys, "argv",
        ["run.py", "--questions", str(FIELD / "questions.json"), "--arms", "alt",
         "--model", "claude-sonnet-5-5", "--reps", "1", "--dry-run"],
    )  # fmt: skip
    run.main()
    out = capsys.readouterr().out
    assert "== arm alt: tools must be exactly ['Bash', 'Read']" in out
    assert "--tools Bash,Read --allowedTools Bash,Read" in out
    assert "WebSearch" not in out and "WebFetch" not in out
    assert out.count("--- alt ") == 30
    assert run.ALT_NOTE.strip() in out  # the library note is unchanged
    assert "model claude-sonnet-5-5; 30 questions; reps 1" in out
    assert not (FIELD / "results.jsonl").exists()


# --- real execution ------------------------------------------------------------------

IQTREE = os.environ.get("IQTREE3") or shutil.which("iqtree3")


@pytest.mark.skipif(IQTREE is None, reason="needs IQ-TREE 3 (set IQTREE3)")
def test_alisim_reproduces_the_stored_alignments_from_the_generator() -> None:
    """Real execution of the data source: AliSim, given each P1 and P2
    question's generator (tree, model, length, seed), writes exactly the
    alignment the question ships. A different seed does not (the check can
    fail)."""
    for q in _questions("lba_sister")[:3] + _questions("supported_clades")[:2]:
        g = q["generator"]
        seqs = B.alisim(g["newick"], g["model"]["alisim"], g["sites"], g["alisim_seed"])
        assert B.sha256(T.fasta(seqs)) == g["alignment_sha256"], q["id"]
    q = _questions("lba_sister")[0]
    g = q["generator"]
    other = B.alisim(
        g["newick"], g["model"]["alisim"], g["sites"], g["alisim_seed"] + 1
    )
    assert B.sha256(T.fasta(other)) != g["alignment_sha256"]


@pytest.mark.skipif(IQTREE is None, reason="needs IQ-TREE 3 (set IQTREE3)")
def test_a_third_seed_of_the_p1_reference_pipeline_is_right_and_jc_is_not() -> None:
    """Real execution of the reference pipeline (IQ-TREE ModelFinder, BIC, then
    ML; never the server) at a seed the builder did not use, beside the JC
    shortcut on the same alignments."""
    for q in _questions("lba_sister"):
        seqs = B.read_fasta(q["files"]["alignment.fasta"])
        (asked,) = re.findall(r"cherry\) with (t\d+)\?", q["question"])
        ref = B.iqtree_ml(seqs, "MFP", 3)
        jc = B.iqtree_ml(seqs, "JC", 3)
        assert _v(q, str(B.cherry_partner(T.parse_newick(ref["newick"]), asked))) == (
            "right"
        ), q["id"]
        assert _v(q, str(B.cherry_partner(T.parse_newick(jc["newick"]), asked))) == (
            "wrong"
        ), q["id"]


def test_the_alt_arms_own_library_recovers_the_p1_cherry_under_gtr_g() -> None:
    """The route an alt run would take (piqtree, the repo's locked engine;
    not the server): GTR+F+G4 recovers every P1 cherry."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import piqtree
    from cogent3 import make_aligned_seqs

    for q in _questions("lba_sister"):
        seqs = B.read_fasta(q["files"]["alignment.fasta"])
        (asked,) = re.findall(r"cherry\) with (t\d+)\?", q["question"])
        tree = piqtree.build_tree(
            make_aligned_seqs(seqs, moltype="dna"), "GTR+F+G4", rand_seed=1
        )
        adj = T.parse_newick(tree.get_newick(with_distances=True))
        assert B.cherry_partner(adj, asked) == q["truth"]["taxon"], q["id"]


def test_the_field_builder_never_imports_the_server_or_its_libraries() -> None:
    """Read from the AST: the key's route is AliSim and the IQ-TREE/MAFFT
    CLIs, never phylokit_mcp (nor piqtree or cogent3, which it is built on)."""
    import ast

    for path in (FIELD / "build_field.py", BENCH / "score.py", BENCH / "controls.py"):
        mods: set[str] = set()
        for n in ast.walk(ast.parse(path.read_text())):
            if isinstance(n, ast.Import):
                mods.update(a.name for a in n.names)
            elif isinstance(n, ast.ImportFrom) and n.module:
                mods.add(n.module)
        banned = {
            m for m in mods if m.split(".")[0] in ("phylokit_mcp", "piqtree", "cogent3")
        }
        assert not banned, path.name
        assert "truth" in mods, path.name  # positive control: the walk sees imports


def test_the_builders_tolerance_step_reproduces_the_stored_bands() -> None:
    """``exact_tolerances`` run again on each stored question, with the bands
    it writes wiped, gives back the key's bands exactly: the code that made
    the key is under test, not only its output."""
    for q in json.loads(json.dumps(_questions())):
        if q["category"] == "lba_sister":
            continue
        want = dict(q["truth"])
        for k in ("nrf_max", "precision_min", "recall_min"):
            q["truth"].pop(k, None)
        q["source"] = q["source"].replace(
            "mean + max(3 sd, one split) over", "mean + 3 sd over"
        )
        B.exact_tolerances(q)
        for k in ("nrf_max", "precision_min", "recall_min"):
            if k in want:
                assert q["truth"][k] == want[k], (q["id"], k)
