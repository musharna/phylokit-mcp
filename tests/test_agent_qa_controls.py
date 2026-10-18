"""Controls for the agent benchmark: the oracle (positive control), the no-LLM
baselines (negative controls), the validity rule and leak scan applied to
every arm, and the statistics in analyze.py.

Each check that something is refused sits beside a check that the legitimate
case passes, in the same test, so a harness that refused everything (or
nothing) fails here.
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

BENCH = Path(__file__).resolve().parent.parent / "benchmarks" / "agent_qa"
sys.path.insert(0, str(BENCH))

import analyze
import controls
import run
from score import score


def _questions() -> list[dict[str, Any]]:
    return json.loads((BENCH / "questions.json").read_text())["questions"]


# --- oracle --------------------------------------------------------------------


def test_the_key_written_as_replies_scores_all_right_and_wrong_replies_all_wrong() -> (
    None
):
    qs = _questions()
    right = [score(q, controls.reply(controls.key_value(q)))["verdict"] for q in qs]
    wrong = [score(q, controls.reply(controls.wrong_value(q)))["verdict"] for q in qs]
    declined = [score(q, "ANSWER: UNKNOWN")["verdict"] for q in qs]
    assert right == ["right"] * 38
    assert wrong == ["wrong"] * 38
    assert declined == ["declined"] * 38
    assert controls.oracle_check(qs) == []


def test_the_oracle_check_names_a_question_whose_key_disagrees_with_the_scorer() -> (
    None
):
    qs = _questions()
    broken = json.loads(json.dumps(qs))
    broken[0]["truth"]["taxa"] = [t for t in broken[0]["truth"]["taxa"] if t != "t21"]
    assert broken[0]["truth"]["taxon"] == "t21"  # the edit removed the key's own answer
    bad = controls.oracle_check(broken)
    assert len(bad) == 1 and bad[0].startswith("q01: key answer 't21' scored wrong")
    assert controls.oracle_check(qs) == []


# --- baselines ---------------------------------------------------------------


def test_the_baselines_on_the_frozen_key() -> None:
    """The numbers the README quotes. Same residue number and variable sites
    get none by construction (the asked residue is >= 5 positions off, and
    every alignment has constant-but-two-state columns). The flanking motif
    solves 5 of 6 homology questions with no alignment at all: the no-gap
    margin that makes the key certain also makes the residue findable by
    context. RF on the groups as written gets none, since each string is
    rooted at a random internal node."""
    rows = controls.baselines(_questions())
    right: dict[str, int] = {}
    total: dict[str, int] = {}
    for r in rows:
        total[r["baseline"]] = total.get(r["baseline"], 0) + 1
        right[r["baseline"]] = right.get(r["baseline"], 0) + (r["verdict"] == "right")
    assert total["best constant"] == 38
    assert right == {
        "best constant": 10,
        "nearest taxon by JC distance": 4,
        "neighbour-joining cherry partner": 5,
        "same residue number": 0,
        "same relative position": 0,
        "flanking-motif match (+-6, no alignment)": 5,
        "RF on the bracketed groups as written": 0,
        "uncorrected p-distance": 3,
        "variable sites": 0,
        "composition test + transition share": 5,
    }


def test_the_controls_never_import_the_server_or_its_libraries() -> None:
    for name in ("controls.py", "analyze.py"):
        tree = ast.parse((BENCH / name).read_text())
        mods: set[str] = set()
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                mods.update(a.name for a in n.names)
            elif isinstance(n, ast.ImportFrom) and n.module:
                mods.add(n.module)
        assert not {
            m for m in mods if m.split(".")[0] in ("phylokit_mcp", "cogent3", "piqtree")
        }
        assert "numpy" in mods, name  # positive control: the walk sees real imports


# --- validity and leaks --------------------------------------------------------

OWN = "/x/tmp/phyloqa-alt-aaaa1111/cwd"


def _events(
    tools: list[str],
    calls: list[tuple[str, dict]] = (),  # type: ignore[assignment]
    outputs: list[str] = (),  # type: ignore[assignment]
    servers: list[dict] | None = None,
    subtype: str = "success",
) -> list[dict[str, Any]]:
    ev: list[dict[str, Any]] = [
        {"type": "system", "subtype": "init", "cwd": OWN, "tools": tools,
         "mcp_servers": servers or []}
    ]  # fmt: skip
    for name, inp in calls:
        ev.append({"type": "assistant", "message": {"content": [
            {"type": "tool_use", "name": name, "input": inp}]}})  # fmt: skip
    for text in outputs:
        ev.append({"type": "user", "message": {"content": [
            {"type": "tool_result", "content": text}]}})  # fmt: skip
    ev.append({"type": "result", "subtype": subtype, "is_error": subtype != "success"})
    return ev


ALT = ["Bash", "Read", "WebSearch", "WebFetch"]
READ_OWN = ("Read", {"file_path": OWN + "/alignment.fasta"})


def test_a_clean_run_of_each_arm_is_valid() -> None:
    assert run.validity("bare", _events([]), 0)["valid"]
    assert run.validity("read", _events(["Read"], [READ_OWN]), 0)["valid"]
    assert run.validity("alt", _events(ALT, [READ_OWN], ["ok " + OWN]), 0)["valid"]
    mcp = _events(
        ["Read", "mcp__phylokit__infer_tree"],
        [READ_OWN, ("mcp__phylokit__infer_tree", {"alignment_fasta": ">t1\nACGT"})],
        servers=[{"name": "phylokit", "status": "connected"}],
    )
    assert run.validity("mcp", mcp, 0)["valid"]


def test_the_same_rule_marks_failures_invalid_for_every_arm() -> None:
    v = run.validity
    assert v("alt", _events(ALT, [READ_OWN]), "timeout")["invalid_reasons"] == [
        "exit timeout"
    ]
    assert (
        "result error_max_turns"
        in v("read", _events(["Read"], [READ_OWN], subtype="error_max_turns"), 0)[
            "invalid_reasons"
        ]
    )
    assert v("bare", _events(["Read"]), 0)["invalid_reasons"] == [
        "tools ['Read'] != []"
    ]
    assert v("alt", _events(ALT), 0)["invalid_reasons"] == ["no tool call"]
    assert v("read", _events(["Read", "Bash"], [READ_OWN]), 0)["valid"] is False
    failed = _events(
        ["Read", "mcp__phylokit__infer_tree"], [READ_OWN],
        servers=[{"name": "phylokit", "status": "failed"}],
    )  # fmt: skip
    assert v("mcp", failed, 0)["invalid_reasons"] == ["server failed"]
    assert v("mcp", _events(["Read"], [READ_OWN]), 0)["invalid_reasons"] == [
        "no phylokit tools", "server None"
    ]  # fmt: skip
    assert v(
        "alt", _events(ALT, [("Bash", {"command": "pip install phylokit-mcp"})]), 0
    )["invalid_reasons"] == ["named phylokit"]
    # Positive control: the same arm, legitimately run, is valid.
    assert v("alt", _events(ALT, [("Bash", {"command": "python3 -c 1"})]), 0)["valid"]


def test_a_leak_in_a_tool_input_or_output_invalidates_but_the_runs_own_dir_does_not() -> (
    None
):
    v = run.validity
    clean = v("alt", _events(ALT, [("Bash", {"command": f"ls {OWN}/.."})], [OWN]), 0)
    assert clean["valid"] and clean["leaks"] == []
    by_input = v(
        "alt", _events(ALT, [("Bash", {"command": "cat ../../questions.json"})]), 0
    )
    assert by_input["leaks"] == ["questions.json"] and not by_input["valid"]
    by_output = v("alt", _events(ALT, [("Bash", {"command": "ls /x/tmp"})],
                                 ["phyloqa-alt-aaaa1111\nphyloqa-mcp-bbbb2222"]), 0)  # fmt: skip
    assert by_output["leaks"] == ["phyloqa-mcp-bbbb2222"]
    repo = "https://raw.githubusercontent.com/musharna/phylokit-mcp/main/benchmarks/agent_qa/x"
    web = v("alt", _events(ALT, [("WebFetch", {"url": repo})]), 0)
    assert not web["valid"] and "agent_qa" in web["leaks"]
    trace = (
        "/home/u/phylokit-mcp/.claude/worktrees/b/benchmarks/agent_qa/.alt-env/lib/x.py"
    )
    tb = v("alt", _events(ALT, [("Bash", {"command": "python3 x.py"})], [trace]), 0)
    assert tb["leaks"] == [".claude/worktrees", "agent_qa"]
    # A web search that only names the server is not a leak (the arm rule
    # handles naming it), and a plain GitHub page of the server is not either.
    home = v(
        "alt",
        _events(
            ALT, [("WebFetch", {"url": "https://github.com/musharna/phylokit-mcp"})]
        ),
        0,
    )
    assert home["leaks"] == []


def test_a_run_directory_that_names_the_benchmark_is_refused_before_the_run() -> None:
    """The 2026-10-17 read-arm restart put run directories under
    ~/phylokit-mcp-artifacts/agent_qa/tmp; every Read path then named both,
    and the first two runs were invalid. The check refuses such a path."""
    import pytest

    for bad in ("/home/u/phylokit-mcp-artifacts/agent_qa/tmp/phyloqa-read-x",
                "/home/u/phylokit/tmp/phyloqa-read-x"):  # fmt: skip
        with pytest.raises(SystemExit, match="names"):
            run.check_run_dir(bad)
    run.check_run_dir("/tmp/phyloqa-read-abc123")  # positive control


# --- statistics ------------------------------------------------------------------


def test_the_exact_sign_test_matches_hand_values() -> None:
    assert analyze.sign_test(0, 2) == 0.5
    assert analyze.sign_test(0, 6) == 2 / 64
    assert analyze.sign_test(6, 0) == 2 / 64
    assert analyze.sign_test(3, 3) == 1.0
    assert analyze.sign_test(0, 0) == 1.0
    assert abs(analyze.sign_test(1, 9) - 2 * 11 / 1024) < 1e-12


def test_bootstrap_intervals_cover_the_estimate_and_collapse_when_every_question_agrees() -> (
    None
):
    rng = np.random.default_rng(1)
    qids = [f"q{i}" for i in range(20)]
    right = {q: i < 15 for i, q in enumerate(qids)}

    def share(qs: list[str]) -> float:
        return sum(right[q] for q in qs) / len(qs)

    lo, hi = analyze.bootstrap(qids, share, rng, b=4000)
    assert lo < 0.75 < hi and 0.5 < lo and hi <= 0.95
    assert analyze.bootstrap(qids, lambda qs: 1.0, rng, b=200) == (1.0, 1.0)


def _row(
    arm: str, qid: str, rep: int, verdict: str, valid: bool = True, s: float = 10.0
) -> dict:
    return {
        "arm": arm, "qid": qid, "category": "rf_distance", "rep": rep, "verdict": verdict,
        "valid": valid, "wrong_id": False, "cost_usd": s / 100, "seconds": s,
        "tool_calls": ["Bash"],
    }  # fmt: skip


def test_an_invalid_row_counts_as_not_right_and_pairs_are_compared_per_question() -> (
    None
):
    rows = []
    for i in range(10):
        for rep in (1, 2):
            rows.append(_row("alt", f"q{i}", rep, "right", s=10))
            # mcp misses q0 once, and q1's second run is a right answer that is invalid.
            verdict = "wrong" if (i == 0 and rep == 1) else "right"
            rows.append(
                _row(
                    "mcp", f"q{i}", rep, verdict, valid=not (i == 1 and rep == 2), s=40
                )
            )
    rng = np.random.default_rng(2)
    arms = {r["arm"]: r for r in analyze.arm_rows(rows, rng) if r["level"] == "arm"}
    assert (arms["mcp"]["right"], arms["mcp"]["invalid"], arms["mcp"]["n"]) == (
        18,
        1,
        20,
    )
    assert arms["mcp"]["questions_right_every_run"] == 8
    assert arms["alt"]["questions_right_every_run"] == 10
    (p,) = analyze.paired_rows(rows, rng)
    assert (p["a"], p["b"], p["a_only"], p["b_only"], p["both"]) == (
        "alt",
        "mcp",
        2,
        0,
        8,
    )
    assert p["sign_test_p"] == 0.5
    assert abs(p["seconds_ratio_b_over_a"] - 4.0) < 1e-9
    assert p["questions_b_slower"] == 10 and p["seconds_sign_test_p"] == 2 / 1024


def test_the_needs_compute_subset_keeps_only_its_kinds_and_pairs_on_them() -> None:
    """The post-hoc subset (sister, model, informative sites) is scored on its
    own questions: here read misses both rf questions (outside the subset) and
    one of the two sister questions, so the subset row shows 1 of 2 questions,
    not 1 of 4, and the paired subset row compares the 2 sister questions only."""
    rows = []
    for i, cat in enumerate(
        ["rf_distance", "rf_distance", "sister_taxon", "sister_taxon"]
    ):
        for rep in (1, 2):
            rows.append({**_row("alt", f"q{i}", rep, "right"), "category": cat})
            verdict = "right" if i == 2 else "wrong"
            rows.append({**_row("read", f"q{i}", rep, verdict), "category": cat})
    rng = np.random.default_rng(3)
    sub = {
        r["arm"]: r for r in analyze.arm_rows(rows, rng) if r["level"] == "arm_x_subset"
    }
    assert (sub["read"]["right"], sub["read"]["n"], sub["read"]["questions"]) == (
        2,
        4,
        2,
    )
    assert sub["read"]["questions_right_every_run"] == 1
    assert sub["read"]["category"] == "needs_compute"
    assert sub["alt"]["accuracy"] == 1.0
    by = {p["subset"]: p for p in analyze.paired_rows(rows, rng)}
    assert (by["all"]["questions"], by["all"]["b_only"]) == (4, 3)
    assert (by["needs_compute"]["questions"], by["needs_compute"]["b_only"]) == (2, 1)
    # Results with none of the subset's kinds get no subset row (no division by 0).
    rf_only = [r for r in rows if r["category"] == "rf_distance"]
    assert not [
        r for r in analyze.arm_rows(rf_only, rng) if r["level"] == "arm_x_subset"
    ]
