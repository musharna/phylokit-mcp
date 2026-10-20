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


# --- server spec ---------------------------------------------------------------


def test_server_spec_takes_an_exact_release_or_a_full_commit_and_nothing_else() -> None:
    sha = "0123456789abcdef0123456789abcdef01234567"
    assert run.check_server_spec(run.SERVER) is None
    assert run.check_server_spec("phylokit-mcp==0.7.0") is None
    assert (
        run.check_server_spec(f"git+https://github.com/musharna/phylokit-mcp@{sha}")
        == sha
    )
    for bad in (
        "git+https://github.com/musharna/phylokit-mcp@main",  # a branch moves
        f"git+https://github.com/musharna/phylokit-mcp@{sha[:12]}",  # short sha
        f"git+https://github.com/someone/phylokit-mcp@{sha}",  # another repo
        "phylokit-mcp>=0.7.0",  # a range moves
        "phylokit-mcp==0.7",  # not one exact release
        "phylokit-mcp",
        "other-mcp==0.7.0",
    ):
        try:
            run.check_server_spec(bad)
        except SystemExit as exc:
            assert "--server-spec must be" in str(exc)
        else:
            raise AssertionError(f"{bad!r} was accepted")


def test_the_win_rule_needs_non_inferiority_and_a_significant_cost_cut() -> None:
    base = {
        "questions": 38, "a_right_every_run": 37, "b_right_every_run": 36,
        "questions_b_costlier": 5, "cost_sign_test_p": 1e-5, "cost_ratio_ci": (0.4, 0.7),
    }  # fmt: skip
    assert analyze.win_rule(base) == "WIN"
    # Two questions down is still non-inferior; three is not.
    assert analyze.win_rule(base | {"b_right_every_run": 35}) == "WIN"
    assert analyze.win_rule(base | {"b_right_every_run": 34}) == "LOSS"
    # Cheaper but not significantly, or a CI that reaches 1: a tie.
    assert analyze.win_rule(base | {"cost_sign_test_p": 0.06}) == "TIE"
    assert analyze.win_rule(base | {"cost_ratio_ci": (0.6, 1.02)}) == "TIE"
    # Significantly costlier is a loss.
    costly = base | {"questions_b_costlier": 33, "cost_ratio_ci": (1.5, 3.0)}
    assert analyze.win_rule(costly) == "LOSS"


def test_api_failures_are_told_apart_from_model_failures() -> None:
    init = {"type": "system", "subtype": "init"}
    done = {"type": "result", "subtype": "success", "terminal_reason": "completed"}
    allowed = {
        "type": "rate_limit_event",
        "rate_limit_info": {"status": "allowed_warning"},
    }
    rejected = {"type": "rate_limit_event", "rate_limit_info": {"status": "rejected"}}
    # A finished run, even a wrong or slow one, or one warned near its limit, is the model's.
    assert run.api_failure([init, allowed, done], 0) is None
    assert run.api_failure([init, done], "timeout") is None
    errored = {"type": "result", "is_error": True, "terminal_reason": "api_error"}
    assert run.api_failure([init, rejected, errored], 1) == "rate_limited"
    assert run.api_failure([init, errored], 1) == "api_error"
    assert run.api_failure([init], "timeout") == "timeout"


# --- model tiers (PREREG_2026-10-20_model-tiers) ---------------------------------


def test_a_non_default_model_writes_its_transcripts_to_its_own_dir() -> None:
    try:
        run.set_model("claude-haiku-5-5")
        assert run.MODEL == "claude-haiku-5-5"
        haiku = run.transcript("mcp", "q01", 1)
        assert haiku == run.ARTIFACTS / "runs-claude-haiku-5-5" / "mcp" / "q01-1.jsonl"
        run.set_model(run.DEFAULT_MODEL)
        # The default model's transcripts stay where the 10-16/10-18 ones are.
        assert run.transcript("mcp", "q01", 1) == run.ARTIFACTS / "runs/mcp/q01-1.jsonl"
    finally:
        run.set_model(run.DEFAULT_MODEL)


def _rl(status: str, five: float, seven: float) -> dict:
    return {
        "type": "rate_limit_event",
        "rate_limit_info": {
            "status": status, "resetsAt": 2000, "rateLimitType": "five_hour",
            "unifiedWindows": {
                "five_hour": {"utilization": five, "resetsAt": 1000},
                "seven_day": {"utilization": seven, "resetsAt": 9000},
            },
        },
    }  # fmt: skip


def test_the_usage_gate_pauses_at_the_five_hour_limit_and_stops_at_the_seven_day() -> (
    None
):
    # Below both limits, or no reading at all: carry on.
    assert run.usage_gate([_rl("allowed", 0.84, 0.89)])[0] is None
    assert run.usage_gate([{"type": "result"}]) == (None, {})
    # Five-hour at the limit: hold until that window resets.
    act, rec = run.usage_gate([_rl("allowed", 0.85, 0.5)])
    assert (act, rec["until"]) == ("pause", 1000)
    # A rejected request: hold until the reset the rejection names.
    act, rec = run.usage_gate([_rl("rejected", 0.5, 0.5)])
    assert (act, rec["until"]) == ("pause", 2000)
    # Seven-day at the limit stops, whatever the five-hour says.
    assert run.usage_gate([_rl("allowed", 0.9, 0.90)])[0] == "stop"
    # The last reading counts, not the first.
    assert (
        run.usage_gate([_rl("allowed", 0.9, 0.5), _rl("allowed", 0.1, 0.5)])[0] is None
    )
    # A pause with nothing to wait for is a stop, not a guess.
    bad = _rl("allowed", 0.9, 0.5)
    del bad["rate_limit_info"]["unifiedWindows"]["five_hour"]["resetsAt"]
    assert run.usage_gate([bad])[0] == "stop"


def test_the_gate_holds_later_runs_and_stops_on_another_model(
    tmp_path: Path, monkeypatch: Any
) -> None:
    calls: list[str] = []

    def fake_run_one(arm: str, q: dict, rep: int, spec: str) -> dict:
        calls.append(q["id"])
        path = run.transcript(arm, q["id"], rep)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(_rl("allowed", 0.1, 0.1)) + "\n")
        return {"model": served}

    monkeypatch.setattr(run, "RUNS", tmp_path)
    monkeypatch.setattr(run, "MODEL", "claude-haiku-5-5")
    monkeypatch.setattr(run, "run_one", fake_run_one)
    gate = run.Gate()
    served = "claude-haiku-5-5"
    # Served by the model asked for, under the limits: it runs, the gate stays open.
    assert run.gated_run(gate, "mcp", {"id": "q01"}, 1, run.SERVER) == {"model": served}
    assert gate.stopped is None and gate.wait() is True and gate.until == 0.0
    # A pause reading moves the hold to the reset time (+60 s), and is logged.
    gate.update([_rl("allowed", 0.9, 0.1)], served)
    assert gate.until == 1060 and gate.log[-1]["action"] == "pause"
    gate.until = 0.0
    # Served by another model: the gate stops, and no later run starts.
    served = "claude-opus-5-5"
    run.gated_run(gate, "mcp", {"id": "q02"}, 1, run.SERVER)
    assert (
        gate.stopped is not None and "served by claude-opus-5-5" in gate.stopped["why"]
    )
    assert run.gated_run(gate, "mcp", {"id": "q03"}, 1, run.SERVER) is None
    assert calls == ["q01", "q02"]


def test_holm_matches_hand_values() -> None:
    assert analyze.holm({"a": 0.01, "b": 0.04}) == {"a": 0.02, "b": 0.04}
    # Step-down is monotone: the larger p never drops below the smaller's adjustment.
    assert analyze.holm({"a": 0.03, "b": 0.02}) == {"b": 0.04, "a": 0.04}
    assert analyze.holm({"a": 0.7, "b": 0.6}) == {"b": 1.0, "a": 1.0}


def _tier(a_only: int, b_only: int, p: float, ci: tuple[float, float]) -> dict:
    return {
        "a_only": a_only, "b_only": b_only, "sign_test_p": p,
        "right_every_run_diff_ci": ci,
    }  # fmt: skip


def test_the_tier_rule_needs_a_holm_significant_lift_and_no_drop() -> None:
    lift = _tier(0, 9, 0.004, (4.0, 13.0))
    drop = _tier(9, 0, 0.004, (-13.0, -4.0))
    flat = _tier(2, 3, 1.0, (-3.0, 4.0))
    assert analyze.tier_rule({"haiku": lift, "sonnet": flat}) == (
        "WIN", {"haiku": "LIFT", "sonnet": "NONE"},
    )  # fmt: skip
    assert analyze.tier_rule({"haiku": drop, "sonnet": flat})[0] == "LOSS"
    assert analyze.tier_rule({"haiku": lift, "sonnet": drop})[0] == "MIXED"
    assert analyze.tier_rule({"haiku": flat, "sonnet": flat})[0] == "TIE"
    # p = 0.03 is significant raw but not after Holm over two tiers.
    near = _tier(0, 6, 0.03, (2.0, 9.0))
    assert analyze.tier_rule({"haiku": near, "sonnet": flat})[1]["haiku"] == "NONE"
    # Significant, but a CI that reaches 0 is not a lift.
    assert analyze.tier_rule({"haiku": lift | {"right_every_run_diff_ci": (0.0, 9.0)},
                              "sonnet": flat})[1]["haiku"] == "NONE"  # fmt: skip
    # One tier in: pending, and judged with the full family's Holm factor.
    assert analyze.tier_rule({"haiku": lift}) == ("PENDING", {"haiku": "LIFT"})
    assert analyze.tier_rule({"haiku": near})[1]["haiku"] == "NONE"


def test_paired_rows_report_the_accuracy_difference_and_the_server_call_count() -> None:
    rows = []
    for i in range(10):
        for rep in (1, 2, 3):
            rows.append(_row("alt", f"q{i}", rep, "right" if i < 4 else "wrong"))
            m = _row("mcp", f"q{i}", rep, "right")
            m["tool_calls"] = (
                ["Read"] if i == 0 else ["Read", "mcp__phylokit__infer_tree"]
            )
            rows.append(m)
    (p,) = analyze.paired_rows(rows, np.random.default_rng(3))
    assert p["right_every_run_diff_b_minus_a"] == 6
    lo, hi = p["right_every_run_diff_ci"]
    assert 0 < lo <= 6 <= hi <= 10
    arms = {
        r["arm"]: r
        for r in analyze.arm_rows(rows, np.random.default_rng(3))
        if r["level"] == "arm"
    }
    assert arms["mcp"]["runs_without_server_call"] == 3
    assert abs(arms["mcp"]["server_calls_mean"] - 0.9) < 1e-9
    assert "runs_without_server_call" not in arms["alt"]


def test_a_release_pin_must_install_that_release(monkeypatch: Any) -> None:
    """check_mcp reads the installed version back; a resolver that hands over
    another release (here 0.6.0 for ==0.7.0) is refused before any run."""
    checked = {"piqtree": "0.8.3"}

    def fake(installed: str) -> Any:
        def sub(cmd: list[str], **kw: Any) -> Any:
            code = cmd[-1]
            out = (
                {"piqtree": "0.8.3", "cogent3": "x", "numpy": "y"}
                if "piqtree" in code
                else {"phylokit-mcp": installed, "commit": None}
            )
            return type("P", (), {"returncode": 0, "stdout": json.dumps(out)})()

        return sub

    monkeypatch.setattr(run, "_uvx_server_args", lambda spec: ["--from", spec, "x"])
    monkeypatch.setattr(run.shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setattr(run.subprocess, "run", fake("0.7.0"))
    assert run.check_mcp(checked, "phylokit-mcp==0.7.0")["phylokit-mcp"] == "0.7.0"
    monkeypatch.setattr(run.subprocess, "run", fake("0.6.0"))
    try:
        run.check_mcp(checked, "phylokit-mcp==0.7.0")
    except SystemExit as exc:
        assert "installed '0.6.0'" in str(exc)
    else:
        raise AssertionError("0.6.0 was accepted for ==0.7.0")


def test_tiers_main_prints_the_holm_p_the_verdict_was_made_on(
    tmp_path: Path, capsys: Any
) -> None:
    rows = []
    for i in range(12):
        for rep in (1, 2, 3):
            rows.append(_row("alt", f"q{i}", rep, "right" if i < 4 else "wrong"))
            rows.append(_row("mcp", f"q{i}", rep, "right"))
    path = tmp_path / "results.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    # 8 discordant questions, all for mcp: raw p = 2/256 = 0.0078; one tier of
    # two planned, so Holm's factor is 2 (0.0156), and the call is a lift.
    analyze.tiers_main(f"haiku={path}")
    out = capsys.readouterr().out
    assert "raw p=0.00781, Holm p=0.0156" in out and ": LIFT" in out
    assert "tier rule (PREREG_2026-10-20_model-tiers): PENDING" in out
    # A full family of two equal tiers: Holm p = 2 * p for both, and a win.
    analyze.tiers_main(f"haiku={path},sonnet={path}")
    out = capsys.readouterr().out
    assert out.count("Holm p=0.0156") == 2 and "): WIN" in out
