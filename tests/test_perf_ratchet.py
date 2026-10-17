"""The work-count ratchet (benchmarks/perf/ratchet.py): every tool does exactly
the work recorded in benchmarks/perf/baseline.json on its fixed input.

Counts, not seconds, so the gate has no headroom: one extra call fails it.
Measured identical on Python 3.12, 3.13 and 3.14 and under different hash
seeds before the baseline was committed.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

PERF = Path(__file__).resolve().parent.parent / "benchmarks" / "perf"
sys.path.insert(0, str(PERF))

import ratchet

# mutmut 3 routes every function through a trampoline, which is itself a call
# under the package root, so under mutmut every count moves and the gate would
# read the trampolines, not the mutant. mutmut sets this variable for its
# clean run and for every mutant.
pytestmark = pytest.mark.skipif(
    "MUTANT_UNDER_TEST" in os.environ,
    reason="mutmut trampolines add calls under the package root",
)


def test_every_tool_does_exactly_the_recorded_work() -> None:
    baseline = json.loads(ratchet.BASELINE.read_text())
    assert set(baseline) == set(ratchet.WORKLOADS)
    current = ratchet.measure()
    problems = ratchet.compare(baseline, current)
    assert not problems, (
        "work counts differ from benchmarks/perf/baseline.json:\n"
        + "\n".join(problems)
        + "\n\n"
        + ratchet.summary(current)
    )


def _row(calls: dict[str, int], witness: object = 1) -> dict:
    return {
        "witness": witness,
        "python_calls": calls,
        "engine_calls": {"piqtree:build_tree": 3},
        "subprocesses": {},
    }


def test_a_rise_names_the_function_and_its_delta_and_an_untouched_row_passes() -> None:
    base = {"a": _row({"x.py:f": 5, "x.py:g": 2}), "b": _row({"y.py:h": 1})}
    more = {"a": _row({"x.py:f": 6, "x.py:g": 2}), "b": _row({"y.py:h": 1})}
    problems = ratchet.compare(base, more)
    assert len(problems) == 1
    assert problems[0].startswith("a: python_calls ROSE")
    assert "7 -> 8 (+1)" in problems[0]
    assert "+1  x.py:f  (5 -> 6)" in problems[0]
    assert "x.py:g" not in problems[0]  # unchanged keys are not listed
    assert ratchet.compare(base, base) == []  # positive control: equal is green


def test_a_drop_also_fails_until_it_is_locked_in() -> None:
    base = {"a": _row({"x.py:f": 5})}
    fewer = {"a": _row({"x.py:f": 4})}
    (problem,) = ratchet.compare(base, fewer)
    assert "FELL" in problem and "--update" in problem and "(-1)" in problem
    assert ratchet.compare(fewer, fewer) == []


def test_a_workload_whose_output_changed_fails_even_with_equal_counts() -> None:
    """The harness assertion: the witness is an output only that workload
    produces, so a run of the wrong workload cannot pass by its name."""
    base = {"a": _row({"x.py:f": 5}, witness={"rf": 2})}
    (problem,) = ratchet.compare(base, {"a": _row({"x.py:f": 5}, witness={"rf": 4})})
    assert "witness" in problem
    assert ratchet.compare(base, base) == []


def _update(
    monkeypatch, tmp_path: Path, base: dict, current: dict, *flags: str
) -> tuple[int, dict]:
    """ratchet.main with the measurement and the baseline file swapped for
    these rows; returns the exit code and the baseline file afterwards."""
    path = tmp_path / "baseline.json"
    path.write_text(json.dumps(base))
    monkeypatch.setattr(ratchet, "BASELINE", path)
    monkeypatch.setattr(ratchet, "WORKLOADS", tuple(current))
    monkeypatch.setattr(ratchet, "measure", lambda names=None: current)
    code = ratchet.main(["--update", *flags])
    return code, json.loads(path.read_text())


def test_update_refuses_a_rise_without_allow_increase_and_takes_a_drop(
    monkeypatch, tmp_path: Path
) -> None:
    base = {"a": _row({"x.py:f": 5})}
    more = {"a": _row({"x.py:f": 6})}
    code, after = _update(monkeypatch, tmp_path, base, more)
    assert code == 1 and after == base  # refused, file untouched
    code, after = _update(monkeypatch, tmp_path, base, more, "--allow-increase")
    assert code == 0 and after == more  # an explicit rise is written
    fewer = {"a": _row({"x.py:f": 4})}
    code, after = _update(monkeypatch, tmp_path, base, fewer)
    assert code == 0 and after == fewer  # positive control: a drop locks in


def test_update_refuses_a_changed_output_even_when_the_counts_are_equal(
    monkeypatch, tmp_path: Path
) -> None:
    base = {"a": _row({"x.py:f": 5}, witness={"sha256": "aaa"})}
    other = {"a": _row({"x.py:f": 5}, witness={"sha256": "bbb"})}
    code, after = _update(monkeypatch, tmp_path, base, other)
    assert code == 1 and after == base
    code, after = _update(monkeypatch, tmp_path, base, base)
    assert code == 0 and after == base  # positive control: same output passes


def test_an_unknown_workload_is_refused_not_run_as_a_default() -> None:
    bad = subprocess.run(
        [sys.executable, str(PERF / "ratchet.py"), "--child", "infer_tre"],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert bad.returncode != 0
    assert "unknown workload 'infer_tre'" in bad.stderr
    good = ratchet.measure(("compare_trees",))
    # The two fixture trees differ by one split each side, worked by hand:
    # {a,b}|rest against {a,c}|rest.
    assert good["compare_trees"]["witness"] == {"rf": 2, "only_in_a": ["c|d|e|f|g|h"]}
