"""Tables and statistics for the agent benchmark, from results.jsonl.

    uv run python benchmarks/agent_qa/analyze.py [results.jsonl] [outdir]

Writes, next to the results:

- ``summary.tsv``: per arm (and arm x kind of question) the verdict counts,
  invalid rows, accuracy, questions right in every run, cost, time and tool
  calls; per baseline (``controls.py``, no LLM) the same verdict columns.
  Arm-level and baseline-level rows carry bootstrap 95% intervals over
  QUESTIONS (resampled with replacement, all runs of a drawn question kept
  together, since the runs of one question are not independent).
- ``paired.tsv``: each pair of arms on the same questions: how many questions
  each got right in every run, the discordant counts and an exact two-sided
  sign test (McNemar's exact test) on them, and the per-question ratio of
  median seconds and median cost with a bootstrap 95% interval.
- ``baselines.tsv``: every baseline's answer and verdict, per question.

An invalid row (``run.py``'s ``validity``) counts as not right: one rule for
every arm. ``accuracy_valid_only`` is the same share over valid rows alone.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from controls import baselines

HERE = Path(__file__).resolve().parent
ARM_ORDER = ["bare", "read", "alt", "mcp"]
CATEGORIES = [
    "sister_taxon",
    "rf_distance",
    "homologous_position",
    "model_selection",
    "pairwise_distance",
    "informative_sites",
]
B = 10_000
SEED = 20261016
PAIRS = [("alt", "mcp"), ("read", "alt"), ("read", "mcp")]


def sign_test(b: int, c: int) -> float:
    """Exact two-sided sign test on b vs c discordant pairs (McNemar exact)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2**n
    return min(1.0, 2 * tail)


def bootstrap(
    qids: list[str],
    stat: Any,
    rng: np.random.Generator,
    b: int = B,
) -> tuple[float, float]:
    """2.5 and 97.5 percentiles of ``stat(drawn_qids)`` over ``b`` resamples
    of the questions; NaN draws (a statistic undefined on that draw) dropped."""
    idx = rng.integers(0, len(qids), size=(b, len(qids)))
    vals = np.array([stat([qids[i] for i in row]) for row in idx], dtype=float)
    vals = vals[~np.isnan(vals)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def _ok(r: dict[str, Any]) -> bool:
    return r["verdict"] == "right" and r.get("valid", True)


def _fmt(x: float | None, nd: int = 4) -> str:
    return (
        "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}g}"
    )


def arm_rows(
    rows: list[dict[str, Any]], rng: np.random.Generator
) -> list[dict[str, Any]]:
    out = []
    arms = [a for a in ARM_ORDER if any(r["arm"] == a for r in rows)]
    for arm in arms:
        mine = [r for r in rows if r["arm"] == arm]
        by_q: dict[str, list[dict[str, Any]]] = {}
        for r in mine:
            by_q.setdefault(r["qid"], []).append(r)
        qids = sorted(by_q)

        def acc(qs: list[str], by_q: dict = by_q) -> float:
            rs = [r for q in qs for r in by_q[q]]
            return sum(map(_ok, rs)) / len(rs)

        def all_right(qs: list[str], by_q: dict = by_q) -> float:
            return sum(all(map(_ok, by_q[q])) for q in qs) / len(qs)

        def med(key: str, by_q: dict = by_q) -> Any:
            return lambda qs: float(
                np.median([r[key] for q in qs for r in by_q[q] if r[key] is not None])
            )

        valid = [r for r in mine if r.get("valid", True)]
        row = _tally("arm", arm, "all", mine)
        row.update(
            accuracy=acc(qids),
            accuracy_ci=bootstrap(qids, acc, rng),
            accuracy_valid_only=(sum(map(_ok, valid)) / len(valid))
            if valid
            else float("nan"),
            questions=len(qids),
            questions_right_every_run=sum(all(map(_ok, by_q[q])) for q in qids),
            questions_right_every_run_ci=tuple(
                len(qids) * x for x in bootstrap(qids, all_right, rng)
            ),
            cost_usd_median=med("cost_usd")(qids),
            cost_usd_median_ci=bootstrap(qids, med("cost_usd"), rng),
            cost_usd_total=sum(r["cost_usd"] or 0 for r in mine),
            seconds_median=med("seconds")(qids),
            seconds_median_ci=bootstrap(qids, med("seconds"), rng),
            seconds_max=max(r["seconds"] for r in mine),
            tool_calls_mean=float(np.mean([len(r["tool_calls"]) for r in mine])),
        )
        out.append(row)
        for cat in CATEGORIES:
            sub = [r for r in mine if r["category"] == cat]
            if sub:
                c = _tally("arm_x_category", arm, cat, sub)
                c["questions"] = len({r["qid"] for r in sub})
                c["questions_right_every_run"] = sum(
                    all(map(_ok, by_q[q])) for q in {r["qid"] for r in sub}
                )
                out.append(c)
    return out


def _tally(level: str, name: str, cat: str, rs: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "level": level,
        "arm": name,
        "category": cat,
        "n": len(rs),
        "right": sum(map(_ok, rs)),
        "wrong": sum(r["verdict"] == "wrong" and r.get("valid", True) for r in rs),
        "declined": sum(
            r["verdict"] == "declined" and r.get("valid", True) for r in rs
        ),
        "no_answer": sum(
            r["verdict"] == "no_answer" and r.get("valid", True) for r in rs
        ),
        "invalid": sum(not r.get("valid", True) for r in rs),
        "wrong_id": sum(bool(r.get("wrong_id")) for r in rs),
    }


def baseline_rows(
    base: list[dict[str, Any]], rng: np.random.Generator
) -> list[dict[str, Any]]:
    out = []
    names = list(dict.fromkeys(b["baseline"] for b in base))
    for name in names:
        mine = [b for b in base if b["baseline"] == name]
        by_q = {b["qid"]: b for b in mine}
        qids = sorted(by_q)

        def acc(qs: list[str], by_q: dict = by_q) -> float:
            return sum(by_q[q]["verdict"] == "right" for q in qs) / len(qs)

        cats = sorted({b["category"] for b in mine}, key=CATEGORIES.index)
        row = _tally("baseline", name, "all" if len(cats) > 1 else cats[0], mine)
        row.update(
            accuracy=acc(qids),
            accuracy_ci=bootstrap(qids, acc, rng),
            questions=len(qids),
            questions_right_every_run=row["right"],
        )
        out.append(row)
        if len(cats) > 1:
            for cat in cats:
                sub = [b for b in mine if b["category"] == cat]
                c = _tally("baseline_x_category", name, cat, sub)
                c["questions"] = c["questions_right_every_run"] = None
                c["questions"] = len(sub)
                c["questions_right_every_run"] = c["right"]
                out.append(c)
    return out


def paired_rows(
    rows: list[dict[str, Any]], rng: np.random.Generator
) -> list[dict[str, Any]]:
    out = []
    arms = {r["arm"] for r in rows}
    for a, b in PAIRS:
        if a not in arms or b not in arms:
            continue
        qa: dict[str, list[dict[str, Any]]] = {}
        qb: dict[str, list[dict[str, Any]]] = {}
        for r in rows:
            if r["arm"] == a:
                qa.setdefault(r["qid"], []).append(r)
            elif r["arm"] == b:
                qb.setdefault(r["qid"], []).append(r)
        qids = sorted(set(qa) & set(qb))
        ra = {q: all(map(_ok, qa[q])) for q in qids}
        rb = {q: all(map(_ok, qb[q])) for q in qids}
        a_only = sum(ra[q] and not rb[q] for q in qids)
        b_only = sum(rb[q] and not ra[q] for q in qids)

        def ratio(key: str, qa: dict = qa, qb: dict = qb) -> Any:
            def f(qs: list[str]) -> float:
                logs = [
                    math.log(
                        np.median([r[key] for r in qb[q]])
                        / np.median([r[key] for r in qa[q]])
                    )
                    for q in qs
                ]
                return math.exp(float(np.mean(logs)))

            return f

        slower = sum(
            np.median([r["seconds"] for r in qb[q]])
            > np.median([r["seconds"] for r in qa[q]])
            for q in qids
        )
        out.append(
            {
                "a": a,
                "b": b,
                "questions": len(qids),
                "a_right_every_run": sum(ra.values()),
                "b_right_every_run": sum(rb.values()),
                "both": sum(ra[q] and rb[q] for q in qids),
                "a_only": a_only,
                "b_only": b_only,
                "neither": sum(not ra[q] and not rb[q] for q in qids),
                "sign_test_p": sign_test(a_only, b_only),
                "seconds_ratio_b_over_a": ratio("seconds")(qids),
                "seconds_ratio_ci": bootstrap(qids, ratio("seconds"), rng),
                "questions_b_slower": slower,
                "seconds_sign_test_p": sign_test(slower, len(qids) - slower),
                "cost_ratio_b_over_a": ratio("cost_usd")(qids),
                "cost_ratio_ci": bootstrap(qids, ratio("cost_usd"), rng),
            }
        )
    return out


def _flatten(row: dict[str, Any]) -> dict[str, str]:
    flat: dict[str, str] = {}
    for k, v in row.items():
        if isinstance(v, tuple):
            flat[f"{k}_lo"], flat[f"{k}_hi"] = _fmt(v[0]), _fmt(v[1])
        elif isinstance(v, float):
            flat[k] = _fmt(v)
        else:
            flat[k] = "NA" if v is None else str(v)
    return flat


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    flat = [_flatten(r) for r in rows]
    cols = list(dict.fromkeys(k for r in flat for k in r))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="NA", lineterminator="\n")
        w.writeheader()
        w.writerows(flat)


def main(argv: list[str]) -> None:
    results = Path(argv[0]) if argv else HERE / "results.jsonl"
    outdir = Path(argv[1]) if len(argv) > 1 else results.parent
    rows = [json.loads(ln) for ln in results.read_text().splitlines()]
    questions = json.loads((HERE / "questions.json").read_text())["questions"]
    rng = np.random.default_rng(SEED)
    base = baselines(questions)
    summary = arm_rows(rows, rng) + baseline_rows(base, rng)
    write_tsv(outdir / "summary.tsv", summary)
    paired = paired_rows(rows, rng)
    write_tsv(outdir / "paired.tsv", paired)
    write_tsv(outdir / "baselines.tsv", base)
    for r in summary:
        if r["level"] in ("arm", "baseline"):
            lo, hi = r["accuracy_ci"]
            print(
                f"{r['level']:8} {r['arm']:38} {r['right']:>3}/{r['n']:<3} "
                f"{r['accuracy']:.3f} [{lo:.3f}, {hi:.3f}] invalid={r['invalid']}"
            )
    for p in paired:
        print(
            f"{p['b']} vs {p['a']}: right every run {p['b_right_every_run']} vs "
            f"{p['a_right_every_run']} of {p['questions']}, discordant {p['b_only']}:{p['a_only']}, "
            f"p={p['sign_test_p']:.3g}; seconds x{p['seconds_ratio_b_over_a']:.2f} "
            f"{tuple(round(x, 2) for x in p['seconds_ratio_ci'])}; cost "
            f"x{p['cost_ratio_b_over_a']:.2f} {tuple(round(x, 2) for x in p['cost_ratio_ci'])}"
        )


if __name__ == "__main__":
    main(sys.argv[1:])
