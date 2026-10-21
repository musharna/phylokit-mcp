"""Tables and statistics for the agent benchmark, from results.jsonl.

    uv run python benchmarks/agent_qa/analyze.py [results.jsonl[,results_read.jsonl]] [outdir]

The first argument is one results file or several, comma-separated (the read
control has its own); a run (arm, question, rep) found twice is an error.
Writes, next to the first results file:

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
- with ``--substitute=<results.jsonl>`` (another model's rows, e.g. the
  2026-10-18 Opus re-run): ``substitute.tsv``, this file's ``mcp`` arm against
  that file's ``alt`` arm (``other_alt``), paired the same way. Not concurrent:
  report it as a secondary comparison only.

PREREG_2026-10-20_model-tiers.md's primary endpoint per model tier is the
paired accuracy columns of ``paired.tsv`` (alt vs mcp): questions right in
every run, McNemar's exact test on the discordant questions, and a bootstrap
95% interval of the difference; ``tier_rule`` combines the tiers (Holm).

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
# The kinds of question the read control (Read only, nothing to compute with)
# did NOT solve by hand. Chosen AFTER the read arm ran, from its per-kind
# results: on RF distance, JC69 distance and homology it worked the answer out
# from the file (splits listed, differences counted, flanking motifs matched),
# so those kinds cannot separate a tool from no tool. A post-hoc subset, and
# labelled as one wherever it is reported.
NEEDS_COMPUTE = ("sister_taxon", "model_selection", "informative_sites")
SUBSET = "needs_compute"


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


def server_calls(r: dict[str, Any]) -> int:
    """Calls this run made to the server's tools (``mcp__phylokit__*``)."""
    return sum(t.startswith("mcp__phylokit__") for t in r["tool_calls"])


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
        if arm == "mcp":
            # The mechanism: did the server arm call the server at all?
            n_srv = [server_calls(r) for r in mine]
            row.update(
                server_calls_mean=float(np.mean(n_srv)),
                server_calls_median=float(np.median(n_srv)),
                runs_without_server_call=sum(n == 0 for n in n_srv),
            )
        out.append(row)
        # No subset row on results that hold none of its kinds (it is undefined).
        sq = sorted({r["qid"] for r in mine if r["category"] in NEEDS_COMPUTE})
        if sq:
            srow = _tally("arm_x_subset", arm, SUBSET, [r for q in sq for r in by_q[q]])
            srow.update(
                accuracy=acc(sq),
                accuracy_ci=bootstrap(sq, acc, rng),
                questions=len(sq),
                questions_right_every_run=sum(all(map(_ok, by_q[q])) for q in sq),
                questions_right_every_run_ci=tuple(
                    len(sq) * x for x in bootstrap(sq, all_right, rng)
                ),
            )
            out.append(srow)
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
    rows: list[dict[str, Any]],
    rng: np.random.Generator,
    pairs: list[tuple[str, str]] = PAIRS,
) -> list[dict[str, Any]]:
    out = []
    arms = {r["arm"] for r in rows}
    cat = {r["qid"]: r["category"] for r in rows}
    for subset in ("all", SUBSET):
        for a, b in pairs:
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
            if subset == SUBSET:
                qids = [q for q in qids if cat[q] in NEEDS_COMPUTE]
                if not qids:
                    continue
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
            costlier = sum(
                np.median([r["cost_usd"] for r in qb[q]])
                > np.median([r["cost_usd"] for r in qa[q]])
                for q in qids
            )

            def diff(qs: list[str], ra: dict = ra, rb: dict = rb) -> float:
                return float(sum(rb[q] for q in qs) - sum(ra[q] for q in qs))

            # Its own generator, seeded alike, so adding this interval leaves
            # every interval drawn from ``rng`` (the 10-16 and 10-18 tables) as
            # it was.
            diff_ci = bootstrap(qids, diff, np.random.default_rng(SEED))
            out.append(
                {
                    "a": a,
                    "b": b,
                    "subset": subset,
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
                    "questions_b_costlier": costlier,
                    "cost_sign_test_p": sign_test(costlier, len(qids) - costlier),
                    "right_every_run_diff_b_minus_a": b_only - a_only,
                    "right_every_run_diff_ci": diff_ci,
                }
            )
    return out


def holm(ps: dict[str, float]) -> dict[str, float]:
    """Holm-adjusted p values (step-down, monotone, capped at 1)."""
    order = sorted(ps, key=ps.__getitem__)
    adj: dict[str, float] = {}
    run = 0.0
    for i, k in enumerate(order):
        run = max(run, min(1.0, (len(order) - i) * ps[k]))
        adj[k] = run
    return adj


def holm_family(ps: dict[str, float], expected: int = 2) -> dict[str, float]:
    """Holm over the planned family: missing tiers are padded with p = 1, so a
    partial family is adjusted as conservatively as the full one could be."""
    pad = {f"_missing{i}": 1.0 for i in range(expected - len(ps))}
    adj = holm(ps | pad)
    return {t: adj[t] for t in ps}


def tier_rule(
    tiers: dict[str, dict[str, Any]], expected: int = 2, alpha: float = 0.05
) -> tuple[str, dict[str, str]]:
    """PREREG_2026-10-20_model-tiers' rule for the server (b) against alt (a)
    on accuracy, from each tier's ``paired.tsv`` row (alt vs mcp, all questions).

    Per tier, with McNemar's exact p Holm-adjusted across the tiers: a LIFT is
    adjusted p < alpha, more questions right in every run only with the server
    (b_only > a_only) and the bootstrap 95% CI of the difference wholly above
    0; a DROP is the mirror image; anything else is NONE. Over the tiers: WIN
    is a lift at one tier or more and no drop; LOSS a drop and no lift; MIXED
    both; TIE neither. Until all ``expected`` tiers are in, the verdict is
    PENDING (the per-tier calls use the Holm factor of the full family, so a
    lone tier is never judged on its raw p).
    """
    adj = holm_family({t: p["sign_test_p"] for t, p in tiers.items()}, expected)
    calls: dict[str, str] = {}
    for t, p in tiers.items():
        lo, hi = p["right_every_run_diff_ci"]
        sig = adj[t] < alpha
        if sig and p["b_only"] > p["a_only"] and lo > 0:
            calls[t] = "LIFT"
        elif sig and p["a_only"] > p["b_only"] and hi < 0:
            calls[t] = "DROP"
        else:
            calls[t] = "NONE"
    if len(tiers) < expected:
        return "PENDING", calls
    lift, drop = "LIFT" in calls.values(), "DROP" in calls.values()
    verdict = {(True, False): "WIN", (False, True): "LOSS", (True, True): "MIXED"}
    return verdict.get((lift, drop), "TIE"), calls


def win_rule(p: dict[str, Any], margin: int = 2, alpha: float = 0.05) -> str:
    """PREREG_2026-10-18's rule for b (the server) against a (alt), on cost.

    Non-inferior: b's questions right in every run >= a's - margin. Better on the
    primary endpoint (per-question median cost): the sign test on questions where
    b costs more vs less has p < alpha AND the bootstrap 95% CI of the geometric
    mean b/a cost ratio lies wholly below 1. Win needs both. Loss: b not
    non-inferior, or b significantly costlier (p < alpha, CI wholly above 1).
    Anything else is a tie.
    """
    non_inferior = p["b_right_every_run"] >= p["a_right_every_run"] - margin
    lo, hi = p["cost_ratio_ci"]
    sig = p["cost_sign_test_p"] < alpha
    cheaper = sig and hi < 1 and p["questions_b_costlier"] < p["questions"] / 2
    costlier = sig and lo > 1 and p["questions_b_costlier"] > p["questions"] / 2
    if non_inferior and cheaper:
        return "WIN"
    if not non_inferior or costlier:
        return "LOSS"
    return "TIE"


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


def read_tier(path: Path) -> dict[str, Any]:
    """A tier's alt-vs-mcp, all-questions paired row, recomputed at full
    precision from its results.jsonl (paired.tsv rounds to 4 digits). The
    difference CI has its own seeded generator, so it equals paired.tsv's."""
    rows = [json.loads(ln) for ln in path.read_text().splitlines()]
    (p,) = [
        p
        for p in paired_rows(rows, np.random.default_rng(SEED), [("alt", "mcp")])
        if p["subset"] == "all"
    ]
    return p


def tiers_main(spec: str) -> None:
    """``--tiers=haiku=<results.jsonl>,sonnet=<results.jsonl>``: Holm and the
    verdict; the Holm p printed is the one ``tier_rule`` decides on."""
    tiers = {
        name: read_tier(Path(path))
        for name, path in (x.split("=", 1) for x in spec.split(","))
    }
    adj = holm_family({t: p["sign_test_p"] for t, p in tiers.items()})
    verdict, calls = tier_rule(tiers)
    for t, p in tiers.items():
        print(
            f"{t}: mcp-only {p['b_only']}, alt-only {p['a_only']}, raw p="
            f"{p['sign_test_p']:.3g}, Holm p={adj[t]:.3g}, diff CI "
            f"{p['right_every_run_diff_ci']}: {calls[t]}"
        )
    print(f"tier rule (PREREG_2026-10-20_model-tiers): {verdict}")


def as_written(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The rows judged by the leak rule as written before 2026-10-21 (the
    secondary result); a row without that judgement is an error."""
    missing = [r for r in rows if "valid_as_written" not in r]
    if missing:
        raise SystemExit(f"{len(missing)} rows lack valid_as_written: run.py --rescore")
    return [r | {"valid": r["valid_as_written"]} for r in rows]


def main(argv: list[str]) -> None:
    written = "--as-written" in argv
    argv = [x for x in argv if x != "--as-written"]
    tiers = [x.split("=", 1)[1] for x in argv if x.startswith("--tiers=")]
    if tiers:
        tiers_main(tiers[0])
        return
    sub = [x.split("=", 1)[1] for x in argv if x.startswith("--substitute=")]
    argv = [x for x in argv if not x.startswith("--substitute=")]
    paths = [Path(p) for p in argv[0].split(",")] if argv else [HERE / "results.jsonl"]
    outdir = Path(argv[1]) if len(argv) > 1 else paths[0].parent
    rows = [json.loads(ln) for p in paths for ln in p.read_text().splitlines()]
    if written:
        rows = as_written(rows)
        print("validity: the leak rule AS WRITTEN (secondary)")
    keys = [(r["arm"], r["qid"], r["rep"]) for r in rows]
    if len(set(keys)) != len(keys):
        raise SystemExit(f"a run appears twice across {[str(p) for p in paths]}")
    questions = json.loads((HERE / "questions.json").read_text())["questions"]
    rng = np.random.default_rng(SEED)
    base = baselines(questions)
    summary = arm_rows(rows, rng) + baseline_rows(base, rng)
    write_tsv(outdir / "summary.tsv", summary)
    paired = paired_rows(rows, rng)
    write_tsv(outdir / "paired.tsv", paired)
    write_tsv(outdir / "baselines.tsv", base)
    for r in summary:
        if r["level"] in ("arm", "arm_x_subset", "baseline"):
            lo, hi = r["accuracy_ci"]
            print(
                f"{r['level']:12} {r['arm']:38} {r['category']:13} {r['right']:>3}/{r['n']:<3} "
                f"{r['accuracy']:.3f} [{lo:.3f}, {hi:.3f}] invalid={r['invalid']}"
            )
    for p in paired:
        print(
            f"[{p['subset']}] {p['b']} vs {p['a']}: right every run {p['b_right_every_run']} vs "
            f"{p['a_right_every_run']} of {p['questions']}, discordant {p['b_only']}:{p['a_only']}, "
            f"p={p['sign_test_p']:.3g}; seconds x{p['seconds_ratio_b_over_a']:.2f} "
            f"{tuple(round(x, 2) for x in p['seconds_ratio_ci'])}; cost "
            f"x{p['cost_ratio_b_over_a']:.2f} {tuple(round(x, 2) for x in p['cost_ratio_ci'])}, "
            f"costlier on {p['questions_b_costlier']} of {p['questions']} "
            f"(p={p['cost_sign_test_p']:.3g})"
        )
        if (p["a"], p["b"], p["subset"]) == ("alt", "mcp", "all"):
            print(f"win rule (PREREG_2026-10-18, cost endpoint): {win_rule(p)}")
            lo, hi = p["right_every_run_diff_ci"]
            print(
                f"accuracy (PREREG_2026-10-20 tier endpoint): mcp - alt = "
                f"{p['right_every_run_diff_b_minus_a']:+d} questions right in every run "
                f"(95% CI {lo:+.0f} to {hi:+.0f}), McNemar exact raw p="
                f"{p['sign_test_p']:.3g}; this tier alone: {tier_rule({'this': p})}"
            )
    for r in summary:
        if r["level"] == "arm" and r["arm"] == "mcp":
            print(
                f"mechanism: server calls per mcp run mean {r['server_calls_mean']:.2f}, "
                f"median {r['server_calls_median']:g}; runs with no server call "
                f"{r['runs_without_server_call']} of {r['n']}"
            )
    if sub:
        other = [
            json.loads(ln) | {"arm": "other_alt"}
            for ln in Path(sub[0]).read_text().splitlines()
            if json.loads(ln)["arm"] == "alt"
        ]
        mine = [r for r in rows if r["arm"] == "mcp"]
        srows = paired_rows(
            mine + other, np.random.default_rng(SEED), [("other_alt", "mcp")]
        )
        write_tsv(outdir / "substitute.tsv", srows)
        for p in srows:
            lo, hi = p["right_every_run_diff_ci"]
            print(
                f"substitute [{p['subset']}] (NOT concurrent) mcp vs {sub[0]} alt: right "
                f"every run {p['b_right_every_run']} vs {p['a_right_every_run']} of "
                f"{p['questions']}, discordant {p['b_only']}:{p['a_only']}, "
                f"p={p['sign_test_p']:.3g}, diff CI ({lo:+.0f}, {hi:+.0f}); cost "
                f"x{p['cost_ratio_b_over_a']:.2f} "
                f"{tuple(round(x, 2) for x in p['cost_ratio_ci'])}"
            )


if __name__ == "__main__":
    main(sys.argv[1:])
