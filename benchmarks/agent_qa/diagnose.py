"""Where each arm's tokens and dollars go, read from the raw transcripts.

    uv run python benchmarks/agent_qa/diagnose.py <runs_dir> [out.tsv]

``<runs_dir>`` holds ``<arm>/<qid>-<rep>.jsonl`` stream-json transcripts
(``run.py`` writes them). For every run, from the transcript alone:

- the tools the init event lists, and the context the first API call already
  carries (system prompt + tool definitions: its input + cache-write +
  cache-read tokens);
- turns (``num_turns``) and tool calls;
- characters the model wrote into tool inputs (payload copied into a call is
  paid as output tokens) and characters tools returned (paid as context on
  every later turn);
- tokens by kind, from the result's ``modelUsage`` (the per-call ``usage`` on
  the result event covers the last call only), and the dollars each kind
  cost. The prices are fitted, not assumed: a least-squares fit of
  ``costUSD`` on the four token counts must reproduce every run's cost to
  within 1e-9, or the script stops.

Prints the per-arm medians and totals; writes one row per run to ``out.tsv``
when given.
"""

from __future__ import annotations

import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Any

import numpy as np

KINDS = ("input", "cache_write", "cache_read", "output")
# The server's inline-payload parameters; each now has a *_path twin.
INLINE = {"fasta", "newick_a", "newick_b"}
MECHANISM = (
    "used_path", "used_inline", "used_map_positions", "used_models", "used_read",
)  # fmt: skip


def _events(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(ln) for ln in path.read_text().splitlines() if ln.startswith("{")
    ]


def run_row(path: Path) -> dict[str, Any] | None:
    ev = _events(path)
    result = next((e for e in reversed(ev) if e.get("type") == "result"), None)
    if result is None or not result.get("modelUsage"):
        return None
    init = next((e for e in ev if e.get("subtype") == "init"), {})
    tools = init.get("tools") or []
    mu: dict[str, int] = {}
    for m in result["modelUsage"].values():
        for k, src in (
            ("input", "inputTokens"),
            ("cache_write", "cacheCreationInputTokens"),
            ("cache_read", "cacheReadInputTokens"),
            ("output", "outputTokens"),
            ("thinking", "thinkingTokens"),
        ):
            mu[k] = mu.get(k, 0) + int(m.get(src) or 0)
    first = next(
        (e["message"]["usage"] for e in ev if e.get("type") == "assistant"), {}
    )
    base = sum(
        int(first.get(k) or 0)
        for k in (
            "input_tokens",
            "cache_creation_input_tokens",
            "cache_read_input_tokens",
        )
    )
    in_chars: dict[str, int] = {}
    calls: dict[str, int] = {}
    params: set[str] = set()
    seen: set[str] = set()
    for e in ev:
        if e.get("type") != "assistant":
            continue
        for b in e["message"].get("content", []):
            if b.get("type") == "tool_use" and b.get("id") not in seen:
                seen.add(b.get("id"))
                name = b["name"].replace("mcp__phylokit__", "")
                calls[name] = calls.get(name, 0) + 1
                in_chars[name] = in_chars.get(name, 0) + len(
                    json.dumps(b.get("input", {}))
                )
                if name != b["name"]:
                    params |= {
                        f"{name}.{k}"
                        for k, v in b.get("input", {}).items()
                        if v is not None
                    }
    out_chars: dict[str, int] = {}
    for e in ev:
        if e.get("type") != "user":
            continue
        for b in (e.get("message") or {}).get("content", []):
            if isinstance(b, dict) and b.get("type") == "tool_result":
                out_chars["all"] = out_chars.get("all", 0) + len(
                    json.dumps(b.get("content", ""))
                )
    return {
        "arm": path.parent.name,
        "run": path.stem,
        "n_tools": len(tools),
        "n_mcp_tools": sum(t.startswith("mcp__") for t in tools),
        "base_context": base,
        "turns": result.get("num_turns"),
        "tool_calls": sum(calls.values()),
        "calls": json.dumps(calls, sort_keys=True),
        "tool_input_chars": sum(in_chars.values()),
        "tool_input_chars_by_tool": json.dumps(in_chars, sort_keys=True),
        "tool_output_chars": out_chars.get("all", 0),
        # Which server parameters this run passed (mechanism readout).
        "used_path": any(p.endswith("_path") for p in params),
        "used_inline": any(p.split(".")[1] in INLINE for p in params),
        "used_map_positions": "map_positions" in calls,
        "used_models": any(p.endswith(".models") for p in params),
        "used_read": "Read" in calls,
        "params": json.dumps(sorted(params)),
        **{f"tok_{k}": mu.get(k, 0) for k in (*KINDS, "thinking")},
        "cost_usd": float(result["total_cost_usd"]),
    }


def fit_prices(rows: list[dict[str, Any]]) -> np.ndarray:
    x = np.array([[r[f"tok_{k}"] for k in KINDS] for r in rows], dtype=float)
    y = np.array([r["cost_usd"] for r in rows])
    price, *_ = np.linalg.lstsq(x, y, rcond=None)
    worst = float(np.abs(x @ price - y).max())
    if worst > 1e-9:
        raise SystemExit(
            f"price fit leaves a residual of ${worst:.3g}; costs are not linear in the counts"
        )
    return price


def main() -> None:
    runs = Path(sys.argv[1])
    rows = [r for p in sorted(runs.glob("*/*.jsonl")) if (r := run_row(p))]
    if not rows:
        raise SystemExit(f"no transcripts under {runs}")
    price = fit_prices(rows)
    print(
        "fitted $ per million tokens:",
        dict(zip(KINDS, (price * 1e6).round(4).tolist())),
    )
    for r in rows:
        for k, p in zip(KINDS, price):
            r[f"usd_{k}"] = r[f"tok_{k}"] * p
    arms = sorted({r["arm"] for r in rows})
    cols = [
        "n_tools", "base_context", "turns", "tool_calls", "tool_input_chars",
        "tool_output_chars", "tok_output", "tok_thinking", "tok_cache_write",
        "tok_cache_read", "cost_usd",
    ]  # fmt: skip
    print("\nmedian per run")
    print("arm   n   " + " ".join(f"{c:>17}" for c in cols))
    for arm in arms:
        rs = [r for r in rows if r["arm"] == arm]
        print(f"{arm:5} {len(rs):3} " + " ".join(
            f"{statistics.median(r[c] for r in rs):>17.4g}" for c in cols))  # fmt: skip
    print("\nshare of the arm's total cost, by token kind (total $)")
    for arm in arms:
        rs = [r for r in rows if r["arm"] == arm]
        tot = sum(r["cost_usd"] for r in rs)
        parts = ", ".join(
            f"{k} {sum(r[f'usd_{k}'] for r in rs) / tot:.0%}" for k in KINDS
        )
        print(f"{arm:5} ${tot:7.2f}: {parts}")
    print("\nmechanism: runs (of the arm's runs) passing each kind of argument")
    for arm in arms:
        rs = [r for r in rows if r["arm"] == arm]
        counts = (f"{m[5:]} {sum(r[m] for r in rs)}/{len(rs)}" for m in MECHANISM)
        print(f"{arm:5} " + ", ".join(counts))
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
