"""Work-count ratchet: how much work each tool does on a fixed tiny input.

Counts, never milliseconds. A wall-clock budget needs headroom to survive a
slow runner, and any regression smaller than the headroom passes. These counts
are the same on every run, so the budget is the baseline itself:

- ``python_calls``: calls into functions under ``src/phylokit_mcp`` (generator
  resumptions included), keyed ``relpath:qualname``, taken with
  ``sys.setprofile``. Nothing in the package is edited to count them.
- ``engine_calls``: calls the package makes directly into piqtree's Python
  API, keyed ``piqtree:qualname``: one ``build_tree`` per tree search (the ML
  tree plus one per bootstrap replicate), ``model_finder``,
  ``simulate_alignment``, and small accessors such as ``Model.__str__``. What
  piqtree does inside them is not counted; IQ-TREE's search reads the wall
  clock.
- ``subprocesses``: processes the package launches (MAFFT), keyed by command.
  MAFFT's own children are not counted.

cogent3 calls are deliberately not gated: they are many small accessor calls
whose number is cogent3's implementation detail, so a cogent3 release would
move the gate with no change here.

The total of each counter is the gate; the per-key tables are the diagnosis
printed when it fails, and are stored in baseline.json for that purpose.

Each workload runs in a fresh interpreter, so no cache or lazy import carried
over from another workload changes its count, and returns a witness: an output
only that workload produces (a clade list, an RF distance, a FASTA hash), so a
harness that ran the wrong thing cannot pass by echoing the name.

    uv run python benchmarks/perf/ratchet.py                  # check against baseline.json
    uv run python benchmarks/perf/ratchet.py --update         # lock in a drop
    uv run python benchmarks/perf/ratchet.py --update --allow-increase   # accept a rise

A total that differs from the baseline in EITHER direction fails the check: a
rise is a regression, and a drop is a gain that must be locked in with
``--update`` in the same change, or it can be given back later unnoticed.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
BASELINE = HERE / "baseline.json"
FIXTURES = HERE / "fixtures.json"
WORKLOADS = (
    "infer_tree",
    "select_substitution_model",
    "compare_trees",
    "simulate_alignment",
    "align_sequences",
    "capabilities",
)
COUNTERS = ("python_calls", "engine_calls", "subprocesses")


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def _workload(name: str, server: Any, fx: dict[str, str]) -> dict[str, Any]:
    """Run one workload and return its witness."""
    if name == "infer_tree":
        r = server.infer_tree(fx["aligned_fasta"], model="JC", replicates=20, seed=1)
        return {
            "clades": sorted(c["clade"] for c in r["support"]["clades"]),
            "replicates": r["support"]["replicates"],
        }
    if name == "select_substitution_model":
        r = server.select_substitution_model(
            fx["aligned_fasta"], criterion="BIC", top_n=3
        )
        return {
            "n_models_compared": r["n_models_compared"],
            "n_ranked": len(r["ranking"]),
        }
    if name == "compare_trees":
        r = server.compare_trees(fx["newick_a"], fx["newick_b"])
        return {"rf": r["robinson_foulds"], "only_in_a": r["only_in_a"]}
    if name == "simulate_alignment":
        r = server.simulate_alignment(
            fx["simulate_newick"], model="JC", length=200, seed=1
        )
        return {"fasta_sha": _sha(r["fasta"]), "n_sites": r["alignment"]["n_sites"]}
    if name == "align_sequences":
        r = server.align_sequences(fx["unaligned_fasta"])
        return {"fasta_sha": _sha(r["fasta"]), "n_sites": r["alignment"]["n_sites"]}
    if name == "capabilities":
        r = server.capabilities()
        return {
            "n_substitution_models": r["n_substitution_models"],
            "mafft_found": r["aligner_version"] is not None,
        }
    raise SystemExit(f"unknown workload {name!r}; known: {', '.join(WORKLOADS)}")


def child(name: str) -> dict[str, Any]:
    """Measure one workload in this (fresh) process."""
    if name not in WORKLOADS:
        raise SystemExit(f"unknown workload {name!r}; known: {', '.join(WORKLOADS)}")
    # The server (and so engine, which pins threads) is imported before
    # piqtree, by import_module so no import sorter can reorder the two;
    # measure() also sets OMP_NUM_THREADS=1.
    import phylokit_mcp
    from phylokit_mcp import server

    piqtree = importlib.import_module("piqtree")
    # Imported up front so that module bodies never run inside the counted
    # region; they are not package code, but this keeps the region to the call.
    for mod in (
        "cogent3.core.alphabet",
        "cogent3.core.tree",
        "cogent3.parse.newick",
        "piqtree.model",
    ):
        importlib.import_module(mod)

    pkg_root = str(Path(phylokit_mcp.__file__).resolve().parent) + os.sep
    piq_root = str(Path(piqtree.__file__).resolve().parent) + os.sep
    fx = json.loads(FIXTURES.read_text())
    calls: collections.Counter[str] = collections.Counter()
    engine: collections.Counter[str] = collections.Counter()
    procs: collections.Counter[str] = collections.Counter()

    def profile(frame: Any, event: str, _arg: Any) -> None:
        if event != "call":
            return
        code = frame.f_code
        path = code.co_filename
        if path.startswith(pkg_root):
            calls[f"{path[len(pkg_root) :]}:{code.co_qualname}"] += 1
        elif path.startswith(piq_root):
            back = frame.f_back
            if back is not None and back.f_code.co_filename.startswith(pkg_root):
                engine[f"piqtree:{code.co_qualname}"] += 1

    real_init = subprocess.Popen.__init__

    def counting_init(self: Any, args: Any, *a: Any, **kw: Any) -> None:
        argv = [args] if isinstance(args, str | os.PathLike) else list(args)
        key = os.path.basename(str(argv[0]))
        if "--version" in argv[1:]:
            key += " --version"
        procs[key] += 1
        real_init(self, args, *a, **kw)

    subprocess.Popen.__init__ = counting_init  # type: ignore[method-assign]
    sys.setprofile(profile)
    try:
        witness = _workload(name, server, fx)
    finally:
        sys.setprofile(None)
        subprocess.Popen.__init__ = real_init  # type: ignore[method-assign]
    return {
        "witness": witness,
        "python_calls": dict(sorted(calls.items())),
        "engine_calls": dict(sorted(engine.items())),
        "subprocesses": dict(sorted(procs.items())),
    }


def measure(names: tuple[str, ...] = WORKLOADS) -> dict[str, dict[str, Any]]:
    """Every workload, each in its own interpreter."""
    out: dict[str, dict[str, Any]] = {}
    env = dict(
        os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1"
    )
    for name in names:
        done = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--child", name],
            capture_output=True,
            text=True,
            env=env,
            timeout=600,
            check=False,
        )
        if done.returncode != 0:
            raise RuntimeError(
                f"workload {name} exited {done.returncode}:\n{done.stderr[-4000:]}"
            )
        out[name] = json.loads(done.stdout.strip().splitlines()[-1])
    return out


def totals(row: dict[str, Any]) -> dict[str, int]:
    return {c: sum(row[c].values()) for c in COUNTERS}


def _deltas(old: dict[str, int], new: dict[str, int], top: int) -> list[str]:
    keys = set(old) | set(new)
    diff = [(new.get(k, 0) - old.get(k, 0), k) for k in keys]
    diff = [d for d in diff if d[0]]
    diff.sort(key=lambda d: (-abs(d[0]), d[1]))
    return [
        f"      {d:+d}  {k}  ({old.get(k, 0)} -> {new.get(k, 0)})"
        for d, k in diff[:top]
    ]


def compare(
    baseline: dict[str, dict[str, Any]],
    current: dict[str, dict[str, Any]],
    top: int = 10,
) -> list[str]:
    """Every way ``current`` differs from ``baseline``, as messages; [] if equal.

    A total is compared exactly. Each message names the workload, the counter,
    both totals, and the top per-key deltas that make up the difference.
    """
    problems: list[str] = []
    for name, row in current.items():
        if name not in baseline:
            problems.append(f"{name}: not in baseline.json; run --update to add it")
            continue
        base = baseline[name]
        if row["witness"] != base["witness"]:
            problems.append(
                f"{name}: witness {row['witness']} != baseline {base['witness']}: the "
                "workload did not produce the output it is recorded with"
            )
        bt, ct = totals(base), totals(row)
        for c in COUNTERS:
            if ct[c] == bt[c]:
                continue
            kind = (
                "ROSE (regression; --update --allow-increase to accept)"
                if ct[c] > bt[c]
                else "FELL (lock it in with --update)"
            )
            lines = [f"{name}: {c} {kind}: {bt[c]} -> {ct[c]} ({ct[c] - bt[c]:+d})"]
            lines += _deltas(base[c], row[c], top)
            problems.append("\n".join(lines))
    return problems


def summary(rows: dict[str, dict[str, Any]]) -> str:
    out = [
        f"{'workload':28} {'python_calls':>12} {'engine_calls':>12} {'subprocesses':>12}"
    ]
    for name, row in rows.items():
        t = totals(row)
        out.append(
            f"{name:28} {t['python_calls']:>12} {t['engine_calls']:>12} {t['subprocesses']:>12}"
        )
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--child", help=argparse.SUPPRESS)
    ap.add_argument("--only", default="", help="comma-separated workloads")
    ap.add_argument("--update", action="store_true", help="write baseline.json")
    ap.add_argument(
        "--allow-increase", action="store_true", help="with --update, accept a rise"
    )
    a = ap.parse_args(argv)
    if a.child:
        print(json.dumps(child(a.child)))
        return 0
    names = tuple(a.only.split(",")) if a.only else WORKLOADS
    current = measure(names)
    print(summary(current))
    baseline = json.loads(BASELINE.read_text()) if BASELINE.exists() else {}
    problems = compare(baseline, current)
    if not a.update:
        for p in problems:
            print(p)
        return 1 if problems else 0
    rises = [
        f"{n}: {c} {totals(baseline[n])[c]} -> {totals(r)[c]}"
        for n, r in current.items()
        if n in baseline
        for c in COUNTERS
        if totals(r)[c] > totals(baseline[n])[c]
    ]
    if rises and not a.allow_increase:
        print("refusing to raise the baseline without --allow-increase:")
        print("\n".join(problems))
        return 1
    witness_changed = [
        n
        for n, r in current.items()
        if n in baseline and r["witness"] != baseline[n]["witness"]
    ]
    if witness_changed and not a.allow_increase:
        print(
            f"refusing: witness changed for {witness_changed}; pass --allow-increase if intended"
        )
        return 1
    merged = {**baseline, **current}
    BASELINE.write_text(
        json.dumps({n: merged[n] for n in WORKLOADS if n in merged}, indent=1) + "\n"
    )
    print(f"wrote {BASELINE.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
