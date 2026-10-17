"""Ask every question to Claude Code under each arm, and score the replies.

Arms (same model, same prompt, same system prompt, isolated: no user settings,
hooks or other MCP servers, a fresh empty working directory per run holding
only that question's data files):

- ``bare``: no tools at all.
- ``alt``:  what a user has without this server: Bash, Read, WebSearch and
  WebFetch, with a Python that has piqtree, cogent3 and numpy (but not
  phylokit-mcp) and MAFFT on PATH.
- ``mcp``:  phylokit-mcp (pinned PyPI release via uvx) only, MAFFT on PATH.

Each finished run is appended to ``results.jsonl`` with its score; a rerun
skips (arm, question, rep) triples already there, so an interrupted sweep
resumes. Raw stream-json transcripts go to ``runs/<arm>/<qid>-<rep>.jsonl``.

    uv run python benchmarks/agent_qa/run.py --setup-alt     # once: builds .alt-env
    uv run python benchmarks/agent_qa/run.py --reps 3 --workers 2
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from score import score

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
MODEL = "claude-opus-5-5"
SERVER = "phylokit-mcp==0.6.0"
RUN_TIMEOUT_S = 900
# The alt arm's Python: the repo's locked piqtree/cogent3/numpy, without the
# server. The worktree's own .venv has phylokit_mcp installed (editable), and a
# run that imported it would be the mcp arm by another route.
ALT_ENV = HERE / ".alt-env"
# A directory holding only a `mafft` symlink, so putting MAFFT on PATH does not
# also put whatever else lives beside it (a conda env's python, for one).
BIN = HERE / ".bin"
SYSTEM_PATH = "/usr/local/bin:/usr/bin:/bin"
# Every arm gets this in place of Claude Code's own system prompt. That prompt
# describes tools even when --tools "" removes them all; in the sibling
# plant-genomics-mcp benchmark's pilot the no-tools arm wrote fake tool calls
# as text instead of answering. Tools still reach the model through the API.
SYSTEM = "You are a helpful assistant. Answer the user's question."
PROMPT = (
    "{question}\n\nGive the answer as {answer_format}. Keep your reply short, and end it "
    "with one line of the form\nANSWER: <answer>\nIf you cannot determine the answer, end "
    "with\nANSWER: UNKNOWN"
)


def mafft_bin() -> Path:
    """BIN/mafft, linked to the MAFFT on PATH (or $MAFFT); fails loud if none."""
    target = os.environ.get("MAFFT") or shutil.which("mafft")
    if not target:
        raise SystemExit("MAFFT not found: put it on PATH or set MAFFT=/path/to/mafft")
    BIN.mkdir(exist_ok=True)
    link = BIN / "mafft"
    if link.is_symlink() and Path(os.readlink(link)) != Path(target):
        link.unlink()
    if not link.exists():
        link.symlink_to(target)
    return BIN


def setup_alt() -> None:
    """Build ALT_ENV from uv.lock: every locked dependency, not the project."""
    uv = shutil.which("uv")
    if uv is None:
        raise SystemExit("uv not found")
    reqs = subprocess.run(
        [uv, "export", "--frozen", "--no-dev", "--no-emit-project", "--no-hashes"],
        cwd=REPO, capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip
    req_file = HERE / ".alt-env.requirements.txt"
    req_file.write_text(reqs)
    subprocess.run([uv, "venv", "--python", sys.executable, str(ALT_ENV)], check=True)
    subprocess.run(
        [
            uv,
            "pip",
            "install",
            "--python",
            str(ALT_ENV / "bin" / "python"),
            "-r",
            str(req_file),
        ],
        check=True,
    )
    check_alt()


def check_alt() -> None:
    """The alt Python has the libraries and NOT the server; both halves checked."""
    py = ALT_ENV / "bin" / "python"
    if not py.exists():
        raise SystemExit(f"{ALT_ENV} missing: run with --setup-alt first")
    ok = subprocess.run(
        [str(py), "-c", "import piqtree, cogent3, numpy; print(piqtree.__version__)"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    if ok.returncode != 0:
        raise SystemExit(f"alt env cannot import piqtree/cogent3: {ok.stderr}")
    leak = subprocess.run(
        [str(py), "-c", "import phylokit_mcp"], capture_output=True, check=False
    )
    if leak.returncode == 0:
        raise SystemExit("alt env can import phylokit_mcp; it must not")


def arm_setup(arm: str, tmp: Path) -> tuple[list[str], dict[str, str]]:
    """The claude flags and the environment for one arm."""
    empty = tmp.parent / f"{tmp.name}-no-mcp.json"
    empty.write_text('{"mcpServers": {}}')
    env = dict(os.environ)
    if arm == "bare":
        return ["--tools", "", "--strict-mcp-config", "--mcp-config", str(empty)], env
    if arm == "alt":
        env.pop("VIRTUAL_ENV", None)
        env["PATH"] = f"{ALT_ENV / 'bin'}:{mafft_bin()}:{SYSTEM_PATH}"
        tools = "Bash,Read,WebSearch,WebFetch"
        return [
            "--tools", tools, "--allowedTools", tools,
            "--strict-mcp-config", "--mcp-config", str(empty),
        ], env  # fmt: skip
    if arm == "mcp":
        uvx = shutil.which("uvx")
        if uvx is None:
            raise SystemExit("uvx not found")
        path = f"{mafft_bin()}:{Path(uvx).parent}:{SYSTEM_PATH}"
        env.pop("VIRTUAL_ENV", None)
        env["PATH"] = path
        cfg = tmp.parent / f"{tmp.name}-mcp.json"
        cfg.write_text(
            json.dumps(
                {
                    "mcpServers": {
                        "phylokit": {
                            "command": uvx,
                            "args": [SERVER],
                            "env": {"PATH": path},
                        }
                    }
                }
            )
        )
        return [
            "--tools", "", "--allowedTools", "mcp__phylokit",
            "--strict-mcp-config", "--mcp-config", str(cfg),
        ], env  # fmt: skip
    raise ValueError(arm)


def _tool_inputs(events: list[dict[str, Any]]) -> list[str]:
    return [
        json.dumps(b.get("input", {}))
        for e in events
        if e.get("type") == "assistant"
        for b in e["message"].get("content", [])
        if b.get("type") == "tool_use"
    ]


def run_one(arm: str, q: dict[str, Any], rep: int) -> dict[str, Any]:
    out = HERE / "runs" / arm / f"{q['id']}-{rep}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    prompt = PROMPT.format(question=q["question"], answer_format=q["answer_format"])
    claude = shutil.which("claude")
    if claude is None:
        raise SystemExit("claude not found on PATH")
    with tempfile.TemporaryDirectory(prefix=f"phyloqa-{arm}-") as d:
        tmp = Path(d) / "cwd"
        tmp.mkdir()
        for name, text in q.get("files", {}).items():
            (tmp / name).write_text(text)
        flags, env = arm_setup(arm, tmp)
        cmd = [
            claude, "-p", prompt, "--model", MODEL, "--setting-sources", "",
            "--output-format", "stream-json", "--verbose", "--system-prompt", SYSTEM,
            *flags,
        ]  # fmt: skip
        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                cmd,
                cwd=tmp,
                env=env,
                capture_output=True,
                text=True,
                timeout=RUN_TIMEOUT_S,
                check=False,
            )
            stdout, rc = proc.stdout, proc.returncode
        except subprocess.TimeoutExpired as exc:
            stdout = (
                (exc.stdout or b"").decode()
                if isinstance(exc.stdout, bytes)
                else exc.stdout or ""
            )
            rc = "timeout"
        seconds = round(time.monotonic() - t0, 1)
    out.write_text(stdout)
    events = [json.loads(ln) for ln in stdout.splitlines() if ln.startswith("{")]
    init = next(
        (e for e in events if e.get("type") == "system" and e.get("subtype") == "init"),
        {},
    )
    result = next((e for e in events if e.get("type") == "result"), {})
    tools = [
        b["name"]
        for e in events
        if e.get("type") == "assistant"
        for b in e["message"].get("content", [])
        if b.get("type") == "tool_use"
    ]
    inputs = _tool_inputs(events)
    reply = result.get("result") or ""
    return {
        "arm": arm,
        "qid": q["id"],
        "category": q["category"],
        "rep": rep,
        "exit": rc,
        "seconds": seconds,
        "claude_code": init.get("claude_code_version"),
        "model": init.get("model"),
        "mcp_servers": init.get("mcp_servers"),
        "tools_available": init.get("tools") or [],
        "tool_calls": tools,
        # Read from the transcript, not inferred: did a tool call name the
        # server's package, or the benchmark's own key?
        "mentions_phylokit": any("phylokit" in s.lower() for s in inputs),
        "mentions_key": any("questions.json" in s or "agent_qa" in s for s in inputs),
        "cost_usd": result.get("total_cost_usd"),
        "reply": reply,
        **score(q, reply),
    }


def rescore(results: Path, by_id: dict[str, dict[str, Any]]) -> None:
    """Score each stored reply again, so a truth or scorer fix reaches runs
    made before it; prints every verdict that changed."""
    rows = [json.loads(ln) for ln in results.read_text().splitlines()]
    for r in rows:
        new = score(by_id[r["qid"]], r["reply"])
        if new["verdict"] != r["verdict"] or new["wrong_id"] != r["wrong_id"]:
            print(f"{r['arm']:4} {r['qid']} rep{r['rep']} {r['verdict']} -> {new['verdict']}"
                  f" ({new['value']!r})")  # fmt: skip
        r.update(new)
    results.write_text("".join(json.dumps(r) + "\n" for r in rows))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="bare,alt,mcp")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", default="", help="comma-separated question ids")
    ap.add_argument("--results", default=str(HERE / "results.jsonl"))
    ap.add_argument(
        "--setup-alt", action="store_true", help="build the alt arm's Python"
    )
    ap.add_argument(
        "--rescore",
        action="store_true",
        help="re-score every stored reply against the current questions.json, run nothing",
    )
    a = ap.parse_args()
    if a.setup_alt:
        setup_alt()
        return
    questions = json.loads((HERE / "questions.json").read_text())["questions"]
    if a.rescore:
        rescore(Path(a.results), {q["id"]: q for q in questions})
        return
    arms = a.arms.split(",")
    if "alt" in arms:
        check_alt()
    if a.only:
        keep = set(a.only.split(","))
        questions = [q for q in questions if q["id"] in keep]
    results = Path(a.results)
    done = set()
    if results.exists():
        for ln in results.read_text().splitlines():
            r = json.loads(ln)
            done.add((r["arm"], r["qid"], r["rep"]))
    todo = [
        (arm, q, rep)
        for rep in range(1, a.reps + 1)
        for q in questions
        for arm in arms
        if (arm, q["id"], rep) not in done
    ]
    print(f"{len(todo)} runs to do ({len(done)} already in {results.name})", flush=True)
    lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(run_one, *t) for t in todo]
        for n, f in enumerate(as_completed(futures), 1):
            r = f.result()
            with lock, results.open("a") as fh:
                fh.write(json.dumps(r) + "\n")
            print(
                f"[{n}/{len(todo)}] {r['arm']:4} {r['qid']} rep{r['rep']} {r['verdict']:9}"
                f" exit={r['exit']} {r['seconds']}s calls={len(r['tool_calls'])}"
                f" ${r['cost_usd']}",
                flush=True,
            )
    sys.exit(0)


if __name__ == "__main__":
    main()
