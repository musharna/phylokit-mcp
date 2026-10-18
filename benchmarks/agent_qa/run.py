"""Ask every question to Claude Code under each arm, and score the replies.

Arms (same model, same system prompt, isolated: no user settings, hooks or
other MCP servers, a fresh empty working directory per run holding only that
question's data files; the data are NOT in the prompt):

- ``bare``: no tools at all. It cannot open the data file, so declining is the
  right answer for it; it is the floor, not a competitor.
- ``alt``:  what a user has without this server: Bash, Read, WebSearch and
  WebFetch, with a Python that has piqtree, cogent3 and numpy (but not
  phylokit-mcp) and MAFFT on PATH. Its prompt says so (ALT_NOTE): a user who
  installed them knows they did, as the mcp arm learns its tools from their
  descriptions. That note is the only difference between the arms' prompts.
- ``mcp``:  phylokit-mcp (pinned PyPI release via uvx, or the commit named
  by ``--server-spec``) plus Read, MAFFT on
  PATH, and no shell. Read is there because the server takes sequences as a
  string, not a path, so the data must be read before it can be passed. Its
  libraries are constrained to the alt arm's (the repo's uv.lock), so both
  arms run the same piqtree, cogent3 and numpy.
- ``read``: a control, Read only (no shell, no server, no web). It tests the
  claim that the inputs are past what can be worked by eye: it can see every
  byte of the data and nothing to compute with.

One validity rule for every arm (``validity``): the run exited 0 with a
result; the tools the init event lists are exactly the arm's; a non-bare arm
made at least one tool call (every question needs its file); and no tool
input or output touched the key, the benchmark's code, its transcripts, or
another run's directory (``leak_hits``). Invalid rows are kept and reported,
never retried or overwritten; no arm is retried.

Each finished run is appended to ``results.jsonl`` with its score; a rerun
skips (arm, question, rep) triples already there, so an interrupted sweep
resumes. Raw stream-json transcripts go to ``<artifacts>/runs/<arm>/<qid>-<rep>.jsonl``,
where ``<artifacts>`` is ``$PHYLOQA_ARTIFACTS`` or
``~/phylokit-mcp-artifacts/agent_qa``: outside the repo, so a shell in a run
does not sit next to the key. Each run's working directory is a fresh
``phyloqa-<arm>-*`` directory in the system temp dir, whose path must name
nothing the leak scan looks for (``check_run_dir``): a run reads its own
files by absolute path, and a path naming the server or the benchmark would
mark every run invalid.

    uv run python benchmarks/agent_qa/run.py --setup-alt     # once: builds .alt-env
    uv run python benchmarks/agent_qa/run.py --reps 3 --workers 2

``--server-spec git+https://github.com/musharna/phylokit-mcp@<40-hex sha>``
runs the mcp arm on an unreleased commit instead of the PyPI pin. The spec is
written into every result row and into ``--versions`` (a JSON file of what each
arm ran on), and ``check_mcp`` confirms the installed package came from that
commit before any run is spent.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from controls import oracle_check
from score import score

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
MODEL = "claude-opus-5-5"
SERVER = "phylokit-mcp==0.6.0"
# A git spec must name a full commit: a branch or tag can move under a sweep.
GIT_SPEC = re.compile(
    r"^git\+https://github\.com/musharna/phylokit-mcp@([0-9a-f]{40})$"
)
RUN_TIMEOUT_S = 1800
ARTIFACTS = Path(
    os.environ.get(
        "PHYLOQA_ARTIFACTS", Path.home() / "phylokit-mcp-artifacts" / "agent_qa"
    )
)
RUNS = ARTIFACTS / "runs"
# The alt arm's Python: the repo's locked piqtree/cogent3/numpy, without the
# server. The worktree's own .venv has phylokit_mcp installed (editable), and a
# run that imported it would be the mcp arm by another route. It lives with
# the artifacts, outside the repo: in the v2 sweep it sat in this directory,
# and one traceback printed its path, which names the benchmark's directory.
ALT_ENV = ARTIFACTS / "alt-env"
ALT_REQS = ARTIFACTS / "alt-env.requirements.txt"
# A directory holding only a `mafft` symlink, so putting MAFFT on PATH does not
# also put whatever else lives beside it (a conda env's python, for one).
BIN = ARTIFACTS / "bin"
ARMS = ("bare", "alt", "mcp", "read")
# The tools each arm's init event must list, exactly (the mcp arm also needs
# at least one mcp__phylokit__ tool and a connected server).
INTENDED = {
    "bare": set(),
    "alt": {"Bash", "Read", "WebSearch", "WebFetch"},
    "mcp": {"Read"},
    "read": {"Read"},
}
# A tool input or output matching any of these touched what a run must not:
# the key and the benchmark's code (by file or directory name), its
# transcripts and artifacts, the coordinator's job files, or the benchmark's
# directory in the GitHub repo. Another run's working directory is checked
# separately (RUN_DIR), since a run's own path is legitimately in its Reads.
LEAK_PATTERNS = (
    r"questions\.json",
    r"build_questions",
    r"\btruth\.py",
    r"\bcontrols\.py",
    r"results\.jsonl",
    r"agent_qa",
    r"pilot-v2",
    r"phylokit_v1_partial",
    r"\.claude/(?:worktrees|jobs)",
    r"phylokit-mcp-artifacts/agent_qa/runs",
    r"(?:github\.com|githubusercontent\.com)/musharna/phylokit-mcp\S*benchmarks",
)
RUN_DIR = re.compile(r"phyloqa-[a-z]+-[\w]+")
SYSTEM_PATH = "/usr/local/bin:/usr/bin:/bin"
# Every arm gets this in place of Claude Code's own system prompt. That prompt
# describes tools even when --tools "" removes them all; in the sibling
# plant-genomics-mcp benchmark's pilot the no-tools arm wrote fake tool calls
# as text instead of answering. Tools still reach the model through the API.
SYSTEM = "You are a helpful assistant. Answer the user's question."
ALT_NOTE = (
    "\n\nYour shell's python3 has piqtree (the IQ-TREE 3 bindings), cogent3 and "
    "numpy installed, and MAFFT is on PATH."
)
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
    BIN.mkdir(parents=True, exist_ok=True)
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
    req_file = ALT_REQS
    req_file.parent.mkdir(parents=True, exist_ok=True)
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
    check_alt(json.loads((HERE / "questions.json").read_text())["checked_with"])


_VERSIONS = (
    "import importlib.metadata as m, json; "
    "print(json.dumps({p: m.version(p) for p in ('piqtree', 'cogent3', 'numpy')}))"
)


def _expect_engine(
    where: str, versions: dict[str, str], checked: dict[str, str]
) -> None:
    """The key's answerability filters ran on ``checked['piqtree']``; an arm on
    another piqtree is not answering the questions the key was checked on."""
    if versions.get("piqtree") != checked["piqtree"]:
        raise SystemExit(
            f"{where} has piqtree {versions.get('piqtree')}, but questions.json was "
            f"checked with {checked['piqtree']}"
        )


def check_alt(checked: dict[str, str]) -> dict[str, str]:
    """The alt Python has the libraries and NOT the server; both halves checked,
    and its piqtree is the one the key was checked with. Returns its versions."""
    py = ALT_ENV / "bin" / "python"
    if not py.exists():
        raise SystemExit(f"{ALT_ENV} missing: run with --setup-alt first")
    ok = subprocess.run(
        [str(py), "-c", _VERSIONS], capture_output=True, text=True, check=False
    )
    if ok.returncode != 0:
        raise SystemExit(f"alt env cannot import piqtree/cogent3: {ok.stderr}")
    versions = json.loads(ok.stdout)
    _expect_engine("the alt env", versions, checked)
    leak = subprocess.run(
        [str(py), "-c", "import phylokit_mcp"], capture_output=True, check=False
    )
    if leak.returncode == 0:
        raise SystemExit("alt env can import phylokit_mcp; it must not")
    return versions


def check_server_spec(spec: str) -> str | None:
    """The PyPI pin, or a git spec naming a full commit; returns the commit."""
    if spec == SERVER:
        return None
    m = GIT_SPEC.match(spec)
    if not m:
        raise SystemExit(
            f"--server-spec must be {SERVER!r} or "
            f"git+https://github.com/musharna/phylokit-mcp@<40-hex sha>; got {spec!r}"
        )
    return m.group(1)


def _uvx_server_args(spec: str = SERVER) -> list[str]:
    """uvx arguments for the server, its libraries constrained to the alt env's
    (the repo's uv.lock), so the two arms differ in the server, not the stack."""
    reqs = ALT_REQS
    if not reqs.exists():
        raise SystemExit(f"{reqs} missing: run with --setup-alt first")
    return ["--constraints", str(reqs), "--from", spec, "phylokit-mcp"]


_SERVER_ORIGIN = (
    "import importlib.metadata as m, json; d = m.distribution('phylokit-mcp'); "
    "u = d.read_text('direct_url.json'); "
    "print(json.dumps({'phylokit-mcp': d.version, "
    "'commit': json.loads(u).get('vcs_info', {}).get('commit_id') if u else None}))"
)


def check_mcp(checked: dict[str, str], spec: str = SERVER) -> dict[str, str]:
    """The server's own environment, resolved as the mcp arm will start it; for
    a git spec, the installed package must say it came from that commit."""
    uvx = shutil.which("uvx")
    if uvx is None:
        raise SystemExit("uvx not found")
    args = _uvx_server_args(spec)[:-1]
    out: dict[str, str] = {}
    for code in (_VERSIONS, _SERVER_ORIGIN):
        ok = subprocess.run(
            [uvx, *args, "python", "-c", code],
            capture_output=True, text=True, check=False, timeout=600,
        )  # fmt: skip
        if ok.returncode != 0:
            raise SystemExit(f"cannot resolve {spec}: {ok.stderr}")
        out |= json.loads(ok.stdout.strip().splitlines()[-1])
    _expect_engine(spec, out, checked)
    commit = check_server_spec(spec)
    if commit is not None and out.get("commit") != commit:
        raise SystemExit(f"{spec} installed commit {out.get('commit')!r}, not {commit}")
    return out


def check_mafft(checked: dict[str, str]) -> None:
    exe = mafft_bin() / "mafft"
    out = subprocess.run(
        [str(exe), "--version"], capture_output=True, text=True, check=False
    )
    got = (out.stdout + out.stderr).strip()
    if checked["mafft"] not in got:
        raise SystemExit(
            f"MAFFT is {got!r}, but questions.json was checked with {checked['mafft']!r}"
        )


def arm_setup(
    arm: str, tmp: Path, spec: str = SERVER
) -> tuple[list[str], dict[str, str]]:
    """The claude flags and the environment for one arm."""
    empty = tmp.parent / f"{tmp.name}-no-mcp.json"
    empty.write_text('{"mcpServers": {}}')
    env = dict(os.environ)
    # A TMPDIR of the run's own, so nothing one run writes to a temp dir is
    # where another run looks.
    own_tmp = tmp.parent / "tmp"
    own_tmp.mkdir(exist_ok=True)
    env["TMPDIR"] = str(own_tmp)
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
                            "args": _uvx_server_args(spec),
                            "env": {"PATH": path, "TMPDIR": str(own_tmp)},
                        }
                    }
                }
            )
        )
        return [
            "--tools", "Read", "--allowedTools", "mcp__phylokit,Read",
            "--strict-mcp-config", "--mcp-config", str(cfg),
        ], env  # fmt: skip
    if arm == "read":
        env.pop("VIRTUAL_ENV", None)
        env["PATH"] = SYSTEM_PATH
        return [
            "--tools", "Read", "--allowedTools", "Read",
            "--strict-mcp-config", "--mcp-config", str(empty),
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


def _tool_outputs(events: list[dict[str, Any]]) -> list[str]:
    return [
        json.dumps(b.get("content", ""))
        for e in events
        if e.get("type") == "user"
        for b in (e.get("message") or {}).get("content", [])
        if isinstance(b, dict) and b.get("type") == "tool_result"
    ]


def leak_hits(events: list[dict[str, Any]], own_cwd: str | None) -> list[str]:
    """Every forbidden string in any tool input or output (LEAK_PATTERNS, plus
    any run directory that is not this run's own)."""
    own = Path(own_cwd).parent.name if own_cwd else None
    hits: list[str] = []
    for text in _tool_inputs(events) + _tool_outputs(events):
        for pat in LEAK_PATTERNS:
            hits += re.findall(pat, text)
        hits += [m for m in RUN_DIR.findall(text) if m != own]
    return sorted(set(hits))


def validity(arm: str, events: list[dict[str, Any]], rc: object) -> dict[str, Any]:
    """One rule for every arm; see the module docstring."""
    init = next(
        (e for e in events if e.get("type") == "system" and e.get("subtype") == "init"),
        {},
    )
    result = next((e for e in events if e.get("type") == "result"), None)
    reasons: list[str] = []
    if rc != 0:
        reasons.append(f"exit {rc}")
    if result is None:
        reasons.append("no result event")
    elif result.get("is_error") or result.get("subtype") != "success":
        reasons.append(f"result {result.get('subtype')}")
    tools = set(init.get("tools") or [])
    plain = {t for t in tools if not t.startswith("mcp__")}
    if plain != INTENDED[arm]:
        reasons.append(f"tools {sorted(plain)} != {sorted(INTENDED[arm])}")
    servers = {s.get("name"): s.get("status") for s in init.get("mcp_servers") or []}
    if arm == "mcp":
        if not any(t.startswith("mcp__phylokit__") for t in tools):
            reasons.append("no phylokit tools")
        if servers.get("phylokit") != "connected":
            reasons.append(f"server {servers.get('phylokit')}")
    elif any(t.startswith("mcp__") for t in tools) or servers:
        reasons.append("unexpected MCP server")
    calls = [
        b["name"]
        for e in events
        if e.get("type") == "assistant"
        for b in e["message"].get("content", [])
        if b.get("type") == "tool_use"
    ]
    if arm != "bare" and not calls:
        reasons.append("no tool call")
    if arm != "mcp" and any("phylokit" in s.lower() for s in _tool_inputs(events)):
        reasons.append("named phylokit")
    leaks = leak_hits(events, init.get("cwd"))
    if leaks:
        reasons.append("leak")
    return {"valid": not reasons, "invalid_reasons": reasons, "leaks": leaks}


def _home_as_tilde(path: str | None) -> str | None:
    """The stored copy of a path, with the home directory written as ~ (the
    results are committed; the leak scan reads the transcript, not this)."""
    home = str(Path.home())
    return "~" + path[len(home) :] if path and path.startswith(home + "/") else path


def check_run_dir(path: str) -> None:
    """A run's own directory must not itself name what the leak scan and the
    arm rule look for; fails loud before a run is spent."""
    bad = [p for p in LEAK_PATTERNS if re.search(p, path)]
    if "phylokit" in path.lower():
        bad.append("phylokit")
    if bad:
        raise SystemExit(f"run directory {path} names {bad}; set TMPDIR elsewhere")


def transcript(arm: str, qid: str, rep: int) -> Path:
    return RUNS / arm / f"{qid}-{rep}.jsonl"


def run_one(
    arm: str, q: dict[str, Any], rep: int, spec: str = SERVER
) -> dict[str, Any]:
    out = transcript(arm, q["id"], rep)
    out.parent.mkdir(parents=True, exist_ok=True)
    question = q["question"] + (ALT_NOTE if arm == "alt" else "")
    prompt = PROMPT.format(question=question, answer_format=q["answer_format"])
    claude = shutil.which("claude")
    if claude is None:
        raise SystemExit("claude not found on PATH")
    with tempfile.TemporaryDirectory(prefix=f"phyloqa-{arm}-") as d:
        check_run_dir(d)
        tmp = Path(d) / "cwd"
        tmp.mkdir()
        for name, text in q.get("files", {}).items():
            (tmp / name).write_text(text)
        flags, env = arm_setup(arm, tmp, spec)
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
        "server_spec": spec if arm == "mcp" else None,
        "cwd": _home_as_tilde(init.get("cwd")),
        "reply": reply,
        **score(q, reply),
        **validity(arm, events, rc),
    }


def rescore(results: Path, by_id: dict[str, dict[str, Any]]) -> None:
    """Score each stored reply again, and re-judge its validity from its
    transcript, so a truth, scorer or validity fix reaches runs made before it;
    prints every verdict or validity that changed. A missing transcript is an
    error, not a pass."""
    rows = [json.loads(ln) for ln in results.read_text().splitlines()]
    for r in rows:
        new = score(by_id[r["qid"]], r["reply"])
        path = transcript(r["arm"], r["qid"], r["rep"])
        if not path.exists():
            raise SystemExit(
                f"no transcript for {r['arm']} {r['qid']} rep{r['rep']}: {path}"
            )
        events = [
            json.loads(ln) for ln in path.read_text().splitlines() if ln.startswith("{")
        ]
        init = next((e for e in events if e.get("subtype") == "init"), {})
        new["cwd"] = _home_as_tilde(init.get("cwd"))
        new.update(validity(r["arm"], events, r["exit"]))
        if new["verdict"] != r["verdict"] or new["wrong_id"] != r["wrong_id"]:
            print(f"{r['arm']:4} {r['qid']} rep{r['rep']} {r['verdict']} -> {new['verdict']}"
                  f" ({new['value']!r})")  # fmt: skip
        if "valid" in r and new["valid"] != r["valid"]:
            print(f"{r['arm']:4} {r['qid']} rep{r['rep']} valid {r.get('valid')} -> "
                  f"{new['valid']} {new['invalid_reasons']} {new['leaks']}")  # fmt: skip
        r.update(new)
    results.write_text("".join(json.dumps(r) + "\n" for r in rows))
    invalid = [r for r in rows if not r["valid"]]
    print(f"{len(rows)} rows rescored; {len(invalid)} invalid", flush=True)
    for r in invalid:
        print(
            f"  invalid: {r['arm']} {r['qid']} rep{r['rep']} {r['invalid_reasons']} {r['leaks']}"
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="bare,alt,mcp")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", default="", help="comma-separated question ids")
    ap.add_argument("--results", default=str(HERE / "results.jsonl"))
    ap.add_argument(
        "--server-spec",
        default=SERVER,
        help=f"the mcp arm's server: {SERVER} or git+https://github.com/musharna/"
        "phylokit-mcp@<full sha> (an unreleased commit)",
    )
    ap.add_argument(
        "--versions",
        default="",
        help="write what each arm ran on to this JSON file (default: versions.json "
        "beside --results)",
    )
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
    key = json.loads((HERE / "questions.json").read_text())
    questions = key["questions"]
    # Positive control before anything is spent: the key's answers score
    # right and known-wrong ones score wrong.
    bad = oracle_check(questions)
    if bad:
        raise SystemExit("oracle check failed:\n" + "\n".join(bad))
    if a.rescore:
        rescore(Path(a.results), {q["id"]: q for q in questions})
        return
    arms = a.arms.split(",")
    unknown = set(arms) - set(ARMS)
    if unknown:
        raise SystemExit(f"unknown arms {sorted(unknown)}; known: {ARMS}")
    check_server_spec(a.server_spec)
    versions: dict[str, Any] = {"arms": arms, "model": MODEL}
    if "alt" in arms:
        versions["alt"] = check_alt(key["checked_with"])
        print("alt env:", versions["alt"], flush=True)
    if "mcp" in arms:
        versions["server_spec"] = a.server_spec
        versions["mcp"] = check_mcp(key["checked_with"], a.server_spec)
        print("server env:", versions["mcp"], flush=True)
    if {"alt", "mcp"} & set(arms):
        check_mafft(key["checked_with"])
        versions["mafft"] = key["checked_with"]["mafft"]
    claude = shutil.which("claude")
    if claude is None:
        raise SystemExit("claude not found on PATH")
    versions["claude_code"] = subprocess.run(
        [claude, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    versions["started"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    vpath = Path(a.versions or Path(a.results).with_name("versions.json"))
    if vpath.exists():
        old = json.loads(vpath.read_text())
        same = {k: old.get(k) for k in ("server_spec", "mcp", "alt", "mafft")}
        if same != {k: versions.get(k) for k in same}:
            raise SystemExit(
                f"{vpath} records another setup ({same}); a resumed sweep must run on "
                "the same one"
            )
        versions["started"] = old.get("started", versions["started"])
        versions["resumed"] = old.get("resumed", []) + [
            time.strftime("%Y-%m-%dT%H:%M:%S%z")
        ]
    vpath.write_text(json.dumps(versions, indent=2) + "\n")
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
        futures = [pool.submit(run_one, *t, a.server_spec) for t in todo]
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
