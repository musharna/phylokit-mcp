"""Drive the server as a real MCP process over stdio.

Every other test in this suite calls the tool functions directly, which imports
the engine into the pytest process and can mask state that only misbehaves under
the real server lifecycle. A sibling project shipped a P0 that its whole
in-process suite passed: the second tool call in a fresh server died, because a
context-local was populated inside whichever request touched the engine first.

So this launches the actual entry point and makes SEVERAL calls in one process.
The second call is the point — a suite that only ever makes one cannot see this
class of bug.
"""

from __future__ import annotations

import json
import os
import select
import subprocess
import sys
import tempfile
import time
import weakref

import pytest

TINY = ">a\nACGTACGTAA\n>b\nACGTACGTAA\n>c\nTCGAACGTAA\n>d\nTCGAACGTAA\n"

HANDSHAKE_TIMEOUT = 180
# The progress test bounds SILENCE, which is what a client's idle timeout acts
# on, not the calls' total run time, which grows with a busy host: two queued
# 100-replicate bootstraps took 292 s on a laptop at load 50 and failed a 180 s
# total. STREAM_SILENCE_S is many heartbeats long; RUN_CAP_S only stops a hang.
STREAM_SILENCE_S = 30
RUN_CAP_S = 1800


class _Lines:
    """Lines from a pipe with a timeout, read from its fd so no line can hide.

    `select()` on the fd next to `readline()` on the text wrapper loses lines:
    one read() can pull two messages into Python's buffer, and select then
    reports the drained pipe idle while the second waits there. If that second
    message is the last one, nothing ever wakes the select and a finished call
    reads as a silent one. Here the only buffer is this one.

    A line the peer cut off by closing is not a line: `next` still reports EOF,
    and the fragment waits in `cut_off` for the caller's failure message.
    """

    def __init__(self, pipe):
        self.fd = pipe.fileno()
        self.buf = b""
        self.cut_off = ""

    def next(self, timeout: float) -> str | None:
        """The next non-blank line; None if none completes within `timeout`
        seconds; "" once the other end has closed."""
        deadline = time.monotonic() + timeout
        while True:
            line, nl, rest = self.buf.partition(b"\n")
            if nl:
                self.buf = rest
                if line.strip():
                    return line.decode().strip()
                continue
            remaining = deadline - time.monotonic()
            if remaining <= 0 or not select.select([self.fd], [], [], remaining)[0]:
                return None
            chunk = os.read(self.fd, 65536)
            if not chunk:
                self.cut_off += self.buf.decode(errors="replace").strip()
                self.buf = b""
                return ""
            self.buf += chunk


# One reader per pipe for its whole life: a second reader would not see the
# bytes the first one has already taken off the fd.
_READERS: weakref.WeakKeyDictionary = weakref.WeakKeyDictionary()


def _lines(pipe) -> _Lines:
    if pipe not in _READERS:
        _READERS[pipe] = _Lines(pipe)
    return _READERS[pipe]


def _read_until(stdout, wanted, stderr_file, timeout=HANDSHAKE_TIMEOUT):
    """Read JSON-RPC lines until every id in `wanted` has answered, or time out.

    Returns as soon as the wanted ids arrive, so the caller can keep stdin open
    until it actually has what it asked for.
    """
    lines = _lines(stdout)
    responses: dict = {}
    deadline = time.monotonic() + timeout
    while not wanted <= responses.keys():
        line = lines.next(deadline - time.monotonic())
        if not line:  # timed out (None) or the server closed stdout ("")
            break
        if not line.startswith("{"):
            continue
        msg = json.loads(line)
        if "id" in msg:
            responses[msg["id"]] = msg
    if not wanted <= responses.keys():
        stderr_file.seek(0)
        missing = sorted(wanted - responses.keys())
        cut = f"cut off at EOF: {lines.cut_off}. " if lines.cut_off else ""
        raise AssertionError(
            f"no response to id(s) {missing} within {timeout}s. "
            f"got ids {sorted(responses)}. {cut}"
            f"server stderr:\n{stderr_file.read()[-2000:]}"
        )
    return responses


def test_the_real_entry_point_completes_an_mcp_handshake():
    """Launch the actual server process and speak JSON-RPC to it over stdio.

    This is the only test here that crosses the process boundary the way a
    client does. Checking `find_spec` — as an earlier version did — proves the
    module is importable, which is not the same claim: a broken entry point, a
    crash during tool registration, or anything written to stdout at import
    (which would corrupt the JSON-RPC stream) all pass an importability check
    and fail a real client.

    Driven as a request/response conversation rather than a batch dump. An
    earlier version wrote all three messages, closed stdin immediately, and
    parsed stdout after the process exited — which made answering `tools/list`
    a race against the server's own EOF-triggered shutdown. It failed roughly
    2 runs in 10 (measured), and CI caught it on a README-only commit. Nothing
    about a *timeout* would have fixed that: the server was not slow, it was
    already shutting down. stdin now stays open until the wanted ids arrive.
    """
    init = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "0"},
        },
    }
    # stderr goes to a file, not a pipe nobody drains: a full stderr pipe would
    # block the server mid-handshake and read as a protocol failure.
    with tempfile.TemporaryFile(mode="w+") as stderr_file:
        proc = subprocess.Popen(
            [sys.executable, "-m", "phylokit_mcp.server"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=stderr_file,
            text=True,
            bufsize=1,
        )
        # Popen with stdin/stdout=PIPE always supplies both; naming them locally
        # is what lets a type checker see that.
        stdin, stdout = proc.stdin, proc.stdout
        assert stdin is not None and stdout is not None
        try:
            stdin.write(json.dumps(init) + "\n")
            stdin.flush()
            responses = _read_until(stdout, {1}, stderr_file)

            assert "result" in responses[1], responses[1]
            assert responses[1]["result"]["serverInfo"]["name"] == "phylokit-mcp"

            for msg in (
                {"jsonrpc": "2.0", "method": "notifications/initialized"},
                {"jsonrpc": "2.0", "id": 2, "method": "tools/list"},
            ):
                stdin.write(json.dumps(msg) + "\n")
            stdin.flush()
            responses |= _read_until(stdout, {2}, stderr_file)

            names = {t["name"] for t in responses[2]["result"]["tools"]}
            assert names == {
                "infer_tree",
                "select_substitution_model",
                "compare_trees",
                "simulate_alignment",
                "align_sequences",
                "map_positions",
                "capabilities",
            }
        finally:
            # Only now is EOF safe — every response we needed is already in hand.
            if not stdin.closed:
                stdin.close()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=30)


def test_every_tool_is_registered_with_a_schema():
    """Registration happens at import, so a bad annotation fails here, not in prod."""
    import anyio

    from phylokit_mcp.server import build_server

    mcp = build_server()

    tools = anyio.run(mcp.list_tools)
    names = {t.name for t in tools}
    assert names == {
        "infer_tree",
        "select_substitution_model",
        "compare_trees",
        "simulate_alignment",
        "align_sequences",
        "map_positions",
        "capabilities",
    }
    for tool in tools:
        # mcp 2.x renamed the Tool model's fields from camelCase to snake_case
        # (inputSchema -> input_schema). The server needed no change for this;
        # only assertions that read the model do.
        assert tool.input_schema, f"{tool.name} has no input schema"
        assert tool.description, f"{tool.name} has no description"


def test_repeated_tool_calls_succeed_in_one_server_process():
    """The regression guard for engine state that survives a request.

    Four calls, engine-touching and not, interleaved. If any call after the
    first fails, the failure is in the server lifecycle rather than the maths.
    """
    import anyio

    from phylokit_mcp.server import build_server

    mcp = build_server()

    async def run():
        results = []
        for _ in range(2):
            results.append(
                await mcp.call_tool("capabilities", {"include_models": False})
            )
            results.append(
                await mcp.call_tool(
                    "infer_tree",
                    {"fasta": TINY, "model": "JC", "replicates": 20, "seed": 1},
                )
            )
        return results

    out = anyio.run(run)
    assert len(out) == 4
    for i, result in enumerate(out):
        assert result is not None, f"call {i + 1} returned nothing"


def test_a_ragged_alignment_surfaces_as_a_tool_error_not_a_crash():
    """A bad input must not take the server down, and must say what was wrong.

    The assertion matches the MESSAGE, not merely 'an exception was raised'.
    A bare `pytest.raises(Exception)` passes on almost any breakage — including
    the server failing for a completely different reason — so it cannot
    distinguish a good error from a bad one. The input here is four sequences of
    unequal length, which reaches the ragged-alignment check specifically; an
    earlier version passed only two sequences and was in fact exercising the
    taxon-count check while claiming to test raggedness.
    """
    import anyio
    from mcp.server.mcpserver.exceptions import ToolError

    from phylokit_mcp.server import build_server

    mcp = build_server()

    ragged = ">a\nACGTACGTAA\n>b\nACGTACGT\n>c\nTCGAACGTAA\n>d\nTCGAACGTAA\n"

    async def run():
        with pytest.raises(ToolError, match="not all the same length"):
            await mcp.call_tool("infer_tree", {"fasta": ragged})
        # The server must still be usable afterwards.
        return await mcp.call_tool("capabilities", {"include_models": False})

    assert anyio.run(run) is not None


def test_build_server_returns_a_fresh_instance_each_call():
    """The property the factory exists for.

    Without this, `build_server` could quietly become a cached singleton again --
    every other test would still pass, because they only ever need *a* server.
    """
    import anyio

    from phylokit_mcp.server import build_server

    a, b = build_server(), build_server()
    assert a is not b, "build_server returned the same object twice"

    # Positive control: `a is not b` is also satisfied by two EMPTY servers, so a
    # constructor that registered nothing would pass the assertion above. Both
    # instances must actually carry the full tool surface.
    expected = {
        "infer_tree",
        "select_substitution_model",
        "compare_trees",
        "simulate_alignment",
        "align_sequences",
        "map_positions",
        "capabilities",
    }
    for srv in (a, b):
        assert {t.name for t in anyio.run(srv.list_tools)} == expected


def test_importing_the_module_does_not_construct_a_server():
    """Importing must have no side effect of building a server.

    This is the regression guard for the module-level `mcp = MCPServer(...)` this
    package used to carry: a singleton built at import time is shared by every
    test that touches the module, so state leaks between them in an order nothing
    declares.
    """
    import phylokit_mcp.server as server_module

    assert not hasattr(server_module, "mcp"), (
        "phylokit_mcp.server has a module-level `mcp` again -- importing the "
        "module now constructs a server, which is what build_server() replaced"
    )


def test_long_calls_send_progress_so_a_client_idle_timeout_does_not_fire(fasta_easy):
    """Heartbeats reach a real client while a call runs AND while it waits.

    Claude Code aborts a tool call that sends "no response or progress" for its
    idle timeout. On the field benchmark every 89-taxon `infer_tree` was silent
    for the whole bootstrap and was aborted at 1822 s, ten runs in ten. The
    server now sends `notifications/progress` every PHYLOKIT_HEARTBEAT_S while a
    tool runs -- and while it queues behind `_ENGINE_LOCK`, which is silent too.

    Real process, real stdio: piqtree points fds 1 and 2 at /dev/null during
    every engine call, so a notification written to fd 1 then would be lost.
    This test is what shows the wire is not fd 1 (mcp's stdio transport serves
    it from a private duplicate) -- a fixture could not.
    """
    env = {**os.environ, "PHYLOKIT_HEARTBEAT_S": "0.25"}
    calls = {
        3: {"fasta": fasta_easy, "model": "JC", "replicates": 100, "seed": 1},
        4: {"fasta": fasta_easy, "model": "JC", "replicates": 100, "seed": 2},
    }
    with tempfile.TemporaryFile(mode="w+") as stderr_file:
        proc = subprocess.Popen(
            [sys.executable, "-m", "phylokit_mcp.server"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=stderr_file,
            text=True,
            bufsize=1,
            env=env,
        )
        stdin, stdout = proc.stdin, proc.stdout
        assert stdin is not None and stdout is not None
        try:
            init = {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": {"name": "test", "version": "0"},
                },
            }
            stdin.write(json.dumps(init) + "\n")
            stdin.flush()
            _read_until(stdout, {1}, stderr_file)
            stdin.write(
                json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"})
                + "\n"
            )
            # Sent together: the second queues behind the engine lock.
            for rid, args in calls.items():
                msg = {
                    "jsonrpc": "2.0",
                    "id": rid,
                    "method": "tools/call",
                    "params": {
                        "name": "infer_tree",
                        "arguments": args,
                        "_meta": {"progressToken": f"tok{rid}"},
                    },
                }
                stdin.write(json.dumps(msg) + "\n")
            stdin.flush()

            progress: dict[str, list[float]] = {"tok3": [], "tok4": []}
            done: dict[int, dict] = {}
            lines = _lines(stdout)
            cap = time.monotonic() + RUN_CAP_S
            silent_for: float | None = None
            last = time.monotonic()
            while set(done) != set(calls) and time.monotonic() < cap:
                line = lines.next(min(STREAM_SILENCE_S, cap - time.monotonic()))
                if line is None:
                    # Near the cap the wait is cut short; that is the hang
                    # guard firing, reported below as unanswered calls.
                    if time.monotonic() < cap:
                        silent_for = time.monotonic() - last
                    break
                if not line:  # the server closed stdout
                    break
                last = time.monotonic()
                if not line.startswith("{"):
                    continue
                msg = json.loads(line)
                if msg.get("method") == "notifications/progress":
                    p = msg["params"]
                    progress[p["progressToken"]].append(p["progress"])
                elif msg.get("id") in calls:
                    done[msg["id"]] = msg
            stderr_file.seek(0)
            assert silent_for is None, (
                f"the stream was silent for {silent_for:.1f} s with calls "
                f"{sorted(set(calls) - set(done))} still running: "
                + stderr_file.read()[-2000:]
            )
            assert set(done) == set(calls), (
                f"cut off at EOF: {lines.cut_off!r}. " if lines.cut_off else ""
            ) + stderr_file.read()[-2000:]
        finally:
            if not stdin.closed:
                stdin.close()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=30)

    # Positive control: both calls still return their trees.
    for rid in calls:
        res = done[rid]["result"]
        assert not res.get("isError"), res
        assert res["structuredContent"]["support"]["replicates"] == 100
    # Each call heartbeat at least twice, and the values only go up (the spec
    # requires progress to increase with each notification).
    for tok, values in progress.items():
        assert len(values) >= 2, (tok, values)
        assert values == sorted(values) and len(set(values)) == len(values), values


def test_a_line_already_taken_off_the_pipe_is_returned_at_once():
    """Two messages in one write. `select()` then `readline()` -- the pattern
    `_Lines` replaced -- reported the second one as 2 s of silence; a last
    response stranded that way failed the progress test with the server done."""
    r, w = os.pipe()
    with os.fdopen(r) as pipe:
        lines = _lines(pipe)
        os.write(w, b'{"a": 1}\n\n{"b": 2}\n{"c"')
        assert lines.next(5) == '{"a": 1}'
        t = time.monotonic()
        assert lines.next(5) == '{"b": 2}'
        assert time.monotonic() - t < 1
        # A real silence still reads as one, and a part-line is not a line.
        assert lines.next(0.2) is None
        os.write(w, b": 3}\n")
        assert lines.next(5) == '{"c": 3}'
        assert _lines(pipe) is lines
        # A server that dies mid-message: EOF, with the fragment kept aside.
        assert lines.cut_off == ""
        os.write(w, b'{"d": ')
        os.close(w)
        assert lines.next(5) == ""
        assert lines.cut_off == '{"d":'


def test_a_server_that_dies_mid_message_is_reported_with_what_it_sent():
    """The caller's side of a crash mid-response: `_read_until` names the
    unanswered id, the cut-off text and the server's stderr, instead of
    stopping on the fragment. Id 0 arrives whole: the positive control."""
    script = (
        "import os, sys; sys.stderr.write('boom\\n'); sys.stderr.flush(); "
        'os.write(1, b\'{"id": 0}\\n{"id": 1, "res\')'
    )
    with tempfile.TemporaryFile(mode="w+") as stderr_file:
        proc = subprocess.Popen(
            [sys.executable, "-c", script],
            stdout=subprocess.PIPE,
            stderr=stderr_file,
            text=True,
            bufsize=1,
        )
        try:
            with pytest.raises(AssertionError) as err:
                _read_until(proc.stdout, {0, 1}, stderr_file, timeout=30)
        finally:
            proc.wait(timeout=30)
    assert str(err.value) == (
        "no response to id(s) [1] within 30s. got ids [0]. "
        'cut off at EOF: {"id": 1, "res. server stderr:\nboom\n'
    )
