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
import selectors
import subprocess
import sys
import tempfile
import time

import pytest

TINY = ">a\nACGTACGTAA\n>b\nACGTACGTAA\n>c\nTCGAACGTAA\n>d\nTCGAACGTAA\n"

HANDSHAKE_TIMEOUT = 180


def _read_until(stdout, wanted, stderr_file, timeout=HANDSHAKE_TIMEOUT):
    """Read JSON-RPC lines until every id in `wanted` has answered, or time out.

    Returns as soon as the wanted ids arrive, so the caller can keep stdin open
    until it actually has what it asked for.
    """
    sel = selectors.DefaultSelector()
    sel.register(stdout, selectors.EVENT_READ)
    responses: dict = {}
    deadline = time.monotonic() + timeout
    try:
        while not wanted <= responses.keys():
            remaining = deadline - time.monotonic()
            if remaining <= 0 or not sel.select(timeout=remaining):
                break
            line = stdout.readline()
            if not line:  # server closed stdout
                break
            line = line.strip()
            if not line.startswith("{"):
                continue
            msg = json.loads(line)
            if "id" in msg:
                responses[msg["id"]] = msg
    finally:
        sel.close()
    if not wanted <= responses.keys():
        stderr_file.seek(0)
        missing = sorted(wanted - responses.keys())
        raise AssertionError(
            f"no response to id(s) {missing} within {timeout}s. "
            f"got ids {sorted(responses)}. server stderr:\n{stderr_file.read()[-2000:]}"
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
