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
import subprocess
import sys

import pytest

TINY = ">a\nACGTACGTAA\n>b\nACGTACGTAA\n>c\nTCGAACGTAA\n>d\nTCGAACGTAA\n"


def test_the_real_entry_point_completes_an_mcp_handshake():
    """Launch the actual server process and speak JSON-RPC to it over stdio.

    This is the only test here that crosses the process boundary the way a
    client does. Checking `find_spec` — as an earlier version did — proves the
    module is importable, which is not the same claim: a broken entry point, a
    crash during tool registration, or anything written to stdout at import
    (which would corrupt the JSON-RPC stream) all pass an importability check
    and fail a real client.
    """
    requests = [
        {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "test", "version": "0"},
            },
        },
        {"jsonrpc": "2.0", "method": "notifications/initialized"},
        {"jsonrpc": "2.0", "id": 2, "method": "tools/list"},
    ]
    payload = "".join(json.dumps(r) + "\n" for r in requests)
    proc = subprocess.run(
        [sys.executable, "-m", "phylokit_mcp.server"],
        # check=False deliberately: the server's exit status when stdin closes
        # is not the thing under test, and raising on it would mask the far more
        # informative assertions on the JSON-RPC responses below.
        check=False,
        input=payload,
        capture_output=True,
        text=True,
        timeout=180,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip().startswith("{")]
    assert lines, f"server produced no JSON-RPC output. stderr:\n{proc.stderr[-2000:]}"

    responses = {}
    for line in lines:
        msg = json.loads(line)
        if "id" in msg:
            responses[msg["id"]] = msg

    assert 1 in responses, "no response to initialize"
    assert "result" in responses[1], responses[1]
    assert responses[1]["result"]["serverInfo"]["name"] == "phylokit-mcp"

    assert 2 in responses, "no response to tools/list"
    names = {t["name"] for t in responses[2]["result"]["tools"]}
    assert names == {
        "infer_tree",
        "select_substitution_model",
        "compare_trees",
        "simulate_alignment",
        "capabilities",
    }


def test_every_tool_is_registered_with_a_schema():
    """Registration happens at import, so a bad annotation fails here, not in prod."""
    import anyio

    from phylokit_mcp.server import mcp

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
        assert tool.inputSchema, f"{tool.name} has no input schema"
        assert tool.description, f"{tool.name} has no description"


def test_repeated_tool_calls_succeed_in_one_server_process():
    """The regression guard for engine state that survives a request.

    Four calls, engine-touching and not, interleaved. If any call after the
    first fails, the failure is in the server lifecycle rather than the maths.
    """
    import anyio

    from phylokit_mcp.server import mcp

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
    from mcp.server.fastmcp.exceptions import ToolError

    from phylokit_mcp.server import mcp

    ragged = ">a\nACGTACGTAA\n>b\nACGTACGT\n>c\nTCGAACGTAA\n>d\nTCGAACGTAA\n"

    async def run():
        with pytest.raises(ToolError, match="not all the same length"):
            await mcp.call_tool("infer_tree", {"fasta": ragged})
        # The server must still be usable afterwards.
        return await mcp.call_tool("capabilities", {"include_models": False})

    assert anyio.run(run) is not None
