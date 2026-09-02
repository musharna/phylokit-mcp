"""Every tool must reach the model through the refusal boundary.

Since mcp 2.1 (python-sdk #3314) MCPServer masks any exception that is not a
ToolError: the caller gets `Error executing tool <name>` and the reason stays in
the server log. `server._surfaces_refusals` converts this server's anticipated
refusals at that boundary so their text survives, because the text IS the
product -- a ragged alignment or an unknown criterion each says what to do
instead.

The wrapper is applied per tool, by hand. `test_mcp_protocol` asserts the text of
ONE refusal from ONE tool, which pins the instance: it passes unchanged on the
day a sixth tool is registered without the wrapper. This pins the class, by
reading the registry the server actually exposes -- so the failure lands here
rather than in a transcript where someone reads `Error executing tool` and has
nothing to act on.

Reading the built server rather than parsing the source also keeps the guard
true for either registration shape: the call-site wrap used here, and the
stacked decorator breedsim-mcp uses.
"""

from __future__ import annotations

import pytest
from mcp.server.mcpserver.exceptions import ToolError

from phylokit_mcp import server


def _boundary_code():
    """The code object of the wrapper `_surfaces_refusals` returns.

    Every wrapper it builds shares one code object, so identity against it tests
    exactly "this callable came out of that boundary". `hasattr(fn, "__wrapped__")`
    would also be satisfied by any unrelated functools.wraps decorator, and
    matching on `__name__` by nothing at all -- functools.wraps copies the
    wrapped function's name onto the wrapper.
    """
    return server._surfaces_refusals(lambda: None).__code__


def _registered_tools():
    tools = server.build_server()._tool_manager.list_tools()
    assert tools, "no tools registered -- this guard would otherwise pass vacuously"
    return tools


def test_every_registered_tool_passes_through_the_refusal_boundary():
    expected = _boundary_code()
    unguarded = sorted(
        t.name
        for t in _registered_tools()
        if getattr(t.fn, "__code__", None) is not expected
    )
    assert not unguarded, (
        f"{unguarded} are registered without _surfaces_refusals, so a refusal "
        "from them reaches the model as `Error executing tool <name>` with the "
        "reason dropped. Wrap them where they are registered in build_server()."
    )


def test_the_boundary_converts_a_refusal_and_leaves_a_real_bug_masked():
    """The negative assertion needs the positive one beside it.

    A boundary that converted every exception would satisfy the refusal half
    while destroying the SDK's crash signal, so both directions are asserted in
    one test rather than trusting a broken harness to look like a pass.
    """

    @server._surfaces_refusals
    def refuses():
        raise ValueError("alignment rows are not all the same length")

    @server._surfaces_refusals
    def crashes():
        raise TypeError("this is a bug, not a refusal")

    with pytest.raises(ToolError, match="not all the same length"):
        refuses()

    with pytest.raises(TypeError):
        crashes()
