"""Concurrent tool calls in one server process.

MCPServer runs a synchronous tool in a worker thread (`anyio.to_thread`), so
two calls a client sends in one turn run at the same time. piqtree's engine
wrapper changes process-global state for the length of a call: it chdirs into
a temporary directory (deleted afterwards) and points file descriptors 1 and 2
at /dev/null. Two of those interleaved leave the process in a deleted
directory, after which every call fails. In the 2026-10-18 benchmark re-run,
the client sent three `infer_tree` calls in one turn and the server failed
every call after them.
"""

from __future__ import annotations

import os

import anyio

from phylokit_mcp.server import build_server

from .conftest import TRUE_NEWICK


def test_parallel_engine_calls_all_succeed_and_leave_the_cwd_alone(
    tmp_path, monkeypatch, fasta_easy
):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "aln.fasta").write_text(fasta_easy)
    mcp = build_server()
    results: dict[int, object] = {}

    async def call(i: int, name: str, args: dict) -> None:
        results[i] = await mcp.call_tool(name, args)

    calls = [
        ("infer_tree", {"fasta_path": "aln.fasta", "model": m, "replicates": 20})
        for m in ("JC", "K2P", "HKY", "F81")
    ] + [
        ("select_substitution_model", {"fasta_path": "aln.fasta", "top_n": 1}),
        ("simulate_alignment", {"newick": TRUE_NEWICK, "length": 100}),
    ]

    async def run() -> None:
        async with anyio.create_task_group() as tg:
            for i, (name, args) in enumerate(calls):
                tg.start_soon(call, i, name, args)

    anyio.run(run)
    assert sorted(results) == list(range(len(calls)))
    assert os.getcwd() == str(tmp_path)
    # And the process can still do work afterwards.
    assert anyio.run(
        mcp.call_tool,
        "compare_trees",
        {"newick_a": TRUE_NEWICK, "newick_b": TRUE_NEWICK},
    )
