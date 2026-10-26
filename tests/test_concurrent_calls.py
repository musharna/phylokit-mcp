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


def test_queued_in_process_calls_survive_the_heartbeat(monkeypatch, fasta_easy):
    """An in-process call has a Context but no request: no heartbeat, no error.

    `mcp.call_tool` (this file's own path, and any embedding client) injects a
    Context with no request behind it, and `report_progress` then raises. With
    the heartbeat at 15 s, a call queued behind the lock for longer died with
    "Context is not available outside of a request" while its work ran on, and
    the orphaned engine call broke the fds of every later test (2026-10-25
    full-suite run). A tiny interval makes every call reach the heartbeat.
    """
    from phylokit_mcp import server

    monkeypatch.setattr(server, "HEARTBEAT_S", 0.01)
    mcp = build_server()
    results: dict[int, object] = {}

    async def call(i: int) -> None:
        results[i] = await mcp.call_tool(
            "infer_tree", {"fasta": fasta_easy, "model": "JC", "replicates": 20}
        )

    async def run() -> None:
        async with anyio.create_task_group() as tg:
            for i in range(3):
                tg.start_soon(call, i)

    anyio.run(run)
    assert sorted(results) == [0, 1, 2]
    assert not server._ENGINE_LOCK.locked()


def test_a_failing_heartbeat_fails_the_call_only_after_the_work_is_done(monkeypatch):
    """Fail loud, but never orphan the work: the lock is free when the error lands.

    The positive control is in the same test: with a context that sends, the
    same call returns its value.
    """
    import threading
    import time

    from phylokit_mcp import server

    monkeypatch.setattr(server, "HEARTBEAT_S", 0.01)
    finished = threading.Event()

    def work() -> str:
        with server._ENGINE_LOCK:
            time.sleep(0.2)
            finished.set()
            return "done"

    class Broken:
        @property
        def request_context(self):
            return object()

        async def report_progress(self, *a, **k):
            raise ConnectionResetError("client went away")

    class Fine(Broken):
        sent = 0

        async def report_progress(self, *a, **k):
            Fine.sent += 1

    async def drive(ctx):
        return await anyio.to_thread.run_sync(
            lambda: server._with_heartbeat(work, ctx, "work")
        )

    try:
        anyio.run(drive, Broken())
    except RuntimeError as exc:
        assert "progress notification failed" in str(exc)
        assert isinstance(exc.__cause__, ConnectionResetError)
    else:
        raise AssertionError("a failing heartbeat must fail the call")
    assert finished.is_set() and not server._ENGINE_LOCK.locked()

    finished.clear()
    assert anyio.run(drive, Fine()) == "done"
    assert Fine.sent >= 1
