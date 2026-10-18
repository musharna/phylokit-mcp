"""File-path inputs: `fasta_path` and `trees_path`.

Each refusal test also drives the legitimate path in the same test, so a
harness that refuses everything cannot pass it. The last test crosses the
process boundary: a real server process, spoken to over stdio, reads a real
file by a path relative to its own working directory.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import anyio
import pytest
from mcp.server.mcpserver.exceptions import ToolError

from phylokit_mcp import inputs
from phylokit_mcp.inputs import InputFileError, read_path, text_or_path
from phylokit_mcp.server import (
    align_sequences,
    build_server,
    compare_trees,
    infer_tree,
    select_substitution_model,
)

from .conftest import TRUE_NEWICK


def test_a_path_gives_the_same_answer_as_the_text(tmp_path, fasta_easy):
    path = tmp_path / "aln.fasta"
    path.write_text(fasta_easy)
    by_text = infer_tree(fasta_easy, model="JC", replicates=20, seed=3)
    by_path = infer_tree(fasta_path=str(path), model="JC", replicates=20, seed=3)
    assert by_path["newick"] == by_text["newick"]
    assert by_path["support"]["clades"] == by_text["support"]["clades"]

    sel_text = select_substitution_model(fasta_easy, criterion="BIC", top_n=3)
    sel_path = select_substitution_model(fasta_path=str(path), criterion="BIC", top_n=3)
    assert sel_path["ranking"] == sel_text["ranking"]


def test_trees_path_reads_tree_a_then_tree_b(tmp_path):
    other = (
        "(((A:0.05,C:0.05):0.05,(B:0.05,D:0.05):0.05):0.05,(E:0.05,F:0.05):0.05,G:0.1);"
    )
    path = tmp_path / "trees.nwk"
    path.write_text(f"{TRUE_NEWICK}\n{other}\n")
    by_path = compare_trees(trees_path=str(path))
    by_text = compare_trees(TRUE_NEWICK, other)
    assert by_path == by_text
    assert by_path["robinson_foulds"] == 4
    # A then B, not sorted or swapped: the clades only in A are TRUE_NEWICK's.
    assert by_path["only_in_a"] == compare_trees(TRUE_NEWICK, other)["only_in_a"]
    assert by_path["only_in_a"] != by_path["only_in_b"]

    path.write_text(TRUE_NEWICK)
    with pytest.raises(ValueError, match="exactly two Newick trees.*found 1"):
        compare_trees(trees_path=str(path))
    with pytest.raises(ValueError, match="not both"):
        compare_trees(newick_a=TRUE_NEWICK, trees_path=str(path))


def test_exactly_one_of_text_and_path(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(">a\nACGT\n")
    assert text_or_path(None, str(path), "fasta", "fasta_path") == ">a\nACGT\n"
    assert text_or_path(">b\nA\n", None, "fasta", "fasta_path") == ">b\nA\n"
    with pytest.raises(InputFileError, match="not both"):
        text_or_path(">b\nA\n", str(path), "fasta", "fasta_path")
    with pytest.raises(InputFileError, match=r"Give fasta .* or fasta_path"):
        text_or_path(None, None, "fasta", "fasta_path")


def test_a_path_that_is_not_a_readable_regular_file_is_refused(tmp_path, monkeypatch):
    good = tmp_path / "good.fasta"
    good.write_text(">a\nACGT\n")
    assert read_path(str(good), "fasta_path") == ">a\nACGT\n"

    with pytest.raises(InputFileError, match="is a directory"):
        read_path(str(tmp_path), "fasta_path")
    fifo = tmp_path / "fifo"
    os.mkfifo(fifo)
    # A FIFO with no writer would block the synchronous server forever if it
    # were opened; it must be refused from its stat alone.
    with pytest.raises(InputFileError, match="not a regular file"):
        read_path(str(fifo), "fasta_path")
    with pytest.raises(InputFileError, match="not a regular file"):
        read_path("/dev/zero", "fasta_path")
    bad = tmp_path / "latin1.fasta"
    bad.write_bytes(b">a\nAC\xe9GT\n")
    with pytest.raises(InputFileError, match="not UTF-8"):
        read_path(str(bad), "fasta_path")
    with pytest.raises(InputFileError, match="empty"):
        read_path("  ", "fasta_path")


def test_a_missing_relative_path_names_the_servers_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "here.fasta").write_text(">a\nACGT\n")
    assert read_path("here.fasta", "fasta_path") == ">a\nACGT\n"
    with pytest.raises(InputFileError, match=f"no such file.*{tmp_path}"):
        read_path("elsewhere.fasta", "fasta_path")


def test_the_size_cap_is_enforced_at_the_boundary(tmp_path, monkeypatch):
    monkeypatch.setattr(inputs, "MAX_INPUT_BYTES", 100)
    ok = tmp_path / "ok.fasta"
    ok.write_text(">a\n" + "A" * 97)  # 100 bytes: at the cap, accepted
    assert len(read_path(str(ok), "fasta_path")) == 100
    big = tmp_path / "big.fasta"
    big.write_text(">a\n" + "A" * 98)
    with pytest.raises(InputFileError, match="101 bytes; the cap is 100"):
        read_path(str(big), "fasta_path")


def test_a_path_refusal_reaches_the_model_through_the_server(tmp_path, fasta_easy):
    path = tmp_path / "aln.fasta"
    path.write_text(fasta_easy)
    mcp = build_server()

    async def run():
        ok = await mcp.call_tool(
            "infer_tree",
            {"fasta_path": str(path), "model": "JC", "replicates": 20},
        )
        with pytest.raises(ToolError, match="no such file"):
            await mcp.call_tool("infer_tree", {"fasta_path": str(tmp_path / "nope")})
        return ok

    assert anyio.run(run) is not None


def test_align_sequences_reads_a_path(tmp_path):
    seqs = ">a\nACGTACGTTTGA\n>b\nACGTACGTGA\n>c\nACGAACGTTTGA\n>d\nACGTTCGTTTGA\n"
    path = tmp_path / "seqs.fasta"
    path.write_text(seqs)
    assert (
        align_sequences(fasta_path=str(path))["fasta"] == align_sequences(seqs)["fasta"]
    )


def _rpc(proc, msg):
    proc.stdin.write(json.dumps(msg) + "\n")
    proc.stdin.flush()
    if "id" not in msg:
        return None
    while True:
        line = proc.stdout.readline()
        if not line:
            raise AssertionError("server closed stdout")
        resp = json.loads(line)
        if resp.get("id") == msg["id"]:
            return resp


def test_a_real_server_process_reads_a_relative_path_from_its_cwd(tmp_path):
    """Over stdio, as a client launches it: the tool reads `trees.nwk` from the
    server process's own working directory."""
    other = (
        "(((A:0.05,C:0.05):0.05,(B:0.05,D:0.05):0.05):0.05,(E:0.05,F:0.05):0.05,G:0.1);"
    )
    (tmp_path / "trees.nwk").write_text(f"{TRUE_NEWICK}\n{other}\n")
    with tempfile.TemporaryFile(mode="w+") as err:
        proc = subprocess.Popen(
            [sys.executable, "-m", "phylokit_mcp.server"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=err,
            text=True,
            bufsize=1,
            cwd=tmp_path,
        )
        try:
            _rpc(proc, {
                "jsonrpc": "2.0", "id": 1, "method": "initialize",
                "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                           "clientInfo": {"name": "test", "version": "0"}},
            })  # fmt: skip
            _rpc(proc, {"jsonrpc": "2.0", "method": "notifications/initialized"})
            good = _rpc(proc, {
                "jsonrpc": "2.0", "id": 2, "method": "tools/call",
                "params": {"name": "compare_trees", "arguments": {"trees_path": "trees.nwk"}},
            })  # fmt: skip
            bad = _rpc(proc, {
                "jsonrpc": "2.0", "id": 3, "method": "tools/call",
                "params": {"name": "compare_trees", "arguments": {"trees_path": "nope.nwk"}},
            })  # fmt: skip
        finally:
            proc.stdin.close()
            proc.wait(timeout=30)
    assert good["result"].get("isError") is not True, good
    assert good["result"]["structuredContent"]["robinson_foulds"] == 4
    assert bad["result"]["isError"] is True
    assert "no such file" in bad["result"]["content"][0]["text"]
    assert str(tmp_path) in bad["result"]["content"][0]["text"]
