"""`--root` / PHYLOKIT_ROOTS: an optional allow-list for every path a tool is given.

Unset, nothing changes: the server opens any path its user can. Set, a
`fasta_path` or `trees_path` outside the roots is refused before it is opened
or stat'ed, and `out_path` creates nothing outside them.

Every refusal test drives the legitimate path in the same test, and checks the
effect, not only the error: an outside file is readable (so the refusal is the
policy, not a missing file) and nothing appears in the outside directory. The
last test starts real server processes, one with `--root` and one without.
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
from phylokit_mcp.inputs import (
    InputFileError,
    PathOutsideRootsError,
    read_path,
    write_new_file,
)
from phylokit_mcp.roots import (
    ENV,
    CannotLocateDescriptorError,
    configured_roots,
    set_roots,
)
from phylokit_mcp.server import build_server

from .conftest import TRUE_NEWICK

OTHER = "(((A:0.05,C:0.05):0.05,(B:0.05,D:0.05):0.05):0.05,(E:0.05,F:0.05):0.05,G:0.1);"
OUTSIDE = "outside the configured roots"


@pytest.fixture
def dirs(tmp_path, monkeypatch):
    """A root and a sibling outside it, each holding a readable file."""
    monkeypatch.delenv(ENV, raising=False)
    root, outside = tmp_path / "root", tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (root / "in.fa").write_text(">a\nACGT\n")
    (outside / "secret.fa").write_text(">s\nTTTT\n")
    yield root, outside
    set_roots(None)


def test_with_no_roots_every_path_is_open_as_before(dirs):
    _root, outside = dirs
    assert configured_roots() is None
    assert read_path(str(outside / "secret.fa"), "fasta_path") == ">s\nTTTT\n"
    written = write_new_file(str(outside / "new.fa"), ">n\nA\n", "out_path")
    assert written == str(outside / "new.fa")
    assert (outside / "new.fa").read_text() == ">n\nA\n"


def test_reads_outside_the_roots_are_refused_by_every_route(dirs, monkeypatch):
    root, outside = dirs
    (root / "link.fa").symlink_to(outside / "secret.fa")
    (root / "linkdir").symlink_to(outside, target_is_directory=True)
    (root / "back.fa").symlink_to(root / "in.fa")
    set_roots([str(root)])

    # Inside, including a path that walks out and back in, and a symlink that
    # stays inside.
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"
    assert read_path(str(root / ".." / "root" / "in.fa"), "fasta_path") == ">a\nACGT\n"
    assert read_path(str(root / "back.fa"), "fasta_path") == ">a\nACGT\n"

    for route in (
        outside / "secret.fa",  # absolute
        root / ".." / "outside" / "secret.fa",  # dot-dot
        root / "link.fa",  # symlinked file
        root / "linkdir" / "secret.fa",  # symlinked directory
    ):
        with pytest.raises(PathOutsideRootsError, match=OUTSIDE):
            read_path(str(route), "fasta_path")

    # Relative to a working directory outside the roots.
    monkeypatch.chdir(outside)
    with pytest.raises(PathOutsideRootsError, match=OUTSIDE):
        read_path("secret.fa", "fasta_path")
    monkeypatch.chdir(root)
    assert read_path("in.fa", "fasta_path") == ">a\nACGT\n"


def test_a_refusal_does_not_say_whether_the_outside_file_exists(dirs):
    root, outside = dirs
    set_roots([str(root)])
    messages = []
    for name in ("secret.fa", "no-such-file.fa"):
        with pytest.raises(PathOutsideRootsError) as exc:
            read_path(str(outside / name), "fasta_path")
        messages.append(str(exc.value).replace(name, "<name>"))
    assert messages[0] == messages[1]
    # Inside the roots a missing file is still named as missing.
    with pytest.raises(InputFileError, match="no such file"):
        read_path(str(root / "no-such-file.fa"), "fasta_path")


def test_inside_the_roots_the_ordinary_refusals_still_hold(dirs, monkeypatch):
    """Under roots the file is opened before it is examined, so each check the
    unconfined read makes from `stat` is made again here from the descriptor."""
    root, _ = dirs
    set_roots([str(root)])
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"

    (root / "adir").mkdir()
    with pytest.raises(InputFileError, match="is a directory"):
        read_path(str(root / "adir"), "fasta_path")
    os.mkfifo(root / "fifo")
    # A FIFO with no writer: opened without blocking, then refused.
    with pytest.raises(InputFileError, match="not a regular file"):
        read_path(str(root / "fifo"), "fasta_path")
    (root / "latin1.fa").write_bytes(b">a\nAC\xe9GT\n")
    with pytest.raises(InputFileError, match="not UTF-8"):
        read_path(str(root / "latin1.fa"), "fasta_path")

    monkeypatch.setattr(inputs, "MAX_INPUT_BYTES", 100)
    (root / "ok.fa").write_text(">a\n" + "A" * 97)  # 100 bytes: at the cap
    assert len(read_path(str(root / "ok.fa"), "fasta_path")) == 100
    (root / "big.fa").write_text(">a\n" + "A" * 98)
    with pytest.raises(InputFileError, match="101 bytes; the cap is 100"):
        read_path(str(root / "big.fa"), "fasta_path")


def test_writes_are_confined_and_nothing_is_created_outside(dirs):
    root, outside = dirs
    (root / "linkdir").symlink_to(outside, target_is_directory=True)
    (root / "dangling.fa").symlink_to(outside / "planted.fa")
    set_roots([str(root)])

    written = write_new_file(str(root / "aln.fa"), ">n\nA\n", "out_path")
    assert written == str(root / "aln.fa")
    assert (root / "aln.fa").read_text() == ">n\nA\n"
    with pytest.raises(InputFileError, match="already exists"):
        write_new_file(str(root / "aln.fa"), ">n\nC\n", "out_path")
    assert (root / "aln.fa").read_text() == ">n\nA\n"

    for route in (
        outside / "new.fa",  # absolute
        root / ".." / "outside" / "new.fa",  # dot-dot
        root / "linkdir" / "new.fa",  # symlinked directory
        root / "dangling.fa",  # a symlink at the file's own name
    ):
        with pytest.raises(PathOutsideRootsError, match=OUTSIDE):
            write_new_file(str(route), ">n\nA\n", "out_path")
    assert sorted(p.name for p in outside.iterdir()) == ["secret.fa"]


def test_a_path_swapped_after_the_name_check_is_caught_at_the_descriptor(
    dirs, monkeypatch
):
    """The name check and the open are two system calls. Here the path is
    re-pointed outside the roots between them, as a second process could: the
    file that was opened, or the directory the new file would land in, is
    what gets judged."""
    root, outside = dirs
    (root / "sub").mkdir()
    (root / "moving.fa").write_text(">a\nACGT\n")
    set_roots([str(root)])
    real_check = inputs._require_name_inside
    swaps: list = []

    def check_then_swap(p, shown, param, roots):
        real_check(p, shown, param, roots)
        while swaps:
            swaps.pop()()

    monkeypatch.setattr(inputs, "_require_name_inside", check_then_swap)

    # Control: with no swap, both go through the patched check untouched.
    assert read_path(str(root / "moving.fa"), "fasta_path") == ">a\nACGT\n"
    assert write_new_file(str(root / "sub" / "a.fa"), ">n\nA\n", "out_path")

    def repoint_file():
        (root / "moving.fa").unlink()
        (root / "moving.fa").symlink_to(outside / "secret.fa")

    swaps.append(repoint_file)
    with pytest.raises(PathOutsideRootsError, match="Nothing was read"):
        read_path(str(root / "moving.fa"), "fasta_path")

    # Re-pointed at an outside DIRECTORY: still the roots refusal, not "is a
    # directory". Nothing about the outside path is looked at, its kind included.
    def repoint_at_a_directory():
        (root / "moving.fa").unlink()
        (root / "moving.fa").symlink_to(outside, target_is_directory=True)

    (root / "moving.fa").unlink()
    (root / "moving.fa").write_text(">a\nACGT\n")  # inside again for the name check
    swaps.append(repoint_at_a_directory)
    with pytest.raises(PathOutsideRootsError, match="Nothing was read") as exc:
        read_path(str(root / "moving.fa"), "fasta_path")
    assert "is a directory" not in str(exc.value)

    def repoint_directory():
        (root / "sub").rename(root / "sub-moved")
        (root / "sub").symlink_to(outside, target_is_directory=True)

    swaps.append(repoint_directory)
    with pytest.raises(PathOutsideRootsError, match="Nothing was written"):
        write_new_file(str(root / "sub" / "b.fa"), ">n\nA\n", "out_path")
    assert sorted(p.name for p in outside.iterdir()) == ["secret.fa"]


def test_where_a_descriptor_cannot_be_located_roots_refuse_rather_than_skip(
    dirs, monkeypatch
):
    """The binding check asks the kernel where an open descriptor lives. A
    platform that cannot answer gets a refusal, never an unchecked read or
    write; with no roots set the same platform is unaffected."""
    root, _outside = dirs
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"  # no roots

    def cannot(fd):
        raise CannotLocateDescriptorError("this platform cannot report it")

    set_roots([str(root)])
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"
    monkeypatch.setattr(inputs, "fd_path", cannot)
    with pytest.raises(InputFileError, match="cannot report it"):
        read_path(str(root / "in.fa"), "fasta_path")
    with pytest.raises(InputFileError, match="cannot report it"):
        write_new_file(str(root / "new.fa"), ">n\nA\n", "out_path")
    assert not (root / "new.fa").exists()
    set_roots(None)
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"


def test_roots_from_the_environment_and_from_set_roots(dirs, monkeypatch, tmp_path):
    root, outside = dirs
    second = tmp_path / "second"
    second.mkdir()
    (second / "two.fa").write_text(">b\nGG\n")
    alias = tmp_path / "alias"
    alias.symlink_to(second, target_is_directory=True)

    # Several roots, one of them named through a symlink.
    monkeypatch.setenv(ENV, f"{root}{os.pathsep}{alias}")
    assert configured_roots() == [str(root), str(second)]
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"
    assert read_path(str(second / "two.fa"), "fasta_path") == ">b\nGG\n"
    with pytest.raises(PathOutsideRootsError, match=OUTSIDE):
        read_path(str(outside / "secret.fa"), "fasta_path")

    # set_roots (the --root flag) overrides the environment; None restores it.
    set_roots([str(second)])
    with pytest.raises(PathOutsideRootsError, match=OUTSIDE):
        read_path(str(root / "in.fa"), "fasta_path")
    set_roots(None)
    assert read_path(str(root / "in.fa"), "fasta_path") == ">a\nACGT\n"

    # An empty variable is no policy, not an empty allow-list.
    monkeypatch.setenv(ENV, "  ")
    assert configured_roots() is None
    assert read_path(str(outside / "secret.fa"), "fasta_path") == ">s\nTTTT\n"


def test_the_refusal_reaches_the_model_through_the_server(dirs):
    root, outside = dirs
    (root / "trees.nwk").write_text(f"{TRUE_NEWICK}\n{OTHER}\n")
    (outside / "trees.nwk").write_text(f"{TRUE_NEWICK}\n{OTHER}\n")
    set_roots([str(root)])
    mcp = build_server()

    async def run():
        ok = await mcp.call_tool(
            "compare_trees", {"trees_path": str(root / "trees.nwk")}
        )
        with pytest.raises(
            ToolError,
            match=f"^Error executing tool compare_trees: trees_path: .*{OUTSIDE}",
        ):
            await mcp.call_tool(
                "compare_trees", {"trees_path": str(outside / "trees.nwk")}
            )
        return ok

    assert anyio.run(run).structured_content["robinson_foulds"] == 4


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


def _compare_over_stdio(extra_args: list[str], cwd, paths: list[str]) -> list[dict]:
    """Start the real entry point with `extra_args`; one compare_trees per path."""
    with tempfile.TemporaryFile(mode="w+") as err:
        proc = subprocess.Popen(
            [sys.executable, "-m", "phylokit_mcp.server", *extra_args],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=err,
            text=True,
            bufsize=1,
            cwd=cwd,
        )
        try:
            _rpc(proc, {
                "jsonrpc": "2.0", "id": 1, "method": "initialize",
                "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                           "clientInfo": {"name": "test", "version": "0"}},
            })  # fmt: skip
            _rpc(proc, {"jsonrpc": "2.0", "method": "notifications/initialized"})
            return [
                _rpc(
                    proc,
                    {
                        "jsonrpc": "2.0",
                        "id": 2 + i,
                        "method": "tools/call",
                        "params": {
                            "name": "compare_trees",
                            "arguments": {"trees_path": p},
                        },
                    },
                )["result"]
                for i, p in enumerate(paths)
            ]
        finally:
            proc.stdin.close()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                # A server that ignores EOF must not outlive the test.
                proc.kill()
                proc.wait()
                raise


def test_a_real_server_started_with_root_reads_only_inside_it(dirs):
    """Over stdio, as a client launches it. The same outside file is read by a
    server started without `--root`, so the refusal is the flag's doing."""
    root, outside = dirs
    (root / "trees.nwk").write_text(f"{TRUE_NEWICK}\n{OTHER}\n")
    (outside / "trees.nwk").write_text(f"{TRUE_NEWICK}\n{OTHER}\n")
    paths = ["trees.nwk", str(outside / "trees.nwk"), "../outside/trees.nwk"]

    inside, absolute, dotdot = _compare_over_stdio(["--root", str(root)], root, paths)
    assert inside.get("isError") is not True, inside
    assert inside["structuredContent"]["robinson_foulds"] == 4
    for refused in (absolute, dotdot):
        assert refused["isError"] is True
        assert OUTSIDE in refused["content"][0]["text"]
        assert str(root) in refused["content"][0]["text"]  # the policy is named

    open_server = _compare_over_stdio([], root, paths)
    assert [r.get("isError") is True for r in open_server] == [False, False, False]
    assert open_server[1]["structuredContent"]["robinson_foulds"] == 4
