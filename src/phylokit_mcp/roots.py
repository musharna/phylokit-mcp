"""Optional allow-list of directories for every path a tool is given.

Unset, the server opens any path its user can; SECURITY.md states that trust
boundary. Set, with `phylokit-mcp --root DIR` (repeatable) or `PHYLOKIT_ROOTS`
(`os.pathsep`-separated), every `fasta_path` and `trees_path` must be a file
inside one of the roots, and `out_path` must create its file inside one.

This module holds the policy; `inputs.py` applies it and words the refusals.
There are two checks, because a name can be re-pointed between a check and an
open:

- `inside(os.path.realpath(path))` judges the NAME, with `..` and symlinks
  followed first. It runs before anything is opened or stat'ed, so a refusal
  says nothing about a file outside the roots, not even whether it exists.
- `fd_path(fd)` asks the kernel where the descriptor that WAS opened lives.
  That check binds: an ancestor directory swapped for a symlink after the
  first check is followed by `open`, and caught there.

The roots are resolved once, when they are set or when the environment
variable is first read with a given value. A policy re-resolved on every call
can itself be re-pointed between the two checks.
"""

from __future__ import annotations

import os

ENV = "PHYLOKIT_ROOTS"

_roots: list[str] | None = None
_env_snapshot: tuple[str, list[str]] | None = None


def set_roots(roots: list[str] | None) -> None:
    """Set the roots, overriding the environment; None goes back to it."""
    global _roots
    _roots = None if roots is None else [os.path.realpath(r) for r in roots]


def configured_roots() -> list[str] | None:
    """The active roots as real paths, or None when paths are unrestricted."""
    global _env_snapshot
    if _roots is not None:
        return list(_roots)
    raw = os.environ.get(ENV)
    if raw is None or not raw.strip():
        return None
    if _env_snapshot is None or _env_snapshot[0] != raw:
        _env_snapshot = (
            raw,
            [os.path.realpath(r) for r in raw.split(os.pathsep) if r.strip()],
        )
    return list(_env_snapshot[1])


def inside(real: str, roots: list[str]) -> bool:
    """Whether the real path `real` is one of `roots` or lies under one."""
    for root in roots:
        try:
            if os.path.commonpath([root, real]) == root:
                return True
        except ValueError:  # different drives on Windows
            continue
    return False


class CannotLocateDescriptorError(RuntimeError):
    """This platform cannot say where an open descriptor lives."""


def fd_path(fd: int) -> str:
    """Where the open descriptor `fd` lives, from the kernel, not from a name.

    Linux answers through /proc, macOS through F_GETPATH. Anywhere else there
    is no way to ask, and a containment check that cannot be made is refused,
    never skipped.
    """
    try:
        return os.readlink(f"/proc/self/fd/{fd}")
    except OSError:
        pass
    try:
        import fcntl
    except ImportError:
        getpath = None
    else:
        getpath = getattr(fcntl, "F_GETPATH", None)  # macOS only
    if getpath is None:
        raise CannotLocateDescriptorError(
            "this platform cannot report where an open file lives, so the "
            f"configured roots cannot be enforced on it. Run without {ENV} / "
            "--root, or on Linux or macOS."
        )
    raw = fcntl.fcntl(fd, getpath, bytes(1024))
    return os.fsdecode(raw.split(b"\0", 1)[0])
