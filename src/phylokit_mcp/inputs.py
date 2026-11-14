"""Tool inputs given as a file path instead of inline text.

An alignment passed inline is paid for twice: the client reads the file into
its context, then writes every residue back out as the tool call's argument.
Measured on this repo's agent benchmark (10-16 transcripts, `diagnose.py`), the
server arm wrote a median of 16,000 characters into tool calls per question,
and output tokens were 55% of its cost. A path is a few dozen characters.

The path is opened by THIS process, so it is resolved against the server's own
working directory when relative, and every refusal names that directory: a
client whose working directory differs from the server's can see why a
relative path missed and send an absolute one.

Each check below is a refusal (a ValueError, which the server surfaces with its
message), never a silent fallback:

- exactly one of the inline text and the path is given;
- the path has no NUL byte, and a leading `~name` names a home directory the
  host can find;
- the path names a regular file (a FIFO or device would block or never end);
- it is at most MAX_INPUT_BYTES, read with a bound so a file that grows after
  the size check cannot slip past it;
- it decodes as UTF-8;
- when roots are configured (`roots.py`), it lies inside one of them: judged
  on the name before anything is opened or stat'ed, and again on the
  descriptor that was opened.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

from .roots import ENV as ROOTS_ENV
from .roots import CannotLocateDescriptorError, configured_roots, fd_path, inside

# The largest alignment the tools accept is MAX_TAXA x MAX_SITES = 20 million
# residues (alignment.py). With headers and line breaks that fits in 32 MiB;
# anything bigger would be refused downstream anyway, after being read.
MAX_INPUT_BYTES = 32 * 1024 * 1024


# POSIX-only open flags, 0 where the platform has none: there the open either
# fails or the descriptor cannot be located, and both are refusals.
_O_NONBLOCK = getattr(os, "O_NONBLOCK", 0)
_O_DIRECTORY = getattr(os, "O_DIRECTORY", 0)


class InputFileError(ValueError):
    """A file-path argument that cannot be used, with the reason named."""


class PathOutsideRootsError(InputFileError):
    """A path argument that leads outside every configured root."""


def _require_name_inside(p: Path, shown: str, param: str, roots: list[str]) -> None:
    """Refuse `p` unless it resolves inside the roots. Opens and stats nothing,
    so the refusal is the same whether or not the file outside exists."""
    real = os.path.realpath(p)
    if not inside(real, roots):
        raise PathOutsideRootsError(
            f"{param}: {shown} resolves to {real}, which is outside the "
            f"configured roots {roots}. This server only opens files under "
            f"those directories ({ROOTS_ENV} / --root)."
        )


def _require_fd_inside(
    fd: int, shown: str, param: str, roots: list[str], nothing: str
) -> str:
    """Refuse unless the file or directory BEHIND `fd` is inside the roots;
    return where it is. This is the check that binds: see roots.py."""
    try:
        real = fd_path(fd)
    except CannotLocateDescriptorError as exc:
        raise InputFileError(f"{param}: {exc}") from None
    if not inside(real, roots):
        raise PathOutsideRootsError(
            f"{param}: {shown} was opened as {real}, which is outside the "
            f"configured roots {roots}: the path changed between the check "
            f"and the open. {nothing}"
        )
    return real


def _as_path(path: str, param: str, what: str) -> tuple[Path, str]:
    """`path` as a Path with `~` expanded, and how to name it in a refusal.

    Two inputs fail inside the standard library before any errno exists, so
    the OSError clauses below never see them (fuzz issue 55): a NUL byte makes
    CPython raise its own ValueError at the system-call boundary, and `~name`
    for a user the host does not have makes `expanduser` raise RuntimeError.
    Both are the caller's argument at fault, so both are refused here by name.
    """
    if not isinstance(path, str) or not path.strip():
        raise InputFileError(f"{param} is empty; give the path of {what}.")
    if "\x00" in path:
        raise InputFileError(
            f"{param} contains a NUL byte, which no file name can; give the path of {what}."
        )
    try:
        p = Path(path).expanduser()
    except RuntimeError:
        head = path.split("/", 1)[0]
        raise InputFileError(
            f"{param}: cannot expand {head!r} in {path!r}: that home directory "
            "cannot be determined on the server's host. Give an absolute path, or "
            f"one relative to the server's working directory, {os.getcwd()}."
        ) from None
    where = f" (relative paths resolve against the server's working directory, {os.getcwd()})"
    return p, str(p) if p.is_absolute() else f"{p}{where}"


def _require_regular_and_small(st: os.stat_result, shown: str, param: str) -> None:
    if not stat.S_ISREG(st.st_mode):
        kind = "a directory" if stat.S_ISDIR(st.st_mode) else "not a regular file"
        raise InputFileError(f"{param}: {shown} is {kind}; give the path of a file.")
    if st.st_size > MAX_INPUT_BYTES:
        raise InputFileError(
            f"{param}: {shown} is {st.st_size:,} bytes; the cap is "
            f"{MAX_INPUT_BYTES:,} (the largest alignment the tools accept)."
        )


def _read_inside(p: Path, shown: str, param: str, roots: list[str]) -> bytes:
    """The bytes of `p` under configured roots, at most one past the cap.

    Open first, then judge the descriptor, then look at what it is: nothing
    about the file (its kind, its size, a byte of it) is examined until the
    thing that was opened is known to be inside the roots. O_NONBLOCK so that
    opening a FIFO returns at once and is refused as not a regular file; it
    changes nothing for a regular file.
    """
    _require_name_inside(p, shown, param, roots)
    try:
        fd = os.open(p, os.O_RDONLY | _O_NONBLOCK)
    except FileNotFoundError:
        raise InputFileError(f"{param}: no such file: {shown}.") from None
    except PermissionError:
        raise InputFileError(f"{param}: permission denied: {shown}.") from None
    except OSError as exc:
        raise InputFileError(f"{param}: cannot read {shown}: {exc}.") from None
    try:
        _require_fd_inside(fd, shown, param, roots, "Nothing was read.")
        _require_regular_and_small(os.fstat(fd), shown, param)
        chunks, left = [], MAX_INPUT_BYTES + 1
        while left > 0 and (chunk := os.read(fd, min(left, 1 << 20))):
            chunks.append(chunk)
            left -= len(chunk)
        return b"".join(chunks)
    except OSError as exc:
        raise InputFileError(f"{param}: cannot read {shown}: {exc}.") from None
    finally:
        os.close(fd)


def read_path(path: str, param: str) -> str:
    """The text of the file at `path`, or a refusal naming what is wrong."""
    p, shown = _as_path(path, param, "a file")
    roots = configured_roots()
    if roots is not None:
        data = _read_inside(p, shown, param, roots)
    else:
        try:
            st = os.stat(p)
        except FileNotFoundError:
            raise InputFileError(f"{param}: no such file: {shown}.") from None
        except PermissionError:
            raise InputFileError(f"{param}: permission denied: {shown}.") from None
        except OSError as exc:
            raise InputFileError(f"{param}: cannot stat {shown}: {exc}.") from None
        _require_regular_and_small(st, shown, param)
        try:
            with open(p, "rb") as fh:
                data = fh.read(MAX_INPUT_BYTES + 1)
        except OSError as exc:
            raise InputFileError(f"{param}: cannot read {shown}: {exc}.") from None
    if len(data) > MAX_INPUT_BYTES:
        raise InputFileError(
            f"{param}: {shown} grew past {MAX_INPUT_BYTES:,} bytes while being read."
        )
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputFileError(
            f"{param}: {shown} is not UTF-8 text ({exc.reason} at byte {exc.start})."
        ) from None


def text_or_path(
    text: str | None, path: str | None, text_name: str, path_name: str
) -> str:
    """Exactly one of an inline argument and its path form, as text."""
    if text is not None and path is not None:
        raise InputFileError(
            f"Give {text_name} or {path_name}, not both: with both there is no "
            "way to know which one you meant."
        )
    if path is not None:
        return read_path(path, path_name)
    if text is None:
        raise InputFileError(
            f"Give {text_name} (the text inline) or {path_name} (the path of a file)."
        )
    return text


def write_new_file(path: str, text: str, param: str) -> str:
    """Write `text` to a file that must not exist yet; return its absolute path.

    The output side of the path inputs above, for a client that cannot write a
    file itself (the field benchmark's server arm had only Read, so it retyped
    the alignment `align_sequences` returned and mis-copied a row each time).

    Never over an existing file: `open(..., "x")` refuses atomically, so a
    path that appears between a check and the write is still not clobbered.
    The parent directory must exist; nothing is created on the caller's behalf.
    Relative paths resolve against this process's working directory, as above.
    """
    p, shown = _as_path(path, param, "a new file")
    roots = configured_roots()
    if roots is not None:
        return _write_new_file_inside(p, shown, text, param, roots)
    if not p.parent.is_dir():
        raise InputFileError(f"{param}: no such directory: {p.parent} for {shown}.")
    try:
        with open(p, "x", encoding="utf-8") as fh:
            fh.write(text)
    except FileExistsError:
        raise InputFileError(
            f"{param}: {shown} already exists; give a new path (it is never overwritten)."
        ) from None
    except OSError as exc:
        raise InputFileError(f"{param}: cannot write {shown}: {exc}.") from None
    return str(p.resolve())


def _write_new_file_inside(
    p: Path, shown: str, text: str, param: str, roots: list[str]
) -> str:
    """`write_new_file` under configured roots.

    The file is created relative to a descriptor of its directory, and that
    descriptor is what gets checked: a directory re-pointed to somewhere
    outside the roots after the name check is the directory that was opened,
    and is refused before anything is created. O_EXCL does not follow a
    symlink at the file's own name, so the new file cannot land elsewhere.
    """
    _require_name_inside(p, shown, param, roots)
    if not p.parent.is_dir():
        raise InputFileError(f"{param}: no such directory: {p.parent} for {shown}.")
    try:
        dir_fd = os.open(p.parent, os.O_RDONLY | _O_DIRECTORY)
    except OSError as exc:
        raise InputFileError(f"{param}: cannot write {shown}: {exc}.") from None
    try:
        where = _require_fd_inside(dir_fd, shown, param, roots, "Nothing was written.")
        try:
            fd = os.open(
                p.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666, dir_fd=dir_fd
            )
            with open(fd, "w", encoding="utf-8") as fh:
                fh.write(text)
        except FileExistsError:
            raise InputFileError(
                f"{param}: {shown} already exists; give a new path (it is never overwritten)."
            ) from None
        except OSError as exc:
            raise InputFileError(f"{param}: cannot write {shown}: {exc}.") from None
        return os.path.join(where, p.name)
    finally:
        os.close(dir_fd)
