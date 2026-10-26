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
- the path names a regular file (a FIFO or device would block or never end);
- it is at most MAX_INPUT_BYTES, read with a bound so a file that grows after
  the size check cannot slip past it;
- it decodes as UTF-8.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

# The largest alignment the tools accept is MAX_TAXA x MAX_SITES = 20 million
# residues (alignment.py). With headers and line breaks that fits in 32 MiB;
# anything bigger would be refused downstream anyway, after being read.
MAX_INPUT_BYTES = 32 * 1024 * 1024


class InputFileError(ValueError):
    """A file-path argument that cannot be used, with the reason named."""


def read_path(path: str, param: str) -> str:
    """The text of the file at `path`, or a refusal naming what is wrong."""
    if not isinstance(path, str) or not path.strip():
        raise InputFileError(f"{param} is empty; give the path of a file.")
    p = Path(path).expanduser()
    where = f" (relative paths resolve against the server's working directory, {os.getcwd()})"
    shown = str(p) if p.is_absolute() else f"{p}{where}"
    try:
        st = os.stat(p)
    except FileNotFoundError:
        raise InputFileError(f"{param}: no such file: {shown}.") from None
    except PermissionError:
        raise InputFileError(f"{param}: permission denied: {shown}.") from None
    except OSError as exc:
        raise InputFileError(f"{param}: cannot stat {shown}: {exc}.") from None
    if not stat.S_ISREG(st.st_mode):
        kind = "a directory" if stat.S_ISDIR(st.st_mode) else "not a regular file"
        raise InputFileError(f"{param}: {shown} is {kind}; give the path of a file.")
    if st.st_size > MAX_INPUT_BYTES:
        raise InputFileError(
            f"{param}: {shown} is {st.st_size:,} bytes; the cap is "
            f"{MAX_INPUT_BYTES:,} (the largest alignment the tools accept)."
        )
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
    if not isinstance(path, str) or not path.strip():
        raise InputFileError(f"{param} is empty; give the path of a new file.")
    p = Path(path).expanduser()
    where = f" (relative paths resolve against the server's working directory, {os.getcwd()})"
    shown = str(p) if p.is_absolute() else f"{p}{where}"
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
