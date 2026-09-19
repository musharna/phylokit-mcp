"""Multiple sequence alignment, by handing unaligned sequences to MAFFT.

IQ-TREE arrives as a wheel (piqtree). MAFFT does not: it is a system binary, so
it is located on PATH at call time and its absence is a refusal that names the
fix, not an import error at startup — the five tree tools must keep working on
a machine that has no aligner.

Sequences reach MAFFT through a file in a private temporary directory and the
command is an argument list. Nothing the caller supplies is ever part of a
command line.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from .alignment import (
    _DNA,
    _NAME_OK,
    _PROTEIN,
    MAX_SITES,
    MAX_TAXA,
    MOLTYPES,
    AlignmentError,
    parse_fasta,
)

# MAFFT needs two sequences to have anything to align. This is deliberately
# lower than alignment.MIN_TAXA (4): a pairwise alignment is a legitimate
# result, it just cannot go on to `infer_tree`, and the tool says so.
MIN_SEQUENCES = 2
# The server is synchronous, so an alignment that does not finish is a hung
# server. The cap is reported as a refusal rather than left to the client.
TIMEOUT_SECONDS = 600
# Gap and missing-data symbols. MAFFT strips '-' from its input silently, which
# would make "the output, degapped, equals the input" false for a reason the
# caller never sees, so they are refused up front instead.
_GAP_CHARS = set("-.?")

_MOLTYPE_FLAG = {"dna": "--nuc", "protein": "--amino"}

_INSTALL_HINT = (
    "align_sequences needs the MAFFT binary on PATH and it was not found. "
    "Install it with `apt install mafft`, `brew install mafft` or "
    "`conda install -c bioconda mafft`, then restart the server. The tree tools "
    "do not need it and keep working without it."
)


class MafftUnavailableError(ValueError):
    """MAFFT is not installed. A refusal: the message is the instruction."""


class MafftFailedError(RuntimeError):
    """MAFFT ran and did not produce a usable alignment. A crash, not a refusal."""


def find_mafft() -> str | None:
    return shutil.which("mafft")


def mafft_version() -> str | None:
    """MAFFT's own version string, or None when it is not installed."""
    exe = find_mafft()
    if exe is None:
        return None
    # `mafft --version` writes to STDERR, and exits 0 on some builds and 1 on
    # others, so neither the stream nor the status is assumed.
    done = subprocess.run(
        [exe, "--version"],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
        stdin=subprocess.DEVNULL,
    )
    text = (done.stderr + done.stdout).strip()
    if not text:
        raise MafftFailedError(
            f"`{exe} --version` printed nothing (exit {done.returncode})."
        )
    return text.splitlines()[0]


def validate_unaligned(seqs: dict[str, str], moltype: str = "dna") -> None:
    """Reject input MAFFT would accept and answer meaninglessly.

    MAFFT exits 0 on a single sequence and returns it unchanged, so the count
    check has to live here: the engine will not make it.
    """
    if moltype not in MOLTYPES:
        raise AlignmentError(
            f"Unknown sequence_type {moltype!r}. Valid: {list(MOLTYPES)}."
        )
    if len(seqs) < MIN_SEQUENCES:
        raise AlignmentError(
            f"Need at least {MIN_SEQUENCES} sequences to align, got {len(seqs)}."
        )
    if len(seqs) > MAX_TAXA:
        raise AlignmentError(
            f"Got {len(seqs)} sequences; this server caps at {MAX_TAXA}."
        )
    for nm, s in seqs.items():
        if not _NAME_OK.match(nm):
            raise AlignmentError(
                f"Sequence name {nm!r} contains characters that are unsafe in Newick "
                "(the tree format uses ,:;() as syntax). Use letters, digits, "
                "underscore, dot or hyphen."
            )
        if not s:
            raise AlignmentError(f"Sequence {nm!r} is empty.")
        if len(s) > MAX_SITES:
            raise AlignmentError(
                f"Sequence {nm!r} is {len(s)} residues; this server caps at "
                f"{MAX_SITES}."
            )

    present = {c for s in seqs.values() for c in s}
    gaps = present & _GAP_CHARS
    if gaps:
        raise AlignmentError(
            f"Input contains gap or missing-data characters {sorted(gaps)}. "
            "align_sequences takes UNALIGNED sequences; remove the gaps first. "
            "If the sequences are already aligned, pass them to infer_tree as "
            "they are."
        )
    alphabet = _DNA if moltype == "dna" else _PROTEIN
    bad = present - alphabet
    if bad:
        hint = (
            " If these are protein sequences, pass sequence_type='protein': the "
            "molecule type is declared, never guessed."
            if moltype == "dna"
            else ""
        )
        raise AlignmentError(
            f"Unrecognised characters for {moltype}: {sorted(bad)[:8]}.{hint}"
        )


def mafft_args(moltype: str) -> list[str]:
    """The options, without the executable or the input path.

    `--thread 1` for the same reason engine.py pins OMP_NUM_THREADS: a result
    that depends on scheduling is not reproducible. The molecule flag is passed
    explicitly so MAFFT never sniffs it. `--preservecase` and `--inputorder`
    keep the output comparable to the input, which `align` then checks.
    """
    return [
        "--auto",
        "--thread",
        "1",
        "--preservecase",
        "--inputorder",
        _MOLTYPE_FLAG[moltype],
    ]


def align(seqs: dict[str, str], moltype: str = "dna") -> dict[str, str]:
    """Align with MAFFT and verify the result is an alignment OF THE INPUT."""
    validate_unaligned(seqs, moltype)
    exe = find_mafft()
    if exe is None:
        raise MafftUnavailableError(_INSTALL_HINT)

    with tempfile.TemporaryDirectory(prefix="phylokit-mafft-") as tmp:
        infile = Path(tmp) / "input.fasta"
        infile.write_text("".join(f">{nm}\n{s}\n" for nm, s in seqs.items()))
        cmd = [exe, *mafft_args(moltype), str(infile)]
        try:
            done = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=TIMEOUT_SECONDS,
                check=False,
                stdin=subprocess.DEVNULL,
                cwd=tmp,
            )
        except subprocess.TimeoutExpired as exc:
            raise AlignmentError(
                f"MAFFT did not finish within {TIMEOUT_SECONDS} s on "
                f"{len(seqs)} sequences. This server is synchronous; align an "
                "input this large outside it."
            ) from exc

    if done.returncode != 0:
        raise MafftFailedError(
            f"MAFFT exited {done.returncode}. Command: {cmd}. "
            f"stderr:\n{done.stderr[-4000:]}"
        )
    try:
        aligned = parse_fasta(done.stdout)
    except AlignmentError as exc:
        raise MafftFailedError(
            f"MAFFT exited 0 but its output is not FASTA ({exc}). "
            f"stderr:\n{done.stderr[-4000:]}"
        ) from exc

    # MAFFT drops characters it does not recognise without saying so. An
    # alignment that is no longer of the caller's sequences would go on to
    # produce a confident tree of something else, so it is checked, not trusted.
    if list(aligned) != list(seqs):
        raise MafftFailedError(
            f"MAFFT returned names {list(aligned)[:6]} for input names "
            f"{list(seqs)[:6]}; names or order were not preserved."
        )
    for nm, row in aligned.items():
        if row.replace("-", "") != seqs[nm]:
            raise MafftFailedError(
                f"MAFFT's row for {nm!r}, with gaps removed, is not the input "
                "sequence: residues were altered or dropped."
            )
    if len({len(r) for r in aligned.values()}) != 1:
        raise MafftFailedError("MAFFT returned rows of unequal length.")
    return aligned
