"""Turning caller-supplied sequences into an alignment, and refusing the cases
where inference would return a confident answer to a question the data cannot
address.

The checks here are the ones whose absence produces a *plausible* tree rather
than an error. Unequal sequence lengths raise loudly in the engine; a
234-taxon alignment of 40 informative sites does not.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

MIN_TAXA = 4  # below this an unrooted tree has no internal edge to support
MAX_TAXA = 200
MAX_SITES = 100_000

_DNA = set("ACGTUN-?RYSWKMBDHVacgtunryswkmbdhv.")
_NAME_OK = re.compile(r"^[A-Za-z0-9_.\-]+$")


class AlignmentError(ValueError):
    """Raised when the input cannot support inference, with the reason named."""


@dataclass(frozen=True)
class AlignmentStats:
    n_taxa: int
    n_sites: int
    moltype: str
    n_parsimony_informative: int
    fraction_gaps: float
    duplicate_sequences: list[list[str]]


def parse_fasta(text: str) -> dict[str, str]:
    """Minimal FASTA reader. Sequence names are taken up to the first whitespace."""
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0] if len(line) > 1 else ""
            if not name:
                raise AlignmentError("A FASTA record has an empty name (a bare '>').")
            if name in seqs:
                raise AlignmentError(
                    f"Duplicate sequence name {name!r}. Names must be unique — "
                    "a tree cannot have two tips with the same label."
                )
            chunks = []
        else:
            if name is None:
                raise AlignmentError(
                    "The alignment does not start with a '>' header line; "
                    "this does not look like FASTA."
                )
            chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    if not seqs:
        raise AlignmentError("No sequences found in the input.")
    return seqs


def validate(seqs: dict[str, str]) -> None:
    """Reject inputs that would otherwise yield a confident, meaningless tree."""
    if len(seqs) < MIN_TAXA:
        raise AlignmentError(
            f"Need at least {MIN_TAXA} sequences, got {len(seqs)}. An unrooted tree "
            "on three or fewer taxa has only one topology, so there is nothing to "
            "infer and no clade that could be supported."
        )
    if len(seqs) > MAX_TAXA:
        raise AlignmentError(
            f"Got {len(seqs)} sequences; this server caps at {MAX_TAXA}. Bootstrap "
            "cost grows with taxon count and the server is synchronous."
        )

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        by_len: dict[int, list[str]] = {}
        for nm, s in seqs.items():
            by_len.setdefault(len(s), []).append(nm)
        detail = "; ".join(
            f"{ln} sites: {', '.join(sorted(names)[:4])}"
            for ln, names in sorted(by_len.items())
        )
        raise AlignmentError(
            "Sequences are not all the same length, so this is not an alignment. "
            f"Lengths present — {detail}. Align the sequences first; this server "
            "infers trees, it does not align."
        )

    (n_sites,) = lengths
    if n_sites == 0:
        raise AlignmentError("Sequences are empty (zero sites).")
    if n_sites > MAX_SITES:
        raise AlignmentError(
            f"Alignment is {n_sites} sites; this server caps at {MAX_SITES}."
        )

    for nm in seqs:
        if not _NAME_OK.match(nm):
            raise AlignmentError(
                f"Sequence name {nm!r} contains characters that are unsafe in Newick "
                "(the tree format uses ,:;() as syntax). Use letters, digits, "
                "underscore, dot or hyphen."
            )

    bad = {c for s in seqs.values() for c in s} - _DNA
    if bad:
        raise AlignmentError(
            f"Unrecognised characters for DNA: {sorted(bad)[:8]}. This server "
            "currently handles nucleotide alignments only."
        )


def parsimony_informative(seqs: dict[str, str]) -> int:
    """Count sites with >=2 states each seen in >=2 taxa.

    This, not alignment length, is the quantity that carries topological signal.
    A 10,000-site alignment of near-identical sequences supports nothing, and
    reporting its length would imply otherwise.
    """
    names = list(seqs)
    n_sites = len(seqs[names[0]])
    count = 0
    for i in range(n_sites):
        tally: dict[str, int] = {}
        for nm in names:
            ch = seqs[nm][i].upper()
            if ch in "ACGTU":
                tally[ch] = tally.get(ch, 0) + 1
        if sum(1 for v in tally.values() if v >= 2) >= 2:
            count += 1
    return count


def duplicate_groups(seqs: dict[str, str]) -> list[list[str]]:
    """Groups of taxa with byte-identical sequences.

    Identical sequences cannot be resolved relative to one another, so any
    branching order the tree shows between them is arbitrary — but it is drawn
    with the same confidence as a real one.
    """
    by_seq: dict[str, list[str]] = {}
    for nm, s in seqs.items():
        by_seq.setdefault(s.upper(), []).append(nm)
    return sorted(
        (sorted(g) for g in by_seq.values() if len(g) > 1), key=lambda g: g[0]
    )


def summarise(seqs: dict[str, str]) -> AlignmentStats:
    n_sites = len(next(iter(seqs.values())))
    total = sum(len(s) for s in seqs.values())
    gaps = sum(s.count("-") + s.count("?") for s in seqs.values())
    return AlignmentStats(
        n_taxa=len(seqs),
        n_sites=n_sites,
        moltype="dna",
        n_parsimony_informative=parsimony_informative(seqs),
        fraction_gaps=round(gaps / total, 4) if total else 0.0,
        duplicate_sequences=duplicate_groups(seqs),
    )


def to_cogent3(seqs: dict[str, str]):
    from cogent3 import make_aligned_seqs

    return make_aligned_seqs(seqs, moltype="dna")
