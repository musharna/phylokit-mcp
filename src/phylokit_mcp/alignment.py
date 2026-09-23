"""Turning caller-supplied sequences into an alignment, and refusing the cases
where inference would return a confident answer to a question the data cannot
address.

The checks here are the ones whose absence produces a *plausible* tree rather
than an error. Unequal sequence lengths raise loudly in the engine; a
234-taxon alignment of 40 informative sites does not.
"""

from __future__ import annotations

import functools
import re
from dataclasses import dataclass

MIN_TAXA = 4  # below this an unrooted tree has no internal edge to support
MAX_TAXA = 200
MAX_SITES = 100_000

MOLTYPES: tuple[str, ...] = ("dna", "protein")
# States that count toward parsimony signal, per molecule type. Ambiguity codes
# and gaps are excluded from both: they are not evidence of a shared state.
_INFORMATIVE_STATES = {
    "dna": set("ACGTU"),
    "protein": set("ACDEFGHIKLMNPQRSTVWY"),
}
_NAME_OK = re.compile(r"^[A-Za-z0-9_.\-]+$")


class AlignmentError(ValueError):
    """Raised when the input cannot support inference, with the reason named."""


@functools.cache
def _accepts(moltype: str, char: str) -> bool:
    """Whether cogent3 accepts `char` in a sequence of `moltype`.

    Asked of cogent3 rather than restated here. The previous hand-written sets
    disagreed with it in both molecule types: `.` for DNA, and `*`, `.`, `J`,
    `O` for protein, were passed by validate() and then rejected by cogent3's
    AlphabetError inside the engine call, which the server could only report as
    a crash. This mirrors what cogent3 does to a sequence on the way in --
    `coerce_to` (case folding, U->T for DNA), then the most degenerate alphabet
    -- so the two cannot drift apart again.
    """
    from cogent3 import get_moltype

    mt = get_moltype(moltype)
    raw = char.encode("utf8")
    coerced = mt.coerce_to(raw) if mt.coerce_to else raw
    return bool(mt.is_valid(coerced))


def unrecognised_characters(chars: set[str], moltype: str) -> set[str]:
    """The members of `chars` that cogent3 would reject for `moltype`."""
    return {c for c in chars if not _accepts(moltype, c)}


def stop_codon_hint(bad: set[str]) -> str:
    """The one rejected character with an obvious remedy, named."""
    if "*" not in bad:
        return ""
    return (
        " '*' is a stop codon, and the engine has no state for it: trim stops "
        "(or replace them with X) before aligning or inferring."
    )


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


def validate(seqs: dict[str, str], moltype: str = "dna") -> None:
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
            f"Lengths present — {detail}. If they are unaligned, run align_sequences "
            "first; the tree tools will not guess an alignment."
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

    bad = unrecognised_characters({c for s in seqs.values() for c in s}, moltype)
    if bad:
        hint = (
            " If this is a protein alignment, pass sequence_type='protein': the "
            "molecule type is declared, never guessed, because an alignment of "
            "only A/C/G/T is a valid protein alignment too and guessing wrong "
            "silently fits the wrong substitution model."
            if moltype == "dna"
            else ""
        )
        raise AlignmentError(
            f"Unrecognised characters for {moltype}: {sorted(bad)[:8]}.{hint}"
            f"{stop_codon_hint(bad)}"
        )


def require_phylogenetic_signal(seqs: dict[str, str], moltype: str = "dna") -> None:
    """Refuse an alignment that cannot support ANY tree.

    Deliberately not part of `validate`. That function answers "is this a
    well-formed alignment", and a signal-free alignment is perfectly well
    formed — `select_substitution_model` and the molecule-type checks have
    legitimate reasons to accept one. This asks the narrower question that only
    tree inference needs answered, so it is called only there.

    Without it, a signal-free alignment reaches IQ-TREE, which cannot fit a
    likelihood and fails inside piqtree with "IQ-TREE output is malformed,
    likelihood not found." That message is upstream's and describes a PARSING
    failure, so a caller reads it as this server being broken and retries — when
    the real answer is that the data cannot support a tree and no retry will
    change it. The condition is knowable before the call.
    """
    if parsimony_informative(seqs, moltype) != 0:
        return
    groups = duplicate_groups(seqs)
    identical = sum(len(g) for g in groups) == len(seqs) and len(groups) == 1
    detail = (
        f"All {len(seqs)} sequences are identical."
        if identical
        else "No site has two different states each appearing in two or more "
        "sequences, so no site can distinguish one grouping from another."
    )
    raise AlignmentError(
        f"This alignment carries no phylogenetic signal. {detail} Every topology "
        "fits it equally well, so there is no tree to infer — this is a property "
        "of the data, not a transient failure, and retrying will not change it. "
        "Supply an alignment with variation shared across taxa."
    )


def parsimony_informative(seqs: dict[str, str], moltype: str = "dna") -> int:
    """Count sites with >=2 states each seen in >=2 taxa.

    This, not alignment length, is the quantity that carries topological signal.
    A 10,000-site alignment of near-identical sequences supports nothing, and
    reporting its length would imply otherwise.
    """
    # Hardcoding "ACGTU" here would count ZERO informative sites in every
    # protein alignment, so the low-signal guard would refuse valid data while
    # appearing to have measured it.
    states = _INFORMATIVE_STATES[moltype]
    names = list(seqs)
    n_sites = len(seqs[names[0]])
    count = 0
    for i in range(n_sites):
        tally: dict[str, int] = {}
        for nm in names:
            ch = seqs[nm][i].upper()
            if ch in states:
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


def summarise(seqs: dict[str, str], moltype: str = "dna") -> AlignmentStats:
    n_sites = len(next(iter(seqs.values())))
    total = sum(len(s) for s in seqs.values())
    gaps = sum(s.count("-") + s.count("?") for s in seqs.values())
    return AlignmentStats(
        n_taxa=len(seqs),
        n_sites=n_sites,
        moltype=moltype,
        n_parsimony_informative=parsimony_informative(seqs, moltype),
        fraction_gaps=round(gaps / total, 4) if total else 0.0,
        duplicate_sequences=duplicate_groups(seqs),
    )


def to_cogent3(seqs: dict[str, str], moltype: str = "dna"):
    """The alignment as cogent3 holds it -- the input-parsing boundary.

    cogent3 reports a character it cannot place with its own AlphabetError,
    which is not a ValueError, so the server would mask it as a crash. It is
    the caller's data that is wrong, so it is translated here into the
    refusal type, with the offending characters named.
    """
    from cogent3 import make_aligned_seqs
    from cogent3.core.alphabet import AlphabetError

    if moltype not in MOLTYPES:
        raise AlignmentError(
            f"Unknown sequence_type {moltype!r}. Valid: {list(MOLTYPES)}."
        )
    try:
        return make_aligned_seqs(seqs, moltype=moltype)
    except AlphabetError as exc:
        bad = unrecognised_characters({c for s in seqs.values() for c in s}, moltype)
        raise AlignmentError(
            f"Unrecognised characters for {moltype}: {sorted(bad)[:8]} ({exc})."
        ) from exc
