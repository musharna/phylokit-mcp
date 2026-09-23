"""Bad caller input must come back as a refusal that says what was wrong.

Findings M2-M5, L6 and L7 of the 2026-09-22 MCP-server audit. Every test here
goes through a real MCP client session against `build_server()`, because the
failure being fixed lives AT that boundary: a tool function that raises a
perfectly clear cogent3 or piqtree error still reaches the caller as a bare
`Error executing tool <name>` if the error is not a ValueError. Calling the
tool functions directly would see the clear error and pass.

Each negative case sits beside a positive control in the same test, so a
harness that failed everything would not read as "refused correctly".
"""

from __future__ import annotations

import string
from typing import Any

import anyio
import pytest
from cogent3.core.alphabet import AlphabetError
from mcp.client.client import Client

from phylokit_mcp import server
from phylokit_mcp.alignment import (
    AlignmentError,
    parse_fasta,
    summarise,
    to_cogent3,
    validate,
)

# Restated rather than imported: the largest seed IQ-TREE's C int can hold, and
# the independent source a changed engine constant would disagree with.
SEED_MAX = 2**31 - 1

A4 = (
    ">a\nACGTACGTAAACGTTT\n>b\nACGTACGTAAACGTTA\n"
    ">c\nTCGAACGTAAACGTAA\n>d\nTCGAACGTAAACGAAA\n"
)
TREE5 = "((a:0.3,b:0.3):0.3,(c:0.3,d:0.3):0.3,e:0.3);"
PROT4 = ">a\nMKVLAAGWKLAA\n>b\nMKVLAGGWKLAC\n>c\nMRVLSAGWRLAC\n>d\nMRVLSGGWRLAA\n"


def _call(tool: str, args: dict[str, Any]):
    async def go():
        async with Client(server.build_server()) as c:
            return await c.call_tool(tool, args)

    return anyio.run(go)


def _text(result) -> str:
    return result.content[0].text if result.content else ""


def _assert_refused(tool: str, args: dict[str, Any], *needles: str) -> str:
    r = _call(tool, args)
    text = _text(r)
    assert r.is_error, f"{tool} accepted {args!r}: {text[:200]}"
    # The masked form is exactly this string with nothing after it.
    assert text != f"Error executing tool {tool}", (
        f"{tool} refused {args!r} as an opaque crash, reason dropped"
    )
    for n in needles:
        assert n in text, f"expected {n!r} in refusal, got: {text[:300]}"
    return text


def _assert_ok(tool: str, args: dict[str, Any]):
    r = _call(tool, args)
    assert not r.is_error, f"{tool} failed on a valid call: {_text(r)[:300]}"
    return r.structured_content


# ---------------------------------------------------------------------------
# The class: third-party input errors surface with their reason
# ---------------------------------------------------------------------------

MASKED_BEFORE = [
    # (tool, args, text the refusal must carry). Each of these returned the
    # bare `Error executing tool <name>` on 657e5ff.
    ("infer_tree", {"fasta": A4, "replicates": 20, "seed": 2**40}, "seed must be"),
    ("infer_tree", {"fasta": A4.replace("AAACGTTA", "AA.CGTTA")}, "'.'"),
    ("infer_tree", {"fasta": A4, "model": "LG", "replicates": 20}, "protein model"),
    (
        "infer_tree",
        {"fasta": PROT4, "model": "GTR+G", "sequence_type": "protein"},
        "dna model",
    ),
    ("compare_trees", {"newick_a": "((a,b),(c,d)", "newick_b": TREE5}, "Newick"),
    ("simulate_alignment", {"newick": "((a,b", "length": 50}, "Newick"),
    ("simulate_alignment", {"newick": TREE5, "seed": 2**40}, "seed must be"),
    ("select_substitution_model", {"fasta": A4, "seed": 2**31}, "seed must be"),
]


@pytest.mark.parametrize(("tool", "args", "needle"), MASKED_BEFORE)
def test_input_errors_from_the_libraries_are_refusals_with_a_reason(tool, args, needle):
    _assert_refused(tool, args, needle)


def test_the_same_tools_still_succeed_on_valid_input():
    """Positive control for the parametrised refusals above."""
    _assert_ok("infer_tree", {"fasta": A4, "replicates": 20, "seed": SEED_MAX})
    _assert_ok("compare_trees", {"newick_a": TREE5, "newick_b": TREE5})
    _assert_ok("simulate_alignment", {"newick": TREE5, "length": 50})


# ---------------------------------------------------------------------------
# M2: validate() accepts exactly what cogent3 accepts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("moltype", ["dna", "protein"])
def test_validate_and_cogent3_agree_on_every_printable_character(moltype):
    """The two alphabets were separate lists and disagreed.

    `.` (both), and `*`, `J`, `O` (protein) passed validate() and then failed
    inside cogent3, which the server could only report as a crash. Every
    printable character is tried against both, so any future disagreement in
    either direction fails here.
    """
    base = "ACGTACGTAC" if moltype == "dna" else "ACDEFGHIKL"
    accepted_by_both = 0
    for ch in sorted(set(string.printable) - set(string.whitespace)):
        seqs = {n: base for n in "abcd"}
        seqs["a"] = ch + base[1:]
        try:
            validate(seqs, moltype)
            ours = True
        except AlignmentError:
            ours = False
        try:
            to_cogent3(seqs, moltype)
            theirs = True
        # cogent3's own error before the boundary translated it; either way
        # this is cogent3 rejecting the character.
        except (AlignmentError, AlphabetError):
            theirs = False
        assert ours == theirs, (
            f"{ch!r} for {moltype}: validate() says {ours}, cogent3 says {theirs}"
        )
        accepted_by_both += ours
    # Positive control: an alphabet that rejected everything would agree too.
    assert accepted_by_both >= (16 if moltype == "dna" else 25)


def test_align_output_is_always_accepted_by_infer_tree():
    """align_sequences said ready_for_infer_tree=True for output with `*`.

    MAFFT aligns stop codons happily; infer_tree then crashed on them. The
    refusal now happens before MAFFT runs, with the same alphabet.
    """
    stop = PROT4.replace("KLAA\n", "KLAA*\n").replace("KLAC\n", "KLAC*\n")
    stop = stop.replace("RLAC\n", "RLAC*\n").replace("RLAA\n", "RLAA*\n")
    _assert_refused(
        "align_sequences", {"fasta": stop, "sequence_type": "protein"}, "stop codon"
    )

    # Positive control: the same sequences without stops align and are ready,
    # and "ready" is borne out by infer_tree accepting the output unchanged.
    out = _assert_ok("align_sequences", {"fasta": PROT4, "sequence_type": "protein"})
    assert out["ready_for_infer_tree"] is True
    _assert_ok(
        "infer_tree",
        {
            "fasta": out["fasta"],
            "sequence_type": "protein",
            "model": "LG",
            "replicates": 20,
        },
    )


def test_ready_for_infer_tree_is_false_when_infer_tree_would_refuse():
    """`ready` is infer_tree's own checks on the output, not a taxon count."""
    same = ">a\nACGTACGTAC\n>b\nACGTACGTAC\n>c\nACGTACGTAC\n>d\nACGTACGTAC\n"
    out = _assert_ok("align_sequences", {"fasta": same})
    assert out["ready_for_infer_tree"] is False
    codes = [w["code"] for w in out["warnings"]]
    assert "infer_tree_would_refuse" in codes
    _assert_refused("infer_tree", {"fasta": out["fasta"]}, "no phylogenetic signal")


# ---------------------------------------------------------------------------
# M4: duplicate tip names
# ---------------------------------------------------------------------------


def test_a_tree_that_names_a_tip_twice_is_refused_not_renamed():
    """cogent3 renames a repeat to `a.2`, and the comparison used to run on it."""
    text = _assert_refused(
        "compare_trees",
        {"newick_a": "((a,a),(c,d),e);", "newick_b": "((a,c),(b,d),e);"},
        "more than once",
        "['a']",
    )
    assert "a.2|" not in text
    _assert_refused(
        "simulate_alignment",
        {"newick": "((a:0.1,a:0.1):0.1,(c:0.1,d:0.1):0.1,e:0.1);"},
        "more than once",
    )
    # Positive control: distinct names that merely share a prefix are fine.
    out = _assert_ok(
        "compare_trees",
        {"newick_a": "((a,a2),(c,d),e);", "newick_b": "((a,c),(a2,d),e);"},
    )
    assert out["shared_taxa"] == 5


# ---------------------------------------------------------------------------
# M5: a protein simulation is summarised as protein
# ---------------------------------------------------------------------------


def test_a_protein_model_simulates_and_summarises_protein():
    out = _assert_ok("simulate_alignment", {"newick": TREE5, "model": "LG"})
    assert out["alignment"]["moltype"] == "protein"
    as_protein = summarise(parse_fasta(out["fasta"]), "protein")
    assert (
        out["alignment"]["n_parsimony_informative"]
        == as_protein.n_parsimony_informative
    )
    # Positive control: a nucleotide model is still summarised as DNA.
    dna = _assert_ok("simulate_alignment", {"newick": TREE5, "model": "JC"})
    assert dna["alignment"]["moltype"] == "dna"


# ---------------------------------------------------------------------------
# L6: simulate_alignment's input limits
# ---------------------------------------------------------------------------


def test_simulate_enforces_the_taxon_cap_infer_tree_enforces():
    from phylokit_mcp.alignment import MAX_TAXA

    def star(n: int) -> str:
        return "(" + ",".join(f"t{i}:0.1" for i in range(n)) + ");"

    _assert_refused(
        "simulate_alignment", {"newick": star(MAX_TAXA + 1), "length": 20}, "tips"
    )
    out = _assert_ok("simulate_alignment", {"newick": star(MAX_TAXA), "length": 20})
    assert out["alignment"]["n_taxa"] == MAX_TAXA


def test_simulate_refuses_a_tree_without_branch_lengths():
    """Without lengths every sequence came back identical, reported as success."""
    _assert_refused(
        "simulate_alignment",
        {"newick": "((a,b),(c,d),e);", "length": 50},
        "no branch length",
    )
    out = _assert_ok("simulate_alignment", {"newick": TREE5, "length": 50})
    assert out["alignment"]["duplicate_sequences"] == []


# ---------------------------------------------------------------------------
# L7: seed and replicates are checked before any engine work
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("args", "needle"),
    [
        ({"seed": -5}, "seed must be"),
        ({"seed": SEED_MAX + 1}, "seed must be"),
        ({"replicates": 5}, "replicates must be"),
    ],
)
def test_bad_seed_or_replicates_is_refused_before_the_ml_search(
    monkeypatch, args, needle
):
    """A negative seed used to be refused by numpy -- after the ML search ran."""
    calls: list[str] = []
    real = server.build_ml_tree

    def spy(*a, **k):
        calls.append("ml")
        return real(*a, **k)

    monkeypatch.setattr(server, "build_ml_tree", spy)
    _assert_refused("infer_tree", {"fasta": A4, "replicates": 20, **args}, needle)
    assert calls == [], "the ML search ran before the argument was refused"

    # Positive control: the spy does see a valid call, so `calls == []` above is
    # evidence and not an unhooked spy.
    _assert_ok("infer_tree", {"fasta": A4, "replicates": 20, "seed": 0})
    assert calls == ["ml"]
