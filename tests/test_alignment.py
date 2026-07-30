"""Input handling. These reject the cases that would otherwise produce a
confident tree from data that cannot support one — the failures that return a
result rather than an error.
"""

from __future__ import annotations

import pytest

from phylokit_mcp.alignment import (
    AlignmentError,
    duplicate_groups,
    parse_fasta,
    parsimony_informative,
    summarise,
    validate,
)
from phylokit_mcp.server import infer_tree

FOUR = {"a": "ACGTACGT", "b": "ACGTACGA", "c": "TCGTACGT", "d": "TCGAACGT"}


def _fasta(seqs: dict[str, str]) -> str:
    return "".join(f">{k}\n{v}\n" for k, v in seqs.items())


def test_parses_multiline_fasta():
    parsed = parse_fasta(">a\nACGT\nACGT\n>b\nTTTT\nGGGG\n")
    assert parsed == {"a": "ACGTACGT", "b": "TTTTGGGG"}


def test_name_is_taken_up_to_first_whitespace():
    assert set(parse_fasta(">a description here\nACGT\n>b\nACGT\n")) == {"a", "b"}


def test_duplicate_names_are_refused():
    with pytest.raises(AlignmentError, match="Duplicate sequence name"):
        parse_fasta(">a\nACGT\n>a\nTTTT\n")


def test_non_fasta_input_is_refused():
    with pytest.raises(AlignmentError, match="does not look like FASTA"):
        parse_fasta("a ACGT\nb ACGT\n")


def test_ragged_sequences_are_refused_with_the_lengths_named():
    """Unequal lengths mean this is not an alignment.

    The message must name the lengths and offending taxa: 'invalid alignment'
    sends the caller looking through the whole file.
    """
    with pytest.raises(AlignmentError, match="not all the same length") as exc:
        validate({"a": "ACGT", "b": "ACG", "c": "ACGT", "d": "ACGT"})
    assert "3 sites" in str(exc.value) and "4 sites" in str(exc.value)
    assert "b" in str(exc.value)


def test_fewer_than_four_taxa_is_refused():
    """Three taxa have exactly one unrooted topology, so there is nothing to infer."""
    with pytest.raises(AlignmentError, match="at least 4 sequences"):
        validate({"a": "ACGT", "b": "ACGT", "c": "ACGT"})


def test_newick_hostile_names_are_refused():
    """Commas and parentheses are Newick SYNTAX; a tip named with them corrupts
    the output tree rather than failing."""
    with pytest.raises(AlignmentError, match="unsafe in Newick"):
        validate({"a,b": "ACGT", "b": "ACGT", "c": "ACGT", "d": "ACGT"})


def test_protein_input_is_refused_rather_than_read_as_dna():
    with pytest.raises(AlignmentError, match="nucleotide alignments only"):
        validate({"a": "MKVLE", "b": "MKVLE", "c": "MKVLE", "d": "MKVLE"})


def test_parsimony_informative_ignores_singletons():
    """A site where one taxon differs carries no topological information.

    Counting it would overstate the evidence, which is the number the
    thin_information advisory is built on.
    """
    singleton = {"a": "A", "b": "A", "c": "A", "d": "C"}
    informative = {"a": "A", "b": "A", "c": "C", "d": "C"}
    assert parsimony_informative(singleton) == 0
    assert parsimony_informative(informative) == 1


def test_parsimony_informative_ignores_gaps_and_ambiguity():
    assert parsimony_informative({"a": "-", "b": "-", "c": "N", "d": "N"}) == 0


def test_constant_sites_are_not_informative():
    assert (
        parsimony_informative({"a": "AAAA", "b": "AAAA", "c": "AAAA", "d": "AAAA"}) == 0
    )


def test_identical_sequences_are_reported():
    groups = duplicate_groups({"a": "ACGT", "b": "ACGT", "c": "TTTT", "d": "GGGG"})
    assert groups == [["a", "b"]]


def test_duplicate_sequences_warn_through_the_server():
    seqs = dict(FOUR)
    seqs["e"] = seqs["a"]
    result = infer_tree(fasta=_fasta(seqs), model="JC", replicates=20, seed=1)
    assert "duplicate_sequences" in {w["code"] for w in result["warnings"]}


def test_summary_reports_length_and_information_separately():
    """Length and evidence are different quantities and must not be conflated."""
    stats = summarise(
        {"a": "AAAAAAAC", "b": "AAAAAAAC", "c": "AAAAAAAG", "d": "AAAAAAAG"}
    )
    assert stats.n_sites == 8
    assert stats.n_parsimony_informative == 1


def test_gap_fraction_is_measured():
    stats = summarise({"a": "AC--", "b": "AC--", "c": "GT--", "d": "GT--"})
    assert stats.fraction_gaps == 0.5
