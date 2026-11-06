"""map_positions and _check_rows: the messages, keys and checks a caller reads.

test_map_positions.py checks the mapping itself. These pin what it returns or
says around the mapping: each refusal's full text, which check refuses a row,
and the keys of `aligned_by` and of the MAFFT warning.
"""

from __future__ import annotations

import pytest

from phylokit_mcp.alignment import AlignmentError
from phylokit_mcp.inputs import InputFileError
from phylokit_mcp.positions import _check_rows
from phylokit_mcp.server import map_positions

ALN = ">a\nAC-GT\n>b\nA-CGT\n"

AMBIGUOUS = (
    "Row 'a' contains ['?']. Whether those hold a residue position or are "
    "padding changes every position after them, so this tool will not guess: "
    "replace them with '-' (padding) or 'N'/'X' (an unknown residue)."
)


def test_check_rows_counts_columns_and_names_each_refusal():
    assert _check_rows({"a": "A", "b": "C"}) == 1
    assert _check_rows({"a": "AC-GT", "b": "A-CGT"}) == 5
    with pytest.raises(AlignmentError) as exc:
        _check_rows({"a": "", "b": ""})
    assert str(exc.value) == "Rows are empty (zero columns)."
    with pytest.raises(AlignmentError) as exc:
        _check_rows({"a": "ACG", "b": "AC"})
    assert str(exc.value) == (
        "Rows are not all the same length, so this is not an alignment. If the "
        "sequences are unaligned, pass aligned=false and they are aligned with "
        "MAFFT first."
    )
    with pytest.raises(AlignmentError) as exc:
        _check_rows({"a": "A?GT", "b": "ACGT"})
    assert str(exc.value) == AMBIGUOUS


def test_an_unknown_sequence_type_lists_the_valid_ones():
    assert map_positions(aligned=True, column=1, fasta=ALN)["column"] == 1
    with pytest.raises(AlignmentError) as exc:
        map_positions(aligned=True, column=1, fasta=ALN, sequence_type="rna")
    assert str(exc.value) == "Unknown sequence_type 'rna'. Valid: ['dna', 'protein']."


def test_missing_and_doubled_input_name_both_arguments(tmp_path):
    path = tmp_path / "aln.fasta"
    path.write_text(ALN)
    assert map_positions(aligned=True, column=1, fasta_path=str(path))["column"] == 1
    with pytest.raises(InputFileError) as exc:
        map_positions(aligned=True, column=1)
    assert str(exc.value) == (
        "Give fasta (the text inline) or fasta_path (the path of a file)."
    )
    with pytest.raises(InputFileError) as exc:
        map_positions(aligned=True, column=1, fasta=ALN, fasta_path=str(path))
    assert str(exc.value) == (
        "Give fasta or fasta_path, not both: with both there is no way to know "
        "which one you meant."
    )


def test_characters_are_checked_for_the_declared_sequence_type():
    # Nine characters DNA rejects; the message names the first eight, sorted.
    bad = ">a\nACGTEFIJLOPQZ\n>b\nACGTACGTACGTA\n"
    with pytest.raises(AlignmentError) as exc:
        map_positions(aligned=True, column=1, fasta=bad)
    assert str(exc.value) == (
        "Unrecognised characters for dna: ['E', 'F', 'I', 'J', 'L', 'O', 'P', 'Q']."
    )
    # The same letters are residues in a protein.
    protein = ">a\nMEFIL\n>b\nMQFIL\n"
    r = map_positions(aligned=True, column=2, fasta=protein, sequence_type="protein")
    assert r["rows"]["b"]["residue"] == "Q"
    with pytest.raises(AlignmentError) as exc:
        map_positions(aligned=True, column=2, fasta=protein)
    # 'M' is not listed: it is the IUPAC code for A or C.
    assert (
        str(exc.value) == "Unrecognised characters for dna: ['E', 'F', 'I', 'L', 'Q']."
    )


def test_a_dot_gets_the_padding_question_not_a_character_refusal():
    # cogent3 rejects '.', so without its exemption the row would be refused
    # as an unknown character instead of asked about as possible padding.
    assert map_positions(aligned=True, column=1, fasta=">a\nA-GT\n>b\nACGT\n")
    with pytest.raises(AlignmentError) as exc:
        map_positions(aligned=True, column=1, fasta=">a\nA.GT\n>b\nACGT\n")
    assert str(exc.value) == AMBIGUOUS.replace("['?']", "['.']")


def test_aligned_by_is_none_for_given_alignments_and_names_mafft_otherwise():
    given = map_positions(aligned=True, sequence="a", position=3, fasta=ALN)
    assert given["aligned_by"] is None
    assert given["warnings"] == []

    seqs = ">a\nACGTACGTTTGACCA\n>b\nACGTACGTGACCA\n"
    r = map_positions(aligned=False, sequence="a", position=13, fasta=seqs)
    assert sorted(r["aligned_by"]) == ["name", "options", "version"]
    assert r["aligned_by"]["name"] == "MAFFT"
    assert r["aligned_by"]["version"]
    assert r["aligned_by"]["options"]
    assert r["warnings"] == [
        {
            "code": "column_is_an_inference",
            "message": (
                "The column comes from this MAFFT alignment. Another aligner or "
                "settings can place a residue elsewhere, most often inside gappy "
                "stretches: check `window` before relying on it."
            ),
        }
    ]
