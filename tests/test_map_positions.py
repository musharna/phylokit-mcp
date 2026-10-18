"""map_positions: residue position <-> alignment column.

The expected answers are planted by construction (each residue is tracked as
gaps are inserted), not recomputed by the code under test. The last test runs
the whole chain, real MAFFT included, on the agent benchmark's homology
questions, whose key comes from an indel simulator that never calls MAFFT.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from phylokit_mcp.alignment import AlignmentError
from phylokit_mcp.server import map_positions

ALN = ">a\nAC-GT\n>b\nA-CGT\n>c\n--CG-\n"


def test_a_position_maps_to_its_column_and_across_rows():
    r = map_positions(aligned=True, sequence="a", position=3, fasta=ALN)
    assert r["column"] == 4
    assert r["rows"]["a"] == {"position": 3, "residue": "G"}
    assert r["rows"]["b"] == {"position": 3, "residue": "G"}
    assert r["rows"]["c"] == {"position": 2, "residue": "G"}
    # Residue 2 of a sits over a gap in b: no position, and the residues
    # either side of the gap.
    r = map_positions(aligned=True, sequence="a", position=2, fasta=ALN)
    assert r["column"] == 2
    assert r["rows"]["b"] == {"position": None, "residue": "-", "between": [1, 2]}
    assert r["rows"]["c"] == {"position": None, "residue": "-", "between": [None, 1]}


def test_a_column_maps_to_every_row():
    r = map_positions(aligned=True, column=3, fasta=ALN)
    assert r["rows"]["a"] == {"position": None, "residue": "-", "between": [2, 3]}
    assert r["rows"]["b"] == {"position": 2, "residue": "C"}
    assert r["rows"]["c"] == {"position": 1, "residue": "C"}
    assert r["n_columns"] == 5


def test_planted_random_alignments():
    """Rows built by inserting gaps into labelled residues: the column of every
    residue is known before the tool sees the alignment."""
    rng = random.Random(7)
    for _ in range(30):
        n_cols = rng.randint(20, 80)
        rows, where = {}, {}
        for nm in ("s1", "s2", "s3", "s4"):
            cols = sorted(rng.sample(range(n_cols), rng.randint(5, n_cols)))
            row = ["-"] * n_cols
            for c in cols:
                row[c] = rng.choice("ACGT")
            rows[nm] = "".join(row)
            where[nm] = {c: i + 1 for i, c in enumerate(cols)}  # column -> position
        fasta = "".join(f">{k}\n{v}\n" for k, v in rows.items())
        src = "s1"
        col, pos = rng.choice(sorted(where[src].items()))
        r = map_positions(aligned=True, sequence=src, position=pos, fasta=fasta)
        assert r["column"] == col + 1
        for nm in rows:
            assert r["rows"][nm]["position"] == where[nm].get(col)


def test_bad_queries_are_refused_and_good_ones_are_not():
    assert (
        map_positions(aligned=True, sequence="a", position=4, fasta=ALN)["column"] == 5
    )
    with pytest.raises(AlignmentError, match="position must be 1-4"):
        map_positions(aligned=True, sequence="a", position=5, fasta=ALN)
    with pytest.raises(AlignmentError, match="position must be 1-4"):
        map_positions(aligned=True, sequence="a", position=0, fasta=ALN)
    with pytest.raises(AlignmentError, match="No sequence named 'z'"):
        map_positions(aligned=True, sequence="z", position=1, fasta=ALN)
    with pytest.raises(AlignmentError, match="exactly one of position"):
        map_positions(aligned=True, sequence="a", position=1, column=1, fasta=ALN)
    with pytest.raises(AlignmentError, match="column must be 1-5"):
        map_positions(aligned=True, column=6, fasta=ALN)
    with pytest.raises(AlignmentError, match="not all the same length"):
        map_positions(aligned=True, column=1, fasta=">a\nACG\n>b\nAC\n")
    # '?' could be a residue or padding; the count after it depends on which.
    with pytest.raises(AlignmentError, match=r"contains \['\?'\]"):
        map_positions(aligned=True, column=1, fasta=">a\nA?GT\n>b\nACGT\n")


def test_unaligned_input_is_aligned_by_mafft_first():
    seqs = ">a\nACGTACGTTTGACCA\n>b\nACGTACGTGACCA\n"
    r = map_positions(aligned=False, sequence="a", position=13, fasta=seqs)
    assert r["aligned_by"]["name"] == "MAFFT"
    assert r["rows"]["a"]["residue"] == "C"
    assert r["rows"]["b"] == {"position": 11, "residue": "C"}
    assert [w["code"] for w in r["warnings"]] == ["column_is_an_inference"]
    # Declared, not sniffed: the same unequal rows called aligned are refused.
    with pytest.raises(AlignmentError, match="not all the same length"):
        map_positions(aligned=True, sequence="a", position=13, fasta=seqs)


QUESTIONS = Path(__file__).parent.parent / "benchmarks" / "agent_qa" / "questions.json"


def _homology_questions():
    key = json.loads(QUESTIONS.read_text())["questions"]
    return [q for q in key if q["category"] == "homologous_position"]


@pytest.mark.parametrize("q", _homology_questions(), ids=lambda q: q["id"])
def test_the_benchmark_homology_key_through_real_mafft(q, tmp_path):
    import re

    m = re.search(r"which residue of (\S+) is in the same column as residue (\d+) of (\S+)\?",
                  q["question"])  # fmt: skip
    assert m, q["question"]
    target, pos, source = m.group(1), int(m.group(2)), m.group(3)
    ((name, text),) = q["files"].items()
    path = tmp_path / name
    path.write_text(text)
    moltype = q["truth"]["moltype"]
    r = map_positions(
        aligned=False, sequence=source, position=pos, fasta_path=str(path),
        sequence_type=moltype,
    )  # fmt: skip
    assert r["rows"][target]["position"] == q["truth"]["position"]
