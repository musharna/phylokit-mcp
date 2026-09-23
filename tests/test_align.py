"""align_sequences, against the real MAFFT binary.

Same rule as the rest of the suite: no mocked engine. These tests FAIL when
MAFFT is not on PATH rather than skipping, because a skipped alignment suite is
green on a machine that cannot align. The two tests that put a stand-in `mafft`
on PATH do so to reach failure paths a working MAFFT never takes, and each one
also drives the real binary in the same test.
"""

from __future__ import annotations

import os
import random
import stat

import pytest

from phylokit_mcp import msa
from phylokit_mcp.alignment import AlignmentError, parse_fasta
from phylokit_mcp.server import (
    align_sequences,
    build_server,
    capabilities,
    compare_trees,
    infer_tree,
)

from .conftest import TRUE_NEWICK

DEL_START, DEL_LEN = 40, 12
INS_AT = 80
INSERT = "GATTACAGG"


def _ancestor() -> str:
    """A 120-base ancestor on which the planted indels have ONE correct placement.

    A deletion whose first base equals the base just after it can be slid one
    column with no change in score, and then "the gap is in the right columns"
    has two right answers. The seed is walked until neither indel has that
    freedom, so the expected columns below are the only optimal ones.
    """
    for seed in range(1000):
        rng = random.Random(seed)
        s = "".join(rng.choice("ACGT") for _ in range(120))
        e = DEL_START + DEL_LEN
        deletion_fixed = s[DEL_START] != s[e] and s[DEL_START - 1] != s[e - 1]
        insertion_fixed = INSERT[0] != s[INS_AT] and INSERT[-1] != s[INS_AT - 1]
        if deletion_fixed and insertion_fixed:
            return s
    raise AssertionError("no seed gave an unambiguous indel placement")


def _substitute(s: str, positions: list[int]) -> str:
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return "".join(swap[c] if i in positions else c for i, c in enumerate(s))


def _planted() -> tuple[str, dict[str, str]]:
    anc = _ancestor()
    seqs = {
        "anc": anc,
        "deleted": anc[:DEL_START] + anc[DEL_START + DEL_LEN :],
        "subs1": _substitute(anc, [5, 60, 110]),
        "inserted": anc[:INS_AT] + INSERT + anc[INS_AT:],
        "subs2": _substitute(anc, [20, 95]),
    }
    return anc, seqs


def _fasta(seqs: dict[str, str]) -> str:
    return "".join(f">{nm}\n{s}\n" for nm, s in seqs.items())


def test_planted_indels_come_back_as_gaps_in_the_right_rows_and_columns():
    """The known-answer control: homology was planted, so it can be checked."""
    anc, seqs = _planted()
    out = align_sequences(_fasta(seqs))
    rows = parse_fasta(out["fasta"])

    assert list(rows) == list(seqs)
    width = len(anc) + len(INSERT)
    assert {len(r) for r in rows.values()} == {width}
    assert out["alignment"]["n_sites"] == width
    assert out["input_lengths"] == {"min": len(anc) - DEL_LEN, "max": width}

    # Columns of the ancestor row, by ancestor position.
    col_of = [i for i, c in enumerate(rows["anc"]) if c != "-"]
    assert len(col_of) == len(anc)
    insert_cols = [i for i, c in enumerate(rows["anc"]) if c == "-"]

    # The insertion: the ancestor is gapped in exactly 9 adjacent columns that
    # sit between ancestor positions 79 and 80, and only `inserted` has
    # residues there -- the planted ones.
    assert insert_cols == list(range(col_of[INS_AT - 1] + 1, col_of[INS_AT]))
    assert "".join(rows["inserted"][i] for i in insert_cols) == INSERT
    for nm in ("anc", "deleted", "subs1", "subs2"):
        assert all(rows[nm][i] == "-" for i in insert_cols), nm

    # The deletion: `deleted` is gapped at exactly the 12 ancestor positions
    # that were removed, plus the insertion columns, and nowhere else.
    deleted_cols = {col_of[p] for p in range(DEL_START, DEL_START + DEL_LEN)}
    gaps = {i for i, c in enumerate(rows["deleted"]) if c == "-"}
    assert gaps == deleted_cols | set(insert_cols)

    # Substitution-only rows carry no gap outside the insertion columns, so
    # every one of their residues is in its ancestor's column.
    for nm in ("subs1", "subs2"):
        assert {i for i, c in enumerate(rows[nm]) if c == "-"} == set(insert_cols)

    assert out["engine"]["name"] == "MAFFT"
    assert out["engine"]["version"].startswith("v")
    # Every substitution here is a singleton, so the alignment is correct AND
    # carries no parsimony signal: infer_tree refuses it. This line used to
    # assert True -- `ready` was a taxon count and said yes to output the next
    # tool would reject. It is now infer_tree's own checks, so it says no, and
    # says why.
    assert out["alignment"]["n_parsimony_informative"] == 0
    assert out["ready_for_infer_tree"] is False
    assert [w["code"] for w in out["warnings"]].count("infer_tree_would_refuse") == 1


def test_ragged_sequences_are_refused_by_infer_tree_then_accepted_once_aligned(
    fasta_easy: str,
):
    """End to end: the output feeds the existing tree tool, and gets the tree right."""
    seqs = parse_fasta(fasta_easy)
    cuts = {"A": (30, 6), "C": (120, 9), "F": (200, 3)}
    ragged = {
        nm: s[: cuts[nm][0]] + s[sum(cuts[nm]) :] if nm in cuts else s
        for nm, s in seqs.items()
    }
    assert len({len(s) for s in ragged.values()}) > 1

    with pytest.raises(AlignmentError, match="run align_sequences first"):
        infer_tree(_fasta(ragged), model="JC", replicates=20)

    aligned = align_sequences(_fasta(ragged))
    assert aligned["ready_for_infer_tree"] is True
    tree = infer_tree(aligned["fasta"], model="JC", replicates=20)
    assert tree["alignment"]["n_taxa"] == 7
    assert tree["alignment"]["fraction_gaps"] > 0
    cmp = compare_trees(tree["newick"], TRUE_NEWICK)
    assert cmp["identical_topology"] is True, cmp


def test_a_single_sequence_is_refused_and_two_are_aligned():
    """MAFFT itself exits 0 on one sequence, so this refusal is ours to make."""
    with pytest.raises(AlignmentError, match="at least 2 sequences to align, got 1"):
        align_sequences(">a\nACGTACGTAC\n")

    out = align_sequences(">a\nACGTACGTACGGTTAACC\n>b\nACGTACGTTTAACC\n")
    rows = parse_fasta(out["fasta"])
    assert len(rows["a"]) == len(rows["b"]) == 18
    assert rows["b"].count("-") == 4
    # Two sequences align, but cannot go on to a tree, and the result says so.
    assert out["ready_for_infer_tree"] is False
    assert [w["code"] for w in out["warnings"]].count("too_few_taxa_for_a_tree") == 1


def test_characters_outside_the_declared_alphabet_are_refused():
    protein = ">p1\nMKVLEEFQW\n>p2\nMKVLFQW\n"
    with pytest.raises(AlignmentError, match="Unrecognised characters for dna"):
        align_sequences(protein)
    with pytest.raises(AlignmentError, match="Unrecognised characters for protein"):
        align_sequences(">p1\nMKV1LEEF\n>p2\nMKVLF\n", sequence_type="protein")
    with pytest.raises(AlignmentError, match="Unknown sequence_type 'rna'"):
        align_sequences(protein, sequence_type="rna")

    # The same sequences under the right declaration align.
    out = align_sequences(protein, sequence_type="protein")
    rows = parse_fasta(out["fasta"])
    assert {len(r) for r in rows.values()} == {9}
    assert rows["p2"].replace("-", "") == "MKVLFQW"
    assert out["engine"]["options"][-1] == "--amino"
    assert out["alignment"]["moltype"] == "protein"


def test_gapped_input_is_refused_rather_than_silently_degapped():
    with pytest.raises(AlignmentError, match="takes UNALIGNED sequences"):
        align_sequences(">a\nACGT--ACGT\n>b\nACGTTTACGT\n")
    assert align_sequences(">a\nACGTACGT\n>b\nACGTTTACGT\n")["fasta"]


def test_size_caps_are_the_servers_existing_ones(monkeypatch):
    ok = ">a\nACGTACGTAC\n>b\nACGTACTAC\n"
    monkeypatch.setattr(msa, "MAX_SITES", 9)
    with pytest.raises(AlignmentError, match="'a' is 10 residues; this server caps"):
        align_sequences(ok)
    monkeypatch.setattr(msa, "MAX_SITES", 10)
    monkeypatch.setattr(msa, "MAX_TAXA", 1)
    with pytest.raises(AlignmentError, match="Got 2 sequences; this server caps at 1"):
        align_sequences(ok)
    monkeypatch.setattr(msa, "MAX_TAXA", 2)
    assert align_sequences(ok)["alignment"]["n_taxa"] == 2


def test_a_missing_mafft_is_a_refusal_that_names_the_install(monkeypatch, tmp_path):
    import anyio
    from mcp.server.mcpserver.exceptions import ToolError

    fasta = ">a\nACGTACGTAC\n>b\nACGTACTAC\n"
    mcp = build_server()
    real_path = os.environ["PATH"]

    async def call():
        return await mcp.call_tool("align_sequences", {"fasta": fasta})

    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(msa.MafftUnavailableError, match="MAFFT binary on PATH"):
        align_sequences(fasta)
    # Through the protocol layer the text must survive, not be masked as a crash.
    with pytest.raises(ToolError, match="apt install mafft"):
        anyio.run(call)
    # The tree tools do not need the aligner, and say None rather than failing.
    assert capabilities()["aligner_version"] is None

    monkeypatch.setenv("PATH", real_path)
    assert anyio.run(call) is not None
    assert capabilities()["aligner_version"].startswith("v")
    assert capabilities()["limits"]["min_sequences_to_align"] == 2


def _stand_in(tmp_path, body: str) -> str:
    exe = tmp_path / "mafft"
    exe.write_text("#!/bin/sh\n" + body)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return str(tmp_path)


def test_a_failing_or_lying_aligner_is_an_error_with_its_stderr(monkeypatch, tmp_path):
    fasta = ">a\nACGTACGTAC\n>b\nACGTACTAC\n"
    real_path = os.environ["PATH"]

    crash = tmp_path / "crash"
    crash.mkdir()
    monkeypatch.setenv(
        "PATH", _stand_in(crash, "echo 'segfault in tbfast' >&2\nexit 3\n")
    )
    with pytest.raises(
        msa.MafftFailedError, match=r"exited 3[\s\S]*segfault in tbfast"
    ):
        align_sequences(fasta)

    # Exits 0 with an equal-length "alignment" that has lost a residue of `b`.
    lying = tmp_path / "lying"
    lying.mkdir()
    monkeypatch.setenv(
        "PATH", _stand_in(lying, "printf '>a\\nACGTACGTAC\\n>b\\nACGTAC-AC-\\n'\n")
    )
    with pytest.raises(msa.MafftFailedError, match="'b'.*altered or dropped"):
        align_sequences(fasta)

    reorder = tmp_path / "reorder"
    reorder.mkdir()
    monkeypatch.setenv(
        "PATH", _stand_in(reorder, "printf '>b\\nACGTAC-TAC\\n>a\\nACGTACGTAC\\n'\n")
    )
    with pytest.raises(msa.MafftFailedError, match="names or order"):
        align_sequences(fasta)

    # A hung aligner is a hung server, so the wall-clock cap is a refusal.
    hang = tmp_path / "hang"
    hang.mkdir()
    monkeypatch.setenv("PATH", _stand_in(hang, "exec /bin/sleep 5\n"))
    monkeypatch.setattr(msa, "TIMEOUT_SECONDS", 0.5)
    with pytest.raises(AlignmentError, match="did not finish within 0.5 s"):
        align_sequences(fasta)
    monkeypatch.setattr(msa, "TIMEOUT_SECONDS", 600)

    # Positive control: the same input through the real binary passes all four.
    monkeypatch.setenv("PATH", real_path)
    rows = parse_fasta(align_sequences(fasta)["fasta"])
    assert rows["b"].replace("-", "") == "ACGTACTAC"


def test_sequences_never_reach_a_command_line():
    """The options are fixed and the only caller-derived argument is a temp path."""
    assert msa.mafft_args("dna") == [
        "--auto",
        "--thread",
        "1",
        "--preservecase",
        "--inputorder",
        "--nuc",
    ]
    # Lower case survives, which --preservecase is there for; without it MAFFT
    # lower-cases nucleotides and the degapped-output check would reject it.
    out = align_sequences(">a\nACGTACGTACGGTT\n>b\nacgtacgtttgg\n")
    assert parse_fasta(out["fasta"])["b"].replace("-", "") == "acgtacgtttgg"


def test_a_broken_aligner_is_reported_by_capabilities_not_a_crash(
    monkeypatch, tmp_path
):
    """Installed-but-broken is a third state: not absent, not working.

    `capabilities` is what a caller runs to find out why something fails, so it
    has to survive the fault and name it. `aligner_version: null` alone would
    read as "not installed".
    """
    real_path = os.environ["PATH"]
    silent = tmp_path / "silent"
    silent.mkdir()
    monkeypatch.setenv("PATH", _stand_in(silent, "exit 0\n"))
    caps = capabilities()
    assert caps["aligner_version"] is None
    assert "printed nothing" in caps["aligner_error"]
    assert caps["engine_version"]  # the tree tools' answer is still there

    hang = tmp_path / "hang"
    hang.mkdir()
    monkeypatch.setenv("PATH", _stand_in(hang, "exec /bin/sleep 5\n"))
    monkeypatch.setattr(msa, "VERSION_TIMEOUT_SECONDS", 0.5)
    caps = capabilities()
    assert caps["aligner_version"] is None
    assert "did not answer within 0.5 s" in caps["aligner_error"]

    # Absent is not an error, and neither is working.
    monkeypatch.setenv("PATH", str(tmp_path / "nowhere"))
    assert capabilities()["aligner_error"] is None
    monkeypatch.setenv("PATH", real_path)
    caps = capabilities()
    assert caps["aligner_version"].startswith("v") and caps["aligner_error"] is None
