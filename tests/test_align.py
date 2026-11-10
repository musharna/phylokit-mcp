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
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from phylokit_mcp import engine, msa
from phylokit_mcp.alignment import (
    AlignmentError,
    parse_fasta,
    require_phylogenetic_signal,
    validate,
)
from phylokit_mcp.inputs import InputFileError
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

    # A hung aligner is a hung server, so the wall-clock cap is a refusal. The
    # short cap is scoped so the real binary below runs under the shipped one.
    hang = tmp_path / "hang"
    hang.mkdir()
    with monkeypatch.context() as m:
        m.setenv("PATH", _stand_in(hang, "exec /bin/sleep 5\n"))
        m.setattr(msa, "TIMEOUT_SECONDS", 0.5)
        with pytest.raises(AlignmentError, match="did not finish within 0.5 s"):
            align_sequences(fasta)

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

    # The short cap is scoped: left in force, it timed the real `mafft
    # --version` below too, which takes ~0.2 s and failed on a busy machine.
    hang = tmp_path / "hang"
    hang.mkdir()
    with monkeypatch.context() as m:
        m.setenv("PATH", _stand_in(hang, "exec /bin/sleep 5\n"))
        m.setattr(msa, "VERSION_TIMEOUT_SECONDS", 0.5)
        caps = capabilities()
    assert caps["aligner_version"] is None
    assert "did not answer within 0.5 s" in caps["aligner_error"]

    # Absent is not an error, and neither is working.
    monkeypatch.setenv("PATH", str(tmp_path / "nowhere"))
    assert capabilities()["aligner_error"] is None
    monkeypatch.setenv("PATH", real_path)
    caps = capabilities()
    assert caps["aligner_version"].startswith("v") and caps["aligner_error"] is None


def test_out_path_writes_the_alignment_infer_tree_reads_by_path(fasta_easy, tmp_path):
    """The handoff without retyping: align to a file, infer from that file.

    A client with no file-writing tool of its own (the field benchmark's server
    arm had only Read) used to copy ~15 x 962 aligned characters from
    `align_sequences` into `infer_tree`, mis-copied one row each time, and was
    refused as ragged 13 times in one run. With `out_path` nothing is copied.
    """
    seqs = parse_fasta(fasta_easy)
    ragged = {nm: s[:30] + s[36:] if nm == "A" else s for nm, s in seqs.items()}
    out = tmp_path / "aligned.fasta"

    res = align_sequences(_fasta(ragged), out_path=str(out))
    assert res["written_to"] == str(out.resolve())
    assert res["fasta"] is None, "the rows live in the file, not in the reply"
    rows = parse_fasta(out.read_text())
    assert set(rows) == set(seqs)
    assert len({len(s) for s in rows.values()}) == 1
    assert {nm: s.replace("-", "") for nm, s in rows.items()} == ragged
    tree = infer_tree(fasta_path=res["written_to"], model="JC", replicates=20)
    assert compare_trees(tree["newick"], TRUE_NEWICK)["identical_topology"] is True

    # Never over an existing file, and the refusal leaves it as it was.
    before = out.read_bytes()
    with pytest.raises(ValueError, match="out_path: .* already exists"):
        align_sequences(_fasta(ragged), out_path=str(out))
    assert out.read_bytes() == before
    with pytest.raises(ValueError, match="out_path: no such directory"):
        align_sequences(_fasta(ragged), out_path=str(tmp_path / "nope" / "a.fasta"))
    # Without out_path the reply is unchanged: inline FASTA, nothing written.
    inline = align_sequences(_fasta(ragged))
    assert inline["written_to"] is None
    assert parse_fasta(inline["fasta"]) == rows


def test_each_refusal_before_mafft_names_exactly_what_to_fix():
    """The whole message, not a fragment of it: the callers matched on a
    substring, so 30 rewordings of these refusals survived mutation, among
    them a protein hint shown for protein and a ninth character listed."""
    from phylokit_mcp.msa import validate_unaligned

    def refusal(*args):
        with pytest.raises(AlignmentError) as exc:
            validate_unaligned(*args)
        return str(exc.value)

    ok = {"a": "ACGTACGT", "b": "ACGTTACGT"}
    # The declaration defaults to DNA, and the same pair passes as protein.
    assert validate_unaligned(ok) is None
    assert validate_unaligned(ok, "protein") is None

    # Unreachable from the tools, which check sequence_type first; msa.align
    # is called directly too, so the guard is its own.
    assert refusal(ok, "rna") == (
        "Unknown sequence_type 'rna'. Valid: ['dna', 'protein']."
    )
    assert refusal({"a,b": "ACGT", "c": "ACGA"}) == (
        "Sequence name 'a,b' contains characters that are unsafe in Newick "
        "(the tree format uses ,:;() as syntax). Use letters, digits, "
        "underscore, dot or hyphen."
    )
    assert refusal({"a": "ACGT", "b": ""}) == "Sequence 'b' is empty."
    assert refusal({"a": "AC-GT", "b": "ACGA?"}) == (
        "Input contains gap or missing-data characters ['-', '?']. "
        "align_sequences takes UNALIGNED sequences; remove the gaps first. "
        "If the sequences are already aligned, pass them to infer_tree as "
        "they are."
    )
    # Ten letters DNA rejects, so the list is cut at its eighth; under DNA the
    # protein remedy is named.
    assert len(msa.unrecognised_characters(set("EFIJLOPQXZ"), "dna")) == 10
    assert refusal({"a": "EFIJLOPQXZ", "b": "ACGT"}) == (
        "Unrecognised characters for dna: ['E', 'F', 'I', 'J', 'L', 'O', 'P', "
        "'Q']. If these are protein sequences, pass sequence_type='protein': "
        "the molecule type is declared, never guessed."
    )
    # Under protein that remedy would be wrong, so nothing follows the list.
    assert refusal({"a": "MKV1", "b": "MKVL"}, "protein") == (
        "Unrecognised characters for protein: ['1']."
    )


def test_the_refusals_built_only_from_f_strings_are_asserted_whole_too():
    """mutmut 3.8.1 leaves f-string text alone, so no survivor list can show
    these messages drifting; each is checked here, with its cap's edge."""
    from phylokit_mcp.msa import validate_unaligned

    def refusal(seqs, moltype="dna"):
        with pytest.raises(AlignmentError) as exc:
            validate_unaligned(seqs, moltype)
        return str(exc.value)

    assert refusal({"a": "ACGT"}) == "Need at least 2 sequences to align, got 1."
    assert validate_unaligned({"a": "ACGT", "b": "ACGA"}) is None

    at_cap = {f"s{i}": "ACGT" for i in range(200)}
    assert validate_unaligned(at_cap) is None
    assert refusal({**at_cap, "s200": "ACGT"}) == (
        "Got 201 sequences; this server caps at 200."
    )
    assert validate_unaligned({"a": "A" * 100_000, "b": "ACGT"}) is None
    assert refusal({"a": "A" * 100_001, "b": "ACGT"}) == (
        "Sequence 'a' is 100001 residues; this server caps at 100000."
    )

    stop = (
        " '*' is a stop codon, and the engine has no state for it: trim stops "
        "(or replace them with X) before aligning or inferring."
    )
    assert refusal({"a": "MKV*", "b": "MKVL"}, "protein") == (
        "Unrecognised characters for protein: ['*']." + stop
    )
    # Under DNA the stop-codon remedy follows the protein one.
    assert refusal({"a": "ACG*", "b": "ACGT"}) == (
        "Unrecognised characters for dna: ['*']. If these are protein sequences, "
        "pass sequence_type='protein': the molecule type is declared, never "
        "guessed." + stop
    )


# --------------------------------------------------------------------------
# What MAFFT is run with, and what its failures report, compared whole
# --------------------------------------------------------------------------


def _eight() -> dict[str, str]:
    """Eight distinct sequences: more names than an error message lists."""
    anc = _ancestor()
    return {f"s{i}": _substitute(anc, [i, 30 + i, 70 + 2 * i]) for i in range(8)}


def _wrapping(tmp_path, name: str, script: str, path: str) -> str:
    """A stand-in that runs the REAL binary and edits its output with sed.

    PATH holds only the stand-in, so both are looked up on the real `path`,
    and the real MAFFT (a shell script) runs with that `path` restored."""
    real, sed = shutil.which("mafft", path=path), shutil.which("sed", path=path)
    assert real and sed, "these tests need the real MAFFT, and sed, on PATH"
    d = tmp_path / name
    d.mkdir()
    return _stand_in(d, f'PATH="{path}" "{real}" "$@" | "{sed}" {script}\n')


def test_mafft_runs_in_its_own_directory_and_never_reads_the_servers_stdin(
    monkeypatch, tmp_path
):
    """The server's stdin is the MCP stream. A child that inherited it could
    read protocol bytes; so MAFFT, and `mafft --version`, get /dev/null. The
    input goes to MAFFT as `input.fasta` in a private `phylokit-mafft-` dir
    that is also its working directory, and is gone afterwards."""
    real = shutil.which("mafft")
    assert real
    log = tmp_path / "log"
    log.mkdir()
    shim = tmp_path / "shim"
    shim.mkdir()
    _stand_in(
        shim,
        f'if [ "$1" = --version ]; then cat > "{log}/stdin_version"; '
        f'exec "{real}" "$@"; fi\n'
        f'cat > "{log}/stdin_align"\npwd > "{log}/cwd"\n'
        f'for a in "$@"; do last="$a"; done\nprintf %s "$last" > "{log}/input"\n'
        f'exec "{real}" "$@"\n',
    )
    # Positive control: the stand-in does record a stdin it is given.
    subprocess.run([str(shim / "mafft"), "--version"], input=b"probe", check=False)
    assert (log / "stdin_version").read_bytes() == b"probe"

    monkeypatch.setenv("PATH", str(shim) + os.pathsep + os.environ["PATH"])
    fasta = ">a\nACGTACGTACGGTT\n>b\nACGTACGTTTGG\n"
    r, w = os.pipe()
    os.write(w, b"PROTOCOL BYTES\n")
    os.close(w)
    saved = os.dup(0)
    os.dup2(r, 0)
    try:
        out = align_sequences(fasta)
    finally:
        os.dup2(saved, 0)
        os.close(saved)
        os.close(r)

    assert (log / "stdin_align").read_bytes() == b""
    assert (log / "stdin_version").read_bytes() == b""
    infile = Path((log / "input").read_text())
    assert infile.name == "input.fasta"
    assert infile.parent.name.startswith("phylokit-mafft-")
    assert Path((log / "cwd").read_text().strip()).resolve() == infile.parent.resolve()
    assert not infile.parent.exists()
    assert parse_fasta(out["fasta"])["b"].replace("-", "") == "ACGTACGTTTGG"


def test_mafft_failures_report_the_last_4000_characters_and_six_names(
    monkeypatch, tmp_path
):
    seqs = _eight()
    names = list(seqs)
    tail = "".join(str(i % 10) for i in range(5000))
    real_path = os.environ["PATH"]

    crash = tmp_path / "crash"
    crash.mkdir()
    monkeypatch.setenv("PATH", _stand_in(crash, f"printf %s '{tail}' >&2\nexit 1\n"))
    with pytest.raises(msa.MafftFailedError) as exc:
        msa.align(seqs)
    assert str(exc.value).startswith("MAFFT exited 1. Command: [")
    assert str(exc.value).endswith(f"]. stderr:\n{tail[-4000:]}")

    garbled = tmp_path / "garbled"
    garbled.mkdir()
    monkeypatch.setenv(
        "PATH",
        _stand_in(garbled, f"printf 'not fasta\\n'\nprintf %s '{tail}' >&2\nexit 0\n"),
    )
    with pytest.raises(AlignmentError) as why:
        parse_fasta("not fasta\n")
    with pytest.raises(msa.MafftFailedError) as exc:
        msa.align(seqs)
    assert str(exc.value) == (
        f"MAFFT exited 0 but its output is not FASTA ({why.value}). "
        f"stderr:\n{tail[-4000:]}"
    )

    monkeypatch.setenv("PATH", _wrapping(tmp_path, "renamed", "'s/^>/>x/'", real_path))
    with pytest.raises(msa.MafftFailedError) as exc:
        msa.align(seqs)
    assert str(exc.value) == (
        f"MAFFT returned names {['x' + n for n in names][:6]} for input names "
        f"{names[:6]}; names or order were not preserved."
    )

    monkeypatch.setenv(
        "PATH", _wrapping(tmp_path, "altered", "'2y/ACGT/CGTA/'", real_path)
    )
    with pytest.raises(msa.MafftFailedError) as exc:
        msa.align(seqs)
    assert str(exc.value) == (
        "MAFFT's row for 's0', with gaps removed, is not the input sequence: "
        "residues were altered or dropped."
    )

    monkeypatch.setenv("PATH", _wrapping(tmp_path, "ragged", "'$s/$/-/'", real_path))
    with pytest.raises(msa.MafftFailedError) as exc:
        msa.align(seqs)
    assert str(exc.value) == "MAFFT returned rows of unequal length."

    hang = tmp_path / "hang"
    hang.mkdir()
    with monkeypatch.context() as m:
        m.setenv("PATH", _stand_in(hang, "exec /bin/sleep 5\n"))
        m.setattr(msa, "TIMEOUT_SECONDS", 0.5)
        with pytest.raises(AlignmentError) as exc:
            msa.align(seqs)
    assert str(exc.value) == (
        "MAFFT did not finish within 0.5 s on 8 sequences. This server is "
        "synchronous; align an input this large outside it."
    )

    # Positive control: the real binary, on the same input and the default
    # molecule type, passes every check above.
    monkeypatch.setenv("PATH", real_path)
    aligned = msa.align(seqs)
    assert list(aligned) == names
    assert {nm: row.replace("-", "") for nm, row in aligned.items()} == seqs


def test_a_version_that_exits_1_is_still_read(monkeypatch, tmp_path):
    """Some builds print the version to stderr and exit 1; that is an answer."""
    real_path = os.environ["PATH"]
    one = tmp_path / "one"
    one.mkdir()
    monkeypatch.setenv(
        "PATH", _stand_in(one, "echo 'v7.999 (2099/01/01)' >&2\nexit 1\n")
    )
    assert msa.mafft_version() == "v7.999 (2099/01/01)"
    monkeypatch.setenv("PATH", real_path)
    assert msa.mafft_version().startswith("v")


def test_align_sequences_refusal_names_its_arguments():
    with pytest.raises(InputFileError) as exc:
        align_sequences()
    assert str(exc.value) == (
        "Give fasta (the text inline) or fasta_path (the path of a file)."
    )
    assert align_sequences(">a\nACGTACGTAC\n>b\nACGTACTAC\n")["fasta"]


def test_align_sequences_advises_on_the_aligned_rows_and_the_thread_state(
    monkeypatch,
):
    dup = ">A\nACGTACGTACGG\n>B\nACGTACGTACGG\n>C\nACGAACGTTCGG\n>D\nTCGAACGTTAGG\n"
    for pinned in (True, False):
        monkeypatch.setattr(engine, "threads_pinned", lambda pinned=pinned: pinned)
        codes = [w["code"] for w in align_sequences(dup)["warnings"]]
        assert "duplicate_sequences" in codes
        assert ("threads_not_pinned" in codes) is not pinned


def test_ready_for_infer_tree_checks_protein_as_protein_and_says_why_not():
    prot = (
        ">A\nMKLVWQERTYAAGG\n>B\nMKLVWQERTFAAGG\n"
        ">C\nMKLAWQDRTFAAGH\n>D\nMRLAWQDKTFAAGH\n"
    )
    out = align_sequences(prot, sequence_type="protein")
    assert out["ready_for_infer_tree"] is True
    assert "infer_tree_would_refuse" not in {w["code"] for w in out["warnings"]}

    same = "".join(f">{t}\nACGTACGTACGTACGT\n" for t in "ABCD")
    out = align_sequences(same)
    aligned = parse_fasta(out["fasta"])
    validate(aligned, "dna")
    with pytest.raises(AlignmentError) as why:
        require_phylogenetic_signal(aligned, "dna")
    assert out["ready_for_infer_tree"] is False
    assert {"code": "infer_tree_would_refuse", "message": str(why.value)} in out[
        "warnings"
    ]

    out = align_sequences(">a\nACGTACGTAC\n>b\nACGTACTAC\n")
    assert out["ready_for_infer_tree"] is False
    assert {
        "code": "too_few_taxa_for_a_tree",
        "message": (
            "2 sequences were aligned, but infer_tree needs at least 4: an "
            "unrooted tree on fewer has only one topology."
        ),
    } in out["warnings"]
