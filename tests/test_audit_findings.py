"""Two failures that were real but unreadable.

Both come from the six-audit sweep of this server. Neither returned a wrong
tree — the server was already honest about topology — but both left the caller
unable to act on what happened.
"""

import pytest

from phylokit_mcp.alignment import AlignmentError, require_phylogenetic_signal
from phylokit_mcp.diagnostics import saturated_branch_lengths

S = "ACGTACGTACGTACGTACGT" * 3


def _seqs(*rows):
    return {chr(ord("A") + i): r for i, r in enumerate(rows)}


# --------------------------------------------------------------------------
# No phylogenetic signal
# --------------------------------------------------------------------------


def test_identical_sequences_are_diagnosed_not_left_to_iqtree():
    """Before: 'IQ-TREE output is malformed, likelihood not found.'

    That message is piqtree's, and it describes a PARSING failure — so a caller
    reads it as this server being broken and retries. The real cause is that the
    data cannot support any tree.
    """
    with pytest.raises(AlignmentError) as exc:
        require_phylogenetic_signal(_seqs(S, S, S, S))
    msg = str(exc.value)
    assert "no phylogenetic signal" in msg
    assert "All 4 sequences are identical" in msg
    # The caller must learn retrying is pointless; that is the whole difference
    # from the upstream message.
    assert "will not change it" in msg
    assert "malformed" not in msg


def test_varying_but_uninformative_sites_are_also_caught():
    """Variation alone is not signal. A site where one taxon differs is a
    singleton: it cannot separate one grouping from another."""
    with pytest.raises(AlignmentError) as exc:
        require_phylogenetic_signal(_seqs(S, S, S, "T" + S[1:]))
    assert "no phylogenetic signal" in str(exc.value)
    assert "identical" not in str(exc.value)  # the other branch of the message


def test_a_real_alignment_still_validates():
    """Positive control. A signal check that rejects real data is worse than none."""
    a = "ACGTACGTACGTACGTACGT"
    require_phylogenetic_signal(_seqs(a + "AAAA", a + "AAAA", a + "GGGG", a + "GGGG"))


# --------------------------------------------------------------------------
# Saturated branch lengths
# --------------------------------------------------------------------------


def test_a_branch_at_the_ceiling_is_flagged():
    """Measured on 300bp of random sequence: lengths of 9.9999989, unflagged."""
    w = saturated_branch_lengths(
        "(A:4.2010685031,((B:9.9999989755,E:9.9999989678):9.999998966,C:2e-06):4.78,D:2.41);"
    )
    assert w is not None
    assert w["code"] == "saturated_branch_lengths"
    assert w["n_saturated"] == 3
    assert w["max_branch_length"] == pytest.approx(9.9999989755)
    # It must say the number is a floor, not a measurement.
    assert "at least" in w["message"]


def test_ordinary_branch_lengths_are_not_flagged():
    """Positive control, and it is not hypothetical: a real tree from this
    server carries lengths like 0.03 and 0.108."""
    assert (
        saturated_branch_lengths(
            "(A:1e-06,B:0.0326350929,(C:1e-06,D:0.0326684687):0.1079926343);"
        )
        is None
    )


def test_scientific_notation_is_parsed():
    """`1e-06` appears in real output from this server; a regex that missed it
    would silently under-count branches and could miss a saturated one."""
    assert saturated_branch_lengths("(A:1e-06,B:9.99e0,C:0.1,D:0.2);") is not None


# --------------------------------------------------------------------------
# Provenance and units on the result itself
# --------------------------------------------------------------------------


def _infer(fasta: str) -> dict:
    import asyncio
    import json

    from phylokit_mcp.server import build_server

    mcp = build_server()

    r = asyncio.run(mcp.call_tool("infer_tree", {"fasta": fasta, "seed": 1}))
    return json.loads(
        next(c.text for c in (r.content or []) if getattr(c, "text", None))
    )


REAL = (
    ">A\nACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n"
    ">B\nACGTTGCAACGTTGCAACGTTGCTACGTTGCA\n"
    ">C\nTCGTTGCAACGTAGCAACGTTGCAACGTTGGA\n"
    ">D\nTCGTTGCAACGTAGCAACGTTGCAACGTTGGT\n"
)


def test_the_tree_carries_the_engine_that_built_it():
    out = _infer(REAL)
    assert out["engine"]["version"], "no engine version on the result"
    assert "piqtree" in out["engine"]["name"]
    # Must agree with the introspection tool; two sources that can disagree are
    # worse than one.
    import asyncio
    import json

    from phylokit_mcp.server import build_server

    mcp = build_server()

    cap = asyncio.run(mcp.call_tool("capabilities", {}))
    caps = json.loads(
        next(c.text for c in (cap.content or []) if getattr(c, "text", None))
    )
    assert out["engine"]["version"] == caps["engine_version"]


def test_branch_lengths_state_their_units():
    """Newick carries bare numbers; 'substitutions per site' is not guessable
    from the string, and reading them as time or percent gives a wrong answer."""
    assert _infer(REAL)["branch_length_units"] == "substitutions per site"
