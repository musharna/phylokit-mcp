"""Shared fixtures. No mocked engine anywhere in this suite — every test runs
real IQ-TREE, because the failures worth catching are in what IQ-TREE actually
returns, not in what a stub would be written to return.
"""

from __future__ import annotations

import pytest

# A 7-taxon tree with two clean cherries and a nested group. Balanced enough
# that the true answer is unambiguous, small enough that a 100-replicate
# bootstrap finishes in seconds.
TRUE_NEWICK = (
    "(((A:0.05,B:0.05):0.05,(C:0.05,D:0.05):0.05):0.05,(E:0.05,F:0.05):0.05,G:0.1);"
)

# Clades of TRUE_NEWICK, as canonical splits anchored on taxon "A".
TRUE_CLADES = {"C|D", "E|F", "E|F|G", "C|D|E|F|G"}


def _simulate(length: int, seed: int = 42, model: str = "JC") -> dict:
    from phylokit_mcp.server import simulate_alignment

    return simulate_alignment(newick=TRUE_NEWICK, model=model, length=length, seed=seed)


@pytest.fixture(scope="session")
def easy_alignment() -> dict:
    """Long enough that inference recovers the true tree with full support."""
    return _simulate(length=300)


@pytest.fixture(scope="session")
def hard_alignment() -> dict:
    """Short enough that maximum likelihood gets the topology WRONG.

    This fixture is the point of the suite. At this length the ML tree contains
    a clade that is not in the true tree and misses one that is, while looking
    exactly as authoritative as the easy case. Any test asserting that support
    values are informative needs a case where the topology is actually wrong;
    otherwise it only ever measures the easy regime.
    """
    return _simulate(length=60)


@pytest.fixture(scope="session")
def fasta_easy(easy_alignment: dict) -> str:
    return easy_alignment["fasta"]


@pytest.fixture(scope="session")
def fasta_hard(hard_alignment: dict) -> str:
    return hard_alignment["fasta"]
