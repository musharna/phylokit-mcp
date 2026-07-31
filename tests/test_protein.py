"""Protein alignments, end to end through the real engine.

piqtree has supported protein models all along (`available_models('protein')`);
this server refused them at the door. The refusal is now a DECLARATION rather
than a limitation: `sequence_type="protein"` opts in, and nothing sniffs.

Why nothing sniffs: an alignment of only A/C/G/T is a perfectly valid protein
alignment (alanine, cysteine, glycine, threonine). A guesser would read that as
DNA and fit a nucleotide model to protein data, returning a tree, a likelihood
and support values — all wrong, none complaining.
"""

import pytest

from phylokit_mcp.alignment import AlignmentError
from phylokit_mcp.bootstrap import bootstrap_support
from phylokit_mcp.inference import build_ml_tree, select_model

# Four taxa, two clades: a/b share K and I, c/d share R and V.
PROTEIN = {
    "a": "MKVLEAIGKTLSDAWKMKVLEAIGKTLSDAWK",
    "b": "MKVLEAIGKTLSDAWKMKVLEAIGKTLSDAWK".replace("W", "Y", 1),
    "c": "MRVLEAVGKTLSDAWKMRVLEAVGKTLSDAWK",
    "d": "MRVLEAVGKTLSDAWKMRVLEAVGKTLSDAWK".replace("K", "R", 1),
}


def test_a_protein_tree_is_built_with_a_protein_model():
    """The end-to-end check. A nucleotide model on protein data fails in the
    engine, so getting a tree back at all proves the moltype reached piqtree."""
    tree = build_ml_tree(PROTEIN, "LG", seed=1, moltype="protein")
    newick = tree.get_newick(with_distances=True)
    for taxon in PROTEIN:
        assert taxon in newick


def test_a_nucleotide_model_on_protein_data_fails_rather_than_silently_fitting():
    """The negative control for the test above.

    If this ever stops raising, the previous test proves much less: it would no
    longer distinguish "the protein moltype reached the engine" from "any model
    happens to work".
    """
    with pytest.raises(Exception):  # noqa: B017 - engine-specific type
        build_ml_tree(PROTEIN, "GTR+G", seed=1, moltype="protein")


def test_bootstrap_resamples_protein_columns():
    tree = build_ml_tree(PROTEIN, "LG", seed=1, moltype="protein")
    support = bootstrap_support(
        PROTEIN, tree, "LG", replicates=20, seed=1, moltype="protein"
    )
    assert support.clades, "no clades scored"
    for clade in support.clades:
        assert 0.0 <= clade.support <= 1.0


def test_model_selection_ranks_within_the_protein_model_set():
    """Nucleotide and protein models are not comparable, so the ranking must
    come from the protein set — GTR winning a protein alignment would mean the
    moltype never reached model_finder."""
    out = select_model(PROTEIN, criterion="AIC", seed=1, top_n=5, moltype="protein")
    assert out["n_models_compared"] > 0
    names = [m["model"] for m in out["ranking"]]
    nucleotide_only = {"GTR", "HKY", "JC", "K80", "F81", "TN"}
    assert not ({n.split("+")[0] for n in names} & nucleotide_only), (
        f"nucleotide models ranked for a protein alignment: {names}"
    )


def test_protein_without_declaring_it_is_refused_not_guessed():
    from phylokit_mcp.alignment import validate

    with pytest.raises(AlignmentError, match="sequence_type='protein'"):
        validate(PROTEIN)
