"""MCP server. Five tools over IQ-TREE.

`engine` is imported FIRST and deliberately: it pins OMP_NUM_THREADS before
piqtree loads, and IQ-TREE reads that variable when its thread pool initialises.
Import order here is load-bearing, not incidental.

Tools are synchronous. IQ-TREE holds process-global state and the bootstrap is
CPU-bound; serialising calls on the event loop keeps concurrent invocations from
interleaving inside the engine.
"""

# isort: off  -- engine MUST come before anything that could pull in piqtree.
from . import engine  # noqa: I001

# isort: on
# stdlib TypedDict is safe here: pydantic only rejects it under Python < 3.12,
# and piqtree ships wheels for 3.12+ only, so this package cannot run on a
# version where the typing_extensions backport would be required.
from typing import Any, NotRequired, TypedDict

# mcp 2.x renamed FastMCP to MCPServer and moved it out of mcp.server.fastmcp,
# which no longer exists. The class is the same one — same decorator, same
# `annotations` and `structured_output` kwargs — so this is a rename, not a
# rewrite. ToolAnnotations below did not move.
from mcp.server.mcpserver import MCPServer
from mcp.types import ToolAnnotations

from . import alignment as aln_mod
from . import diagnostics
from .alignment import AlignmentError, parse_fasta, summarise, validate
from .bootstrap import (
    DEFAULT_REPLICATES,
    MAX_REPLICATES,
    MIN_REPLICATES,
    STRONG_SUPPORT,
    WEAK_SUPPORT,
    bootstrap_support,
)
from .inference import (
    CRITERIA,
    DELTA_INDISTINGUISHABLE,
    build_ml_tree,
    select_model,
    tree_log_likelihood,
)
from .splits import (
    canonical_splits,
    format_split,
    normalised_robinson_foulds,
    robinson_foulds,
)

INSTRUCTIONS = f"""\
Phylogenetic inference over IQ-TREE 2.

The rule this server is built around: **a topology without support is not a
result.** `infer_tree` always runs a bootstrap and always returns per-clade
support. There is no way to ask it for a bare tree, because a bare tree from a
60-site alignment and one from a 3000-site alignment are indistinguishable to
look at and one of them is wrong.

Read the output in this order:
1. `support.fraction_resolved` — the share of clades clearing {WEAK_SUPPORT:.2f}.
   Below ~0.7 the tree is partially resolved and should be described that way.
2. `warnings` with code `conflicting_support` — clades the DATA support that the
   reported tree does NOT contain. A Newick string cannot express these, so they
   are the finding most often missed.
3. Individual clade support before repeating any specific grouping as a fact.

`select_substitution_model` returns the ranking and the margin, not just a
winner. Models within {DELTA_INDISTINGUISHABLE} of the best are conventionally
indistinguishable; when the tool says so, do not report the winner as "the
best-fitting model" without that caveat.

Alignment LENGTH is not evidence — `n_parsimony_informative` is. This server
infers trees; it does not align sequences, and it will refuse ragged input
rather than guess.
"""

mcp = MCPServer("phylokit-mcp", instructions=INSTRUCTIONS)

_READ_ONLY = ToolAnnotations(
    readOnlyHint=True, idempotentHint=True, openWorldHint=False
)


class CladeDict(TypedDict):
    taxa: list[str]
    clade: str
    support: float
    in_ml_tree: bool


class SupportDict(TypedDict):
    replicates: int
    method: str
    seed: int
    n_clades: int
    n_supported: int
    n_strongly_supported: int
    fraction_resolved: float
    mean_support: float
    min_support: float | None
    support_thresholds: dict[str, float]
    clades: list[CladeDict]
    conflicting_clades: list[CladeDict]


class AlignmentDict(TypedDict):
    n_taxa: int
    n_sites: int
    moltype: str
    n_parsimony_informative: int
    fraction_gaps: float
    duplicate_sequences: list[list[str]]


class TreeResult(TypedDict):
    newick: str
    newick_with_support: str
    model: str
    log_likelihood: float | None
    alignment: AlignmentDict
    support: SupportDict
    reproducibility: dict[str, Any]
    warnings: list[dict[str, Any]]


class ModelResult(TypedDict):
    criterion: str
    best_model: str
    n_models_compared: int
    ranking: list[dict[str, Any]]
    indistinguishable_from_best: list[str]
    best_by_criterion: dict[str, str]
    criteria_agree: bool
    seed: int
    alignment: AlignmentDict
    warnings: list[dict[str, Any]]


class SimulationResult(TypedDict):
    fasta: str
    true_newick: str
    model: str
    alignment: AlignmentDict
    seed: int
    warnings: list[dict[str, Any]]


class CompareResult(TypedDict):
    robinson_foulds: int
    normalised_robinson_foulds: float | None
    identical_topology: bool
    n_shared_clades: int
    only_in_a: list[str]
    only_in_b: list[str]
    shared_taxa: int
    warnings: list[dict[str, Any]]


class CapabilitiesResult(TypedDict):
    engine: str
    engine_version: str
    substitution_models: NotRequired[list[str]]
    n_substitution_models: int
    criteria: list[str]
    limits: dict[str, int]
    support_thresholds: dict[str, float]
    threads_pinned: bool


def _load(fasta: str) -> dict[str, str]:
    seqs = parse_fasta(fasta)
    validate(seqs)
    return seqs


def _annotate(tree, support_by_clade: dict[str, float], tips: set[str]) -> str:
    """Newick with support values on internal nodes.

    Emitted alongside the structured clade list, not instead of it: the string
    is what other tools consume, but it cannot carry the conflicting clades, so
    it is the lossy view by construction.
    """
    for edge in tree.get_edge_vector(include_root=False):
        if edge.is_tip():
            continue
        below = frozenset(t.name for t in edge.tips())
        ref = min(tips)
        side = below if ref not in below else frozenset(tips) - below
        key = format_split(side)
        if key in support_by_clade:
            edge.name = f"{support_by_clade[key]:.2f}"
    return tree.get_newick(with_node_names=True, with_distances=True)


@mcp.tool(
    title="Infer a phylogenetic tree with bootstrap support",
    annotations=_READ_ONLY,
    structured_output=True,
)
def infer_tree(
    fasta: str,
    model: str = "GTR+G",
    replicates: int = DEFAULT_REPLICATES,
    seed: int = 1,
) -> TreeResult:
    """Build a maximum-likelihood tree and measure how well the data support it.

    Always bootstraps. There is deliberately no option to skip it: an
    unsupported topology is the failure mode this server exists to prevent.

    Args:
        fasta: Aligned nucleotide sequences in FASTA. All sequences must be the
            same length — this server does not align.
        model: Substitution model, e.g. "JC", "HKY", "GTR+G". Run `select_model`
            first if you do not have a reason to prefer one.
        replicates: Bootstrap replicates (20-1000). Cost is roughly linear in
            this, so 100 is a reasonable default and 1000 is for a final answer.
        seed: Fixes both the resampling and the engine's search.
    """
    seqs = _load(fasta)
    stats = summarise(seqs)
    tree = build_ml_tree(seqs, model, seed=seed)
    support = bootstrap_support(seqs, tree, model, replicates=replicates, seed=seed)
    by_clade = {c.as_dict()["clade"]: c.support for c in support.clades}
    return TreeResult(
        newick=tree.get_newick(with_distances=True),
        newick_with_support=_annotate(tree, by_clade, set(seqs)),
        model=model,
        log_likelihood=tree_log_likelihood(tree),
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        support=support.as_dict(),  # type: ignore[typeddict-item]
        reproducibility=engine.reproducibility(),
        warnings=diagnostics.collect(
            stats=stats, result=support, pinned=engine.threads_pinned()
        ),
    )


@mcp.tool(
    title="Rank substitution models, with the margin over the runners-up",
    annotations=_READ_ONLY,
    structured_output=True,
)
def select_substitution_model(
    fasta: str, criterion: str = "AIC", seed: int = 1, top_n: int = 5
) -> ModelResult:
    """Compare substitution models and report how much the winner won by.

    A single model name reads as a finding. The ranking, the delta to the next
    model, and whether AIC/AICc/BIC agree are what make it one.

    Args:
        fasta: Aligned nucleotide sequences in FASTA.
        criterion: "AIC", "AICc" or "BIC". BIC penalises parameters more heavily.
        seed: Fixes the engine's search.
        top_n: How many ranked models to return.
    """
    seqs = _load(fasta)
    stats = summarise(seqs)
    selection = select_model(seqs, criterion=criterion, seed=seed, top_n=top_n)
    return ModelResult(
        **selection,  # type: ignore[typeddict-item]
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        warnings=diagnostics.collect(
            stats=stats, selection=selection, pinned=engine.threads_pinned()
        ),
    )


@mcp.tool(
    title="Compare two tree topologies",
    annotations=_READ_ONLY,
    structured_output=True,
)
def compare_trees(newick_a: str, newick_b: str) -> CompareResult:
    """Robinson-Foulds distance between two trees, and which clades differ.

    Compares SPLITS, not strings: the same topology has many valid Newick
    representations, so string equality answers a different question.

    Args:
        newick_a: First tree in Newick format.
        newick_b: Second tree in Newick format.
    """
    from cogent3 import make_tree

    tree_a, tree_b = make_tree(newick_a.strip()), make_tree(newick_b.strip())
    tips_a, tips_b = set(tree_a.get_tip_names()), set(tree_b.get_tip_names())
    shared = tips_a & tips_b
    warnings: list[dict[str, Any]] = []
    if tips_a != tips_b:
        warnings.append(
            {
                "code": "taxon_sets_differ",
                "message": (
                    "The two trees do not have the same taxa, so the distance is "
                    "computed on the shared subset only and is not comparable to "
                    "an RF distance between trees on the full taxon set."
                ),
                "only_in_a": sorted(tips_a - tips_b),
                "only_in_b": sorted(tips_b - tips_a),
            }
        )
    if len(shared) < 4:
        raise ValueError(
            f"Only {len(shared)} shared taxa; at least 4 are needed for an "
            "unrooted tree to have any internal edge to compare."
        )

    sa = canonical_splits(
        tree_a.get_sub_tree(shared) if tips_a != shared else tree_a, shared
    )
    sb = canonical_splits(
        tree_b.get_sub_tree(shared) if tips_b != shared else tree_b, shared
    )
    return CompareResult(
        robinson_foulds=robinson_foulds(sa, sb),
        normalised_robinson_foulds=normalised_robinson_foulds(sa, sb),
        identical_topology=sa == sb,
        n_shared_clades=len(sa & sb),
        only_in_a=sorted(format_split(s) for s in sa - sb),
        only_in_b=sorted(format_split(s) for s in sb - sa),
        shared_taxa=len(shared),
        warnings=warnings,
    )


@mcp.tool(
    title="Simulate an alignment from a known tree",
    annotations=_READ_ONLY,
    structured_output=True,
)
def simulate_alignment(
    newick: str, model: str = "JC", length: int = 500, seed: int = 1
) -> SimulationResult:
    """Generate sequences along a tree you specify, so the true answer is known.

    This is the positive control for everything else here: infer a tree from the
    output and compare it back with `compare_trees`. If inference cannot recover
    a topology you generated from, the problem is the data or the settings, not
    the biology.

    Args:
        newick: The true tree, with branch lengths.
        model: Substitution model to simulate under.
        length: Number of sites.
        seed: Fixes the simulation.
    """
    from cogent3 import make_tree

    if not 1 <= length <= aln_mod.MAX_SITES:
        raise ValueError(f"length must be between 1 and {aln_mod.MAX_SITES}.")
    tree = make_tree(newick.strip())
    aln = engine.piqtree().simulate_alignment(
        tree=tree, model=model, length=length, rand_seed=seed
    )
    seqs = {nm: str(aln.get_seq(nm)) for nm in aln.names}
    stats = summarise(seqs)
    return SimulationResult(
        fasta="".join(f">{nm}\n{s}\n" for nm, s in seqs.items()),
        true_newick=tree.get_newick(with_distances=True),
        model=model,
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        seed=seed,
        warnings=diagnostics.collect(stats=stats, pinned=engine.threads_pinned()),
    )


@mcp.tool(
    title="Engine capabilities and limits",
    annotations=_READ_ONLY,
    structured_output=True,
)
def capabilities(include_models: bool = False) -> CapabilitiesResult:
    """What this server can do, and the bounds it enforces.

    Args:
        include_models: Include the full substitution-model list (long).
    """
    # available_models() returns a cogent3 Table, not a mapping. Reading it
    # under a bare `except` previously turned that into a silent 0, which is a
    # worse answer than an error: it reads as "this engine has no models".
    models = sorted(
        str(v) for v in engine.piqtree().available_models().columns["Abbreviation"]
    )
    out = CapabilitiesResult(
        engine="IQ-TREE 2 via piqtree",
        engine_version=engine.engine_version(),
        n_substitution_models=len(models),
        criteria=list(CRITERIA),
        limits={
            "min_taxa": aln_mod.MIN_TAXA,
            "max_taxa": aln_mod.MAX_TAXA,
            "max_sites": aln_mod.MAX_SITES,
            "min_replicates": MIN_REPLICATES,
            "max_replicates": MAX_REPLICATES,
        },
        support_thresholds={"supported": WEAK_SUPPORT, "strong": STRONG_SUPPORT},
        threads_pinned=engine.threads_pinned(),
    )
    if include_models:
        out["substitution_models"] = models
    return out


def main() -> None:
    mcp.run()


# `python -m phylokit_mcp.server` must work, not just the console script — it is
# what a client config most often points at, and it is what the protocol test
# launches.
if __name__ == "__main__":
    main()


__all__ = ["AlignmentError", "main", "mcp"]
