"""MCP server. Five tools over IQ-TREE, one over MAFFT, and a position mapper.

`engine` is imported FIRST and deliberately: it pins OMP_NUM_THREADS before
piqtree loads, and IQ-TREE reads that variable when its thread pool initialises.
Import order here is load-bearing, not incidental.

Tools are synchronous, and they run one at a time under `_ENGINE_LOCK`.
MCPServer runs a sync tool in a worker thread (`anyio.to_thread`), so calls a
client sends together run concurrently unless something stops them, and
piqtree's engine wrapper changes process-global state for the length of a call:
it chdirs into a temporary directory it then deletes, and points file
descriptors 1 and 2 at /dev/null. Interleaved, two calls left the process in a
deleted directory and every later call failed (the 2026-10-18 benchmark re-run,
three `infer_tree` calls in one turn). The fd swap does not reach the JSON-RPC
stream: mcp's stdio transport serves it from a private duplicate of fd 1 and
points fd 1 itself at stderr, so a progress notification sent during an engine
call still arrives (`tests/test_mcp_protocol.py` checks it on a real process).

Every tool sends `notifications/progress` every `HEARTBEAT_S` while it runs or
queues for the lock, when the client asked for progress. A client with an idle
timeout otherwise aborts a long call that is working: on the field benchmark
every 89-taxon `infer_tree` went silent for the whole bootstrap and Claude Code
aborted it at 1822 s, ten runs in ten.
"""

# isort: off  -- engine MUST come before anything that could pull in piqtree.
from . import engine  # noqa: I001

# isort: on
# stdlib TypedDict is safe here: pydantic only rejects it under Python < 3.12,
# and piqtree ships wheels for 3.12+ only, so this package cannot run on a
# version where the typing_extensions backport would be required.
import functools
import inspect
import os
import threading
import time
import typing
from typing import Any, Literal, NotRequired, TypedDict

import anyio.from_thread

# mcp 2.x renamed FastMCP to MCPServer and moved it out of mcp.server.fastmcp,
# which no longer exists. The class is the same one — same decorator, same
# `annotations` and `structured_output` kwargs — so this is a rename, not a
# rewrite. ToolAnnotations below did not move.
from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from . import __version__, diagnostics, msa, positions
from . import alignment as aln_mod
from .inputs import text_or_path, write_new_file
from .alignment import (
    MOLTYPES,
    AlignmentError,
    parse_fasta,
    require_phylogenetic_signal,
    summarise,
    validate,
)
from .bootstrap import (
    MAX_REPLICATES,
    MIN_REPLICATES,
    STRONG_SUPPORT,
    UFBOOT_MAX_REPLICATES,
    UFBOOT_MIN_REPLICATES,
    UFBOOT_SUPPORT,
    WEAK_SUPPORT,
    bootstrap_support,
    default_replicates,
    ufboot_tree,
    validate_replicates,
)
from .inference import (
    CRITERIA,
    DELTA_INDISTINGUISHABLE,
    build_ml_tree,
    model_moltype,
    require_model_for,
    select_model,
    tree_log_likelihood,
)
from .newick import parse_newick
from .splits import (
    canonical_splits,
    format_split,
    normalised_robinson_foulds,
    robinson_foulds,
    splits_to_newick,
)

INSTRUCTIONS = f"""\
Phylogenetic inference over IQ-TREE 3, through piqtree.

The rule this server is built around: **a topology without support is not a
result.** `infer_tree` always runs a bootstrap and always returns per-clade
support. There is no way to ask it for a bare tree, because a bare tree from a
60-site alignment and one from a 3000-site alignment are indistinguishable to
look at and one of them is wrong.

Two bootstraps. `support_method="bootstrap"` (the default) is the nonparametric
one: each replicate is a full ML search, so its cost is replicates x one search
(on 89 taxa x 435 sites, 25-93 s per replicate: 100 replicates is over an hour).
`support_method="ufboot"` is IQ-TREE's UFBoot2: 1000 replicates in about the
time of two or three searches (83 s on the same data). Their values are on
different scales: nonparametric clades are supported at {WEAK_SUPPORT:.2f},
UFBoot2 clades at {UFBOOT_SUPPORT:.2f} (Hoang et al. 2018). Each result names its
method and its own thresholds in `support`. Only the nonparametric method can
report `conflicting_support`.

Read the output in this order:
1. `support.fraction_resolved` — the share of clades clearing the method's line.
   Below ~0.7 the tree is partially resolved and should be described that way.
2. `warnings` with code `conflicting_support` — clades the DATA support that the
   reported tree does NOT contain. A Newick string cannot express these, so they
   are the finding most often missed.
3. Individual clade support before repeating any specific grouping as a fact.

To report only the supported clades, use `newick_supported`: the ML tree with
every clade below the method's line collapsed into a polytomy. For the full ML
topology without branch lengths, use `newick_topology`. Copy either as is; do
not rebuild a tree by hand.

`select_substitution_model` returns the ranking and the margin, not just a
winner. Models within {DELTA_INDISTINGUISHABLE} of the best are conventionally
indistinguishable; when the tool says so, do not report the winner as "the
best-fitting model" without that caveat.

Alignment LENGTH is not evidence — `n_parsimony_informative` is. The tree tools
refuse ragged input rather than guess. `align_sequences` is the step before
them: it runs MAFFT on unaligned sequences and returns FASTA that `infer_tree`
accepts as is. An alignment is an inference too — read its `fraction_gaps` and
`n_parsimony_informative` before building on it.

Data in a file: pass its path (`fasta_path`, `trees_path`) instead of the text.
The server reads the file itself, so there is no need to read it first or to
paste it into the call. `align_sequences(out_path=...)` writes the alignment
to a new file instead of returning it; pass that path to `infer_tree` or
`map_positions` as `fasta_path`, rather than copying rows between calls.
`map_positions` answers which residue of one sequence
shares an alignment column with a residue of another; do not count columns by
hand.
"""

# snake_case spellings. The camelCase aliases are still accepted as constructor
# kwargs and build an EQUAL object -- measured, not assumed -- so this is a
# consistency change with the sibling servers, not a bug fix.
_READ_ONLY = ToolAnnotations(
    read_only_hint=True, idempotent_hint=True, open_world_hint=False
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
    # Topology only: the ML tree with every clade below `support`'s own
    # "supported" line collapsed, which is what a caller reports.
    newick_supported: str
    # The ML tree with no branch lengths, for an answer that wants topology.
    newick_topology: str
    model: str
    log_likelihood: float | None
    alignment: AlignmentDict
    support: SupportDict
    reproducibility: dict[str, Any]
    # Which engine produced this, travelling WITH the tree. The version was
    # already reachable via `capabilities`, but a stored or forwarded newick
    # could not say what built it, and two trees could not be told apart.
    engine: dict[str, Any]
    # Newick carries bare numbers after the colons. They are substitutions per
    # site — not time, not percent divergence — and the distinction decides
    # whether a distance read off this tree means anything.
    branch_length_units: str
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
    # Only when `models` was given: the IQ-TREE names ranked, and any textbook
    # name (JC69, K80, HKY85) that was read as one of them.
    restricted_to: NotRequired[list[str]]
    resolved_names: NotRequired[dict[str, str]]


class SimulationResult(TypedDict):
    fasta: str
    true_newick: str
    model: str
    alignment: AlignmentDict
    seed: int
    warnings: list[dict[str, Any]]


class AlignResult(TypedDict):
    # None when `out_path` was given: the rows are in the file at `written_to`.
    fasta: str | None
    written_to: str | None
    alignment: AlignmentDict
    sequence_type: str
    # Longest and shortest INPUT sequence. Against `alignment.n_sites` this is
    # how much of the alignment is gap the aligner introduced.
    input_lengths: dict[str, int]
    engine: dict[str, Any]
    ready_for_infer_tree: bool
    warnings: list[dict[str, Any]]


class PositionsResult(TypedDict):
    column: int
    n_columns: int
    query: dict[str, Any]
    # Per sequence: its residue position and residue at `column`, or at a gap
    # position None, residue "-" and `between`, the residues either side.
    rows: dict[str, dict[str, Any]]
    window: dict[str, Any]
    # None when the input was already aligned; MAFFT's version and options
    # when this call aligned it.
    aligned_by: dict[str, Any] | None
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
    # piqtree's version, not IQ-TREE's; the key name is kept for existing
    # callers. IQ-TREE's own version is `iqtree_version`.
    engine_version: str
    iqtree_version: str
    substitution_models: NotRequired[list[str]]
    n_substitution_models: int
    criteria: list[str]
    # MAFFT's version string, or None when it is not installed. None is an
    # answer, not a failure: only `align_sequences` needs it.
    aligner_version: str | None
    # Set only when MAFFT is installed and broken, so that state does not read
    # as "not installed". The tree tools are unaffected either way.
    aligner_error: str | None
    limits: dict[str, int]
    support_thresholds: dict[str, float]
    threads_pinned: bool


def _load(fasta: str, sequence_type: str = "dna") -> dict[str, str]:
    if sequence_type not in MOLTYPES:
        raise AlignmentError(
            f"Unknown sequence_type {sequence_type!r}. Valid: {list(MOLTYPES)}."
        )
    seqs = parse_fasta(fasta)
    validate(seqs, sequence_type)
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


def infer_tree(
    fasta: str | None = None,
    model: str = "GTR+G",
    replicates: int | None = None,
    seed: int = 1,
    sequence_type: str = "dna",
    fasta_path: str | None = None,
    support_method: Literal["bootstrap", "ufboot"] = "bootstrap",
) -> TreeResult:
    """Build a maximum-likelihood tree and measure how well the data support it.

    Always bootstraps. There is deliberately no option to skip it: an
    unsupported topology is the failure mode this server exists to prevent.
    `support_method` picks which bootstrap; see its entry for the cost.

    Args:
        fasta: Aligned sequences in FASTA, nucleotide or protein as declared by
            `sequence_type`. All sequences must be the same length — run
            `align_sequences` first if they are not.
        model: Substitution model, e.g. "JC", "HKY", "GTR+G" (or "LG+G" for
            protein). Run `select_substitution_model` first if you do not have a
            reason to prefer one.
        replicates: Bootstrap replicates. Nonparametric: 20-1000, default 100,
            and each one costs a full ML search, so the call takes about
            replicates x one search. UFBoot2: 1000-10000, default 1000.
        seed: 0 to 2**31-1. Fixes the column resampling exactly; the engine's
            search is seeded too but is not bit-exact (see `reproducibility`).
        sequence_type: "dna" (default) or "protein". DECLARED, never sniffed: an
            alignment of only A/C/G/T is a valid protein alignment too, so
            guessing would silently fit a nucleotide model to protein data. A
            protein alignment also needs a protein model — "LG+G" or "WAG",
            not the nucleotide default — so run select_substitution_model with
            the same sequence_type first.
        fasta_path: Path of a FASTA file to read instead of `fasta`.
        support_method: "bootstrap" (default): the nonparametric bootstrap,
            clades supported at 0.70, and the only method that reports
            well-supported clades the ML tree lacks (`conflicting_support`).
            "ufboot": IQ-TREE's ultrafast bootstrap (UFBoot2, Hoang et al.
            2018), clades supported at 0.95; far faster, use it when the
            nonparametric cost (replicates x one ML search) is too high, e.g.
            more than a few dozen taxa.
    """
    seqs = _load(text_or_path(fasta, fasta_path, "fasta", "fasta_path"), sequence_type)
    if replicates is None:
        replicates = default_replicates(support_method)
    # Every argument check runs BEFORE the maximum-likelihood search. The
    # bootstrap used to be the first thing to look at `replicates` and `seed`,
    # so a bad value was refused only after the most expensive single step had
    # already run.
    engine.validate_seed(seed)
    validate_replicates(replicates, support_method)
    require_model_for(model, sequence_type)
    # Only tree inference needs this. `_load` is shared with tools that can
    # legitimately work on a signal-free alignment, so the check lives here
    # rather than inside validate().
    require_phylogenetic_signal(seqs, sequence_type)
    stats = summarise(seqs, sequence_type)
    if support_method == "ufboot":
        tree, support = ufboot_tree(
            seqs, model, replicates=replicates, seed=seed, moltype=sequence_type
        )
    else:
        tree = build_ml_tree(seqs, model, seed=seed, moltype=sequence_type)
        support = bootstrap_support(
            seqs, tree, model, replicates=replicates, seed=seed, moltype=sequence_type
        )
    by_clade = {c.as_dict()["clade"]: c.support for c in support.clades}
    newick = tree.get_newick(with_distances=True)
    supported = {
        frozenset(c.taxa) for c in support.clades if c.support >= support.supported_at
    }
    return TreeResult(
        newick=newick,
        newick_with_support=_annotate(tree, by_clade, set(seqs)),
        newick_supported=splits_to_newick(set(seqs), supported),
        newick_topology=tree.get_newick(with_distances=False),
        model=model,
        log_likelihood=tree_log_likelihood(tree),
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        support=support.as_dict(),  # type: ignore[typeddict-item]
        reproducibility=engine.reproducibility(),
        engine={
            "name": "IQ-TREE 3 via piqtree",
            "version": engine.engine_version(),
            "iqtree_version": engine.iqtree_version(),
        },
        branch_length_units="substitutions per site",
        warnings=diagnostics.collect(
            stats=stats,
            result=support,
            pinned=engine.threads_pinned(),
            newick=newick,
        ),
    )


def select_substitution_model(
    fasta: str | None = None,
    criterion: str = "AIC",
    seed: int = 1,
    top_n: int = 5,
    sequence_type: str = "dna",
    fasta_path: str | None = None,
    models: list[str] | None = None,
) -> ModelResult:
    """Compare substitution models and report how much the winner won by.

    A single model name reads as a finding. The ranking, the delta to the next
    model, and whether AIC/AICc/BIC agree are what make it one.

    Args:
        fasta: Aligned sequences in FASTA, nucleotide or protein as declared by
            `sequence_type`.
        criterion: "AIC", "AICc" or "BIC". BIC penalises parameters more heavily.
        seed: Fixes the engine's search.
        top_n: How many ranked models to return.
        sequence_type: "dna" (default) or "protein". Ranks within that molecule
            type's model set — nucleotide and protein models are not comparable.
        fasta_path: Path of a FASTA file to read instead of `fasta`.
        models: Rank only these models (IQ-TREE names, e.g. ["JC", "K2P",
            "F81+F", "HKY+F"]), all returned. Bare "F81"/"HKY" have equal base
            frequencies in IQ-TREE; JC69, K80 and HKY85 are accepted.
    """
    seqs = _load(text_or_path(fasta, fasta_path, "fasta", "fasta_path"), sequence_type)
    engine.validate_seed(seed)
    if criterion not in CRITERIA:
        raise ValueError(f"criterion must be one of {CRITERIA}, got {criterion!r}")
    stats = summarise(seqs, sequence_type)
    selection = select_model(
        seqs,
        criterion=criterion,
        seed=seed,
        top_n=top_n,
        moltype=sequence_type,
        models=models,
    )
    extra = selection.pop("model_warnings", [])
    return ModelResult(
        **selection,  # type: ignore[typeddict-item]
        alignment=stats.__dict__,
        warnings=diagnostics.collect(
            stats=stats, selection=selection, pinned=engine.threads_pinned()
        )
        + extra,
    )


def _two_trees(text: str) -> tuple[str, str]:
    """The two Newick trees in a `trees_path` file, A then B."""
    trees = [t.strip() + ";" for t in text.split(";") if t.strip()]
    if len(trees) != 2:
        raise ValueError(
            f"trees_path must hold exactly two Newick trees (each ending in ';'), "
            f"tree A then tree B; found {len(trees)}."
        )
    return trees[0], trees[1]


def compare_trees(
    newick_a: str | None = None,
    newick_b: str | None = None,
    trees_path: str | None = None,
) -> CompareResult:
    """Robinson-Foulds distance between two trees, and which clades differ.

    Compares SPLITS, not strings: the same topology has many valid Newick
    representations, so string equality answers a different question.

    Args:
        newick_a: First tree in Newick format.
        newick_b: Second tree in Newick format.
        trees_path: Instead of the two strings, the path of a file holding
            exactly two Newick trees, A then B.
    """
    if trees_path is not None:
        if newick_a is not None or newick_b is not None:
            raise ValueError("Give newick_a and newick_b, or trees_path, not both.")
        newick_a, newick_b = _two_trees(
            text_or_path(None, trees_path, "", "trees_path")
        )
    if newick_a is None or newick_b is None:
        raise ValueError("Give both newick_a and newick_b, or trees_path.")
    tree_a = parse_newick(newick_a, "newick_a")
    tree_b = parse_newick(newick_b, "newick_b")
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


def simulate_alignment(
    newick: str, model: str = "JC", length: int = 500, seed: int = 1
) -> SimulationResult:
    """Generate sequences along a tree you specify, so the true answer is known.

    This is the positive control for everything else here: infer a tree from the
    output and compare it back with `compare_trees`. If inference cannot recover
    a topology you generated from, the problem is the data or the settings, not
    the biology.

    Args:
        newick: The true tree, with a branch length on every edge.
        model: Substitution model to simulate under. A protein model ("LG",
            "WAG", ...) simulates protein; pass the output to `infer_tree` with
            sequence_type="protein". `alignment.moltype` says which it is.
        length: Number of sites.
        seed: Fixes the simulation.
    """
    if not 1 <= length <= aln_mod.MAX_SITES:
        raise ValueError(f"length must be between 1 and {aln_mod.MAX_SITES}.")
    engine.validate_seed(seed)
    # The molecule type comes from the model, as piqtree classifies it. It was
    # hard-wired to DNA, so an LG simulation came back summarised in the DNA
    # alphabet: 3 parsimony-informative sites where the protein count was 48.
    moltype = model_moltype(model)
    tree = parse_newick(newick)
    tips = tree.get_tip_names()
    if not aln_mod.MIN_TAXA <= len(tips) <= aln_mod.MAX_TAXA:
        raise ValueError(
            f"The tree has {len(tips)} tips; simulate_alignment takes "
            f"{aln_mod.MIN_TAXA}-{aln_mod.MAX_TAXA}, the same range infer_tree "
            "accepts, so that its output can be fed straight back in."
        )
    missing = [
        e.name for e in tree.get_edge_vector(include_root=False) if e.length is None
    ]
    if missing:
        raise ValueError(
            f"{len(missing)} edge(s) have no branch length (e.g. "
            f"{sorted(missing)[:4]}). Branch lengths are the amount of change to "
            "simulate; without them every sequence comes back identical to the "
            "root, which is an alignment with no signal, not a simulation."
        )
    aln = engine.piqtree().simulate_alignment(
        tree=tree, model=model, length=length, rand_seed=seed
    )
    seqs = {nm: str(aln.get_seq(nm)) for nm in aln.names}
    stats = summarise(seqs, moltype)
    return SimulationResult(
        fasta="".join(f">{nm}\n{s}\n" for nm, s in seqs.items()),
        true_newick=tree.get_newick(with_distances=True),
        model=model,
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        seed=seed,
        warnings=diagnostics.collect(stats=stats, pinned=engine.threads_pinned()),
    )


def align_sequences(
    fasta: str | None = None,
    sequence_type: str = "dna",
    fasta_path: str | None = None,
    out_path: str | None = None,
) -> AlignResult:
    """Align unaligned sequences with MAFFT, ready for `infer_tree`.

    The returned `fasta` has every row the same length and can be passed to
    `infer_tree` or `select_substitution_model` unchanged. Each output row, with
    its gaps removed, is verified to equal the input sequence before it is
    returned.

    Args:
        fasta: UNALIGNED sequences in FASTA, 2-200 of them. Gap characters are
            refused: input that is already aligned does not need this tool.
        sequence_type: "dna" (default) or "protein". Declared, never sniffed,
            and passed to MAFFT explicitly so it does not guess either.
        fasta_path: Path of a FASTA file to read instead of `fasta`.
        out_path: Write the aligned FASTA to this NEW file (never overwritten;
            its directory must exist) instead of returning it: `fasta` comes
            back None and `written_to` holds the absolute path, which
            `infer_tree`, `select_substitution_model` and `map_positions` take
            as `fasta_path`. Use it rather than copying rows between calls.
    """
    if sequence_type not in MOLTYPES:
        raise AlignmentError(
            f"Unknown sequence_type {sequence_type!r}. Valid: {list(MOLTYPES)}."
        )
    seqs = parse_fasta(text_or_path(fasta, fasta_path, "fasta", "fasta_path"))
    aligned = msa.align(seqs, sequence_type)
    stats = summarise(aligned, sequence_type)
    warnings = diagnostics.collect(stats=stats, pinned=engine.threads_pinned())
    # `ready_for_infer_tree` is infer_tree's own input checks run on the output,
    # not a second opinion about them. It used to be a taxon count alone, and
    # said True for output infer_tree then could not parse.
    ready = len(aligned) >= aln_mod.MIN_TAXA
    if ready:
        try:
            validate(aligned, sequence_type)
            require_phylogenetic_signal(aligned, sequence_type)
        except AlignmentError as exc:
            ready = False
            warnings.append({"code": "infer_tree_would_refuse", "message": str(exc)})
    else:
        warnings.append(
            {
                "code": "too_few_taxa_for_a_tree",
                "message": (
                    f"{len(aligned)} sequences were aligned, but infer_tree needs "
                    f"at least {aln_mod.MIN_TAXA}: an unrooted tree on fewer has "
                    "only one topology."
                ),
            }
        )
    lengths = [len(s) for s in seqs.values()]
    text = "".join(f">{nm}\n{s}\n" for nm, s in aligned.items())
    written = None if out_path is None else write_new_file(out_path, text, "out_path")
    return AlignResult(
        fasta=None if written else text,
        written_to=written,
        alignment=stats.__dict__,  # type: ignore[typeddict-item]
        sequence_type=sequence_type,
        input_lengths={"min": min(lengths), "max": max(lengths)},
        engine={
            "name": "MAFFT",
            "version": msa.mafft_version(),
            "options": msa.mafft_args(sequence_type),
        },
        ready_for_infer_tree=ready,
        warnings=warnings,
    )


def map_positions(
    aligned: bool,
    sequence: str | None = None,
    position: int | None = None,
    column: int | None = None,
    fasta: str | None = None,
    fasta_path: str | None = None,
    sequence_type: str = "dna",
) -> PositionsResult:
    """Map a residue to its alignment column, and the column to every sequence.

    Answers "which residue of B is in the same column as residue N of A" without
    counting by hand. Positions count residues from 1, gaps not counted; columns
    count from 1. At a gap, a row's position is null and `between` gives the
    residues either side.

    Args:
        aligned: true if the input is an alignment; false if unaligned, and it
            is aligned with MAFFT first, as `align_sequences` would.
        sequence: The sequence whose residue `position` is meant.
        position: A residue position in `sequence`. Give this or `column`.
        column: An alignment column. Give this or `position`.
        fasta: The sequences in FASTA.
        fasta_path: Path of a FASTA file to read instead of `fasta`.
        sequence_type: "dna" (default) or "protein".
    """
    if sequence_type not in MOLTYPES:
        raise AlignmentError(
            f"Unknown sequence_type {sequence_type!r}. Valid: {list(MOLTYPES)}."
        )
    seqs = parse_fasta(text_or_path(fasta, fasta_path, "fasta", "fasta_path"))
    warnings: list[dict[str, Any]] = []
    aligned_by = None
    if aligned:
        bad = aln_mod.unrecognised_characters(
            {c for s in seqs.values() for c in s} - {positions.GAP, "?", "."},
            sequence_type,
        )
        if bad:
            raise AlignmentError(
                f"Unrecognised characters for {sequence_type}: {sorted(bad)[:8]}."
            )
    else:
        seqs = msa.align(seqs, sequence_type)
        aligned_by = {
            "name": "MAFFT",
            "version": msa.mafft_version(),
            "options": msa.mafft_args(sequence_type),
        }
        warnings.append(
            {
                "code": "column_is_an_inference",
                "message": (
                    "The column comes from this MAFFT alignment. Another aligner or "
                    "settings can place a residue elsewhere, most often inside "
                    "gappy stretches: check `window` before relying on it."
                ),
            }
        )
    out = positions.map_column(seqs, sequence, position, column)
    return PositionsResult(**out, aligned_by=aligned_by, warnings=warnings)  # type: ignore[typeddict-item]


def capabilities(include_models: bool = False) -> CapabilitiesResult:
    """What this server can do, and the bounds it enforces.

    `engine_version` is the installed piqtree version (e.g. "0.8.3");
    `iqtree_version` is the IQ-TREE build inside it (e.g. "3.1.2").

    Args:
        include_models: Include the full substitution-model list (long).
    """
    # available_models() returns a cogent3 Table, not a mapping. Reading it
    # under a bare `except` previously turned that into a silent 0, which is a
    # worse answer than an error: it reads as "this engine has no models".
    models = sorted(
        str(v) for v in engine.piqtree().available_models().columns["Abbreviation"]
    )
    aligner_version, aligner_error = None, None
    try:
        aligner_version = msa.mafft_version()
    except msa.MafftFailedError as exc:
        aligner_error = str(exc)
    out = CapabilitiesResult(
        engine="IQ-TREE 3 via piqtree",
        engine_version=engine.engine_version(),
        iqtree_version=engine.iqtree_version(),
        n_substitution_models=len(models),
        criteria=list(CRITERIA),
        aligner_version=aligner_version,
        aligner_error=aligner_error,
        limits={
            "min_sequences_to_align": msa.MIN_SEQUENCES,
            "min_taxa": aln_mod.MIN_TAXA,
            "max_taxa": aln_mod.MAX_TAXA,
            "max_sites": aln_mod.MAX_SITES,
            "min_replicates": MIN_REPLICATES,
            "max_replicates": MAX_REPLICATES,
            "ufboot_min_replicates": UFBOOT_MIN_REPLICATES,
            "ufboot_max_replicates": UFBOOT_MAX_REPLICATES,
        },
        support_thresholds={
            "supported": WEAK_SUPPORT,
            "strong": STRONG_SUPPORT,
            "ufboot_supported": UFBOOT_SUPPORT,
        },
        threads_pinned=engine.threads_pinned(),
    )
    if include_models:
        out["substitution_models"] = models
    return out


# One tool call at a time per process; see the module docstring.
_ENGINE_LOCK = threading.Lock()

# Seconds between progress notifications. Well under any client idle timeout
# worth having; the env override exists for the protocol test.
HEARTBEAT_S = float(os.environ.get("PHYLOKIT_HEARTBEAT_S", "15"))


def _with_heartbeat(call, ctx: Context | None, name: str):
    """Run `call`, sending progress every HEARTBEAT_S until it returns.

    The work runs on a helper thread so this worker thread stays free to send.
    It covers the wait for `_ENGINE_LOCK` too, because `call` takes the lock.
    `report_progress` is a no-op when the client sent no progress token, so a
    client that did not ask gets nothing extra. The value is elapsed seconds,
    which only increases, as the protocol requires.

    A Context with no request behind it (`mcp.call_tool` in-process) has no
    client to notify, so it gets no heartbeat: `report_progress` raises there.
    A heartbeat that fails while sending fails the call, but only after the
    work has finished and released the lock. Raising at once left the engine
    call running on its helper thread, and the next call raced it.
    """
    if ctx is None or HEARTBEAT_S <= 0 or not _has_request(ctx):
        return call()
    box: dict[str, Any] = {}

    def target() -> None:
        try:
            box["value"] = call()
        except BaseException as exc:  # noqa: BLE001 - re-raised on the caller's thread
            box["error"] = exc

    start = time.monotonic()
    worker = threading.Thread(target=target, name=f"phylokit-{name}", daemon=True)
    worker.start()
    while True:
        worker.join(HEARTBEAT_S)
        if not worker.is_alive():
            break
        elapsed = round(time.monotonic() - start, 1)
        try:
            anyio.from_thread.run(
                ctx.report_progress,
                elapsed,
                None,
                f"{name} still running, {elapsed:.0f} s",
            )
        except Exception as exc:
            worker.join()
            raise RuntimeError(
                f"{name}: progress notification failed after {elapsed:.0f} s "
                f"({type(exc).__name__}: {exc}); the call ran to completion "
                "but its result was not returned."
            ) from exc
    if "error" in box:
        raise box["error"]
    return box["value"]


def _has_request(ctx: Context) -> bool:
    """Whether `ctx` belongs to a client request (the SDK raises if not)."""
    try:
        _ = ctx.request_context
    except ValueError:
        return False
    return True


def _surfaces_refusals(fn):
    """Re-raise an anticipated refusal as ToolError so its text reaches the model.

    Since mcp 2.1 (python-sdk #3314) MCPServer treats any exception other than
    ToolError as a crash: the model sees only `Error executing tool <name>` and
    the reason stays in the server log. Under mcp 2.0 the text went through
    regardless, which is why nothing here needed to say so. Every refusal this
    server raises on purpose is a ValueError — AlignmentError and the argument
    checks in bootstrap and inference — and the message is the product: a ragged
    alignment, an unknown criterion, too few replicates, each says what to do
    instead. The conversion happens once, at registration, so the tool functions
    keep raising their own types for the unit tests that import them directly.
    Anything else IS a crash and stays masked as the SDK intends.

    The libraries underneath do not share that convention: cogent3 reports bad
    input as AlphabetError or TreeParseError, piqtree an out-of-range seed as a
    pybind TypeError and a model for the wrong molecule as a RuntimeError. None
    is a ValueError, so each reached the caller as a bare crash. The fix is not
    to widen this net -- a TypeError or RuntimeError is also what a real bug
    raises -- but to classify at the input boundary, where it is known that the
    caller's argument is at fault: `alignment.to_cogent3` and `newick.parse_newick`
    translate the parse errors into ValueError subclasses, and the seed and
    model/molecule checks run before the engine is called at all.
    """

    def locked(*args, **kwargs):
        with _ENGINE_LOCK:
            return fn(*args, **kwargs)

    @functools.wraps(fn)
    def wrapper(*args, ctx: Context | None = None, **kwargs):
        try:
            return _with_heartbeat(
                functools.partial(locked, *args, **kwargs), ctx, fn.__name__
            )
        except ValueError as exc:
            raise ToolError(str(exc)) from exc

    # MCPServer injects a Context into the parameter annotated with it and
    # leaves it out of the tool's input schema. `functools.wraps` copies the
    # tool's own signature and annotations, which have no such parameter, so
    # both are extended here rather than adding `ctx` to every tool function.
    sig = inspect.signature(fn)
    ctx_param = inspect.Parameter(
        "ctx", inspect.Parameter.KEYWORD_ONLY, default=None, annotation=Context | None
    )
    wrapper.__signature__ = sig.replace(
        parameters=[*sig.parameters.values(), ctx_param]
    )
    wrapper.__annotations__ = {**typing.get_type_hints(fn), "ctx": Context | None}
    return wrapper


def build_server() -> MCPServer:
    """Construct a fresh server.

    A factory rather than a module-level singleton, matching breedsim-mcp and
    plantcv-mcp: a test that needs a server gets its own, and importing this
    module no longer has the side effect of constructing one.

    The tool functions stay at module level and are registered here. The sibling
    servers nest their tool definitions inside this function, which they can do
    because their tools are thin wrappers over logic that lives in other modules.
    Here the tool functions ARE the implementation and are unit-tested by direct
    import, so nesting them would make the logic unreachable. Registering instead
    of nesting keeps both properties, and puts the whole tool surface in one place.
    """
    # `version` defaults to "" in MCPServer, which is what serverInfo carried.
    mcp = MCPServer("phylokit-mcp", version=__version__, instructions=INSTRUCTIONS)
    mcp.tool(
        title="Infer a phylogenetic tree with bootstrap support",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(infer_tree))
    mcp.tool(
        title="Rank substitution models, with the margin over the runners-up",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(select_substitution_model))
    mcp.tool(
        title="Compare two tree topologies",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(compare_trees))
    mcp.tool(
        title="Simulate an alignment from a known tree",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(simulate_alignment))
    mcp.tool(
        title="Align sequences with MAFFT, ready for infer_tree",
        # Not read-only since `out_path`: it can create a file. Never destructive:
        # an existing file is refused, not overwritten.
        annotations=ToolAnnotations(
            read_only_hint=False,
            destructive_hint=False,
            idempotent_hint=False,
            open_world_hint=False,
        ),
        structured_output=True,
    )(_surfaces_refusals(align_sequences))
    mcp.tool(
        title="Map a residue to its alignment column and across sequences",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(map_positions))
    mcp.tool(
        title="Engine capabilities and limits",
        annotations=_READ_ONLY,
        structured_output=True,
    )(_surfaces_refusals(capabilities))
    return mcp


def main() -> None:
    build_server().run()


# `python -m phylokit_mcp.server` must work, not just the console script — it is
# what a client config most often points at, and it is what the protocol test
# launches.
if __name__ == "__main__":
    main()


__all__ = ["AlignmentError", "build_server", "main"]
