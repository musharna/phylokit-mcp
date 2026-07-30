# Changelog

All notable changes to `phylokit-mcp` are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Phase 1. Phylogenetic inference over MCP, driving IQ-TREE 2 through piqtree 0.8.3.

### Added

- **The structural rule: `infer_tree` will not return a bare topology.** It
  always bootstraps and always reports per-clade support. There is no flag to
  skip it, because the failure this server exists to prevent is a confident tree
  from data that does not support one.

  Measured on alignments simulated from a known 7-taxon tree: at 300 sites
  inference recovers the truth exactly with every clade at 1.00, while at 60
  sites it returns a tree containing a clade that does not exist (`C,D,G`) and
  missing one that does (`E,F,G`) — RF 2 from the truth. Both are fully resolved
  Newick strings of the same shape. The only thing distinguishing them is that
  the false clade carries **0.57** support against 0.80–1.00 for the true ones.

- **Bootstrap support computed in-server** (Felsenstein 1985), because piqtree
  0.8.3 runs `bootstrap_replicates` but does not expose the values: internal
  nodes come back carrying only `mprobs`, node names are `edge.N`, and nothing
  is written to disk.

  Computing it here has a benefit worth keeping regardless — the full split
  frequency distribution is available, including **conflicting clades**: groups
  the data support at ≥0.70 that are absent from the reported tree. A
  support-annotated Newick string cannot express those, because such a clade has
  nowhere to attach, so the standard output format drops them silently.

- **Model selection with the runners-up.** `select_substitution_model` ranks
  every candidate with ΔAIC/AICc/BIC rather than naming a winner. On a 400-site
  alignment simulated under **JC**, the AIC winner is **F81** — a model the data
  were not generated under — with several models inside the conventional ±2
  margin. A winner without its margin is a claim the numbers do not support.
  The server also reports whether the three criteria agree.

- **Length and evidence reported separately.** `n_parsimony_informative`
  alongside `n_sites`, since a long alignment of near-identical sequences
  supports nothing. Sites where only one taxon differs are excluded — they carry
  no topological signal, and counting them would overstate the evidence.

- **`compare_trees` compares splits, not strings.** The same topology has many
  valid Newick representations, so string equality answers a different question
  from the one a caller means.

- **`simulate_alignment`** — generate sequences along a specified tree, so the
  true answer is known. This is the positive control for everything else here,
  and it is what the test suite is built on.

- Five MCP tools with typed returns and read-only annotations, six advisories,
  62 tests against real IQ-TREE with no mocked engine, and
  `docs/MUTATION-CHECKS.md` recording ten mutants.

- **A real-process protocol test.** The server is launched as a subprocess and
  driven with JSON-RPC over stdio, because every other test calls the tool
  functions in-process and would miss a broken entry point, a crash during tool
  registration, or anything written to stdout at import — which would corrupt
  the JSON-RPC stream. Confirmed to fail when the entry point is removed.

### Notes

- **Reproducibility is reported, not claimed.** Across fresh processes the same
  request reproduces byte-identically. Within one long-lived process it does
  not: passing the same `rand_seed` does not fully reset IQ-TREE's internal
  state — the same tree built three times gave call 1 == call 2 but call 3
  different. Measured over six repeated 50-replicate calls, three of four clades
  were bit-identical and one moved 0.02, one replicate flipping, well inside the
  bootstrap's own sampling error. An earlier version returned `reproducible:
true`, which a caller would reasonably read as a promise it could not keep.

- **The `thin_information` advisory explains rather than predicts.** It fired
  once on a tree whose every clade sat at 1.00 support, because 7.3 informative
  sites per taxon was below a threshold picked a priori. Where a bootstrap has
  run, support has been _measured_, and a proxy that disagrees with the
  measurement is simply wrong; the advisory is now suppressed whenever no clade
  is unsupported. No choice of threshold fixes that class of error, since the
  proxy and the measurement are different quantities.

- **Licensed GPL-2.0-only**, and the "only" is load-bearing: piqtree declares
  `GPL-2.0-only`, which is incompatible with GPL-3.0, so the distributed
  combination cannot be GPL-3. cogent3 is BSD and imposes nothing.

- Tool functions are synchronous. The bootstrap is CPU-bound and IQ-TREE holds
  process-global state; serialising calls keeps concurrent invocations from
  interleaving inside the engine.

- Threads are pinned to 1 before piqtree is imported, since IQ-TREE reads
  `OMP_NUM_THREADS` when its pool initialises and setting it afterwards is a
  no-op that looks like it worked.

### Not included

Protein and codon models, aLRT / approximate-Bayes / UFBoot support measures,
rooting and divergence-time estimation, ancestral state reconstruction,
partitioned models, and tree rearrangement tests.
