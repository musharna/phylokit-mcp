# Changelog

All notable changes to `phylokit-mcp` are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Protein alignments.** `infer_tree` and `select_substitution_model` take
  `sequence_type="protein"`; piqtree has supported protein models all along
  (`available_models('protein')`) and this server refused them at the door.

  **The molecule type is declared, never inferred.** An alignment of only A/C/G/T
  is a perfectly valid protein alignment — alanine, cysteine, glycine, threonine
  — so no amount of sniffing can separate the two cases. Guessing wrong fits a
  nucleotide model to protein data and returns a tree, a likelihood and bootstrap
  support, all wrong and none complaining. A test pins that ambiguity by
  validating the same alignment successfully under BOTH types.

### Fixed

- **`parsimony_informative` counted states from a hardcoded `ACGTU`.** On a
  protein alignment that scores every site as uninformative, so the
  `thin_information` advisory would have refused valid protein data while
  appearing to have measured it. The alphabet now follows the declared molecule
  type. This was latent rather than user-visible before, because protein input
  could not get that far.

## [0.2.2] — 2026-07-31

### Added

- **Zenodo archival.** This release exists to be archived: the Zenodo↔GitHub
  integration mints a DOI from the tag's tarball, and the previous tag predated
  `.zenodo.json` and `CITATION.cff` entirely — those files were added after it was
  cut. Zenodo archives the tag, not the default branch, so a release was the only
  way to get the metadata into an archived snapshot.

### Changed

- **`.zenodo.json` now uses Zenodo's lowercase licence identifier**
  (`gpl-2.0-only` rather than `GPL-2.0-only`). That is the canonical spelling —
  `zenodo.org/api/vocabularies/licenses/<id>` returns 200 for the lowercase form
  and 404 for the SPDX-cased one. See the correction below: it fixed nothing.

### Correction — added after this release was published

This release was originally described here as **fixing** a defect in which the
SPDX casing "silently dropped the licence from the published record". **That was
wrong, and the entry is corrected rather than quietly deleted.**

Zenodo normalises the licence identifier on ingest. The sibling `ldraw-mcp`
archived with `"MIT"` still in place and its record reads `license: mit-license`;
this project's record reads `license: gpl-2.0-only`. The licence was never dropped.

The apparent evidence was two of my own measurement errors, both the same
mistake — probing a proxy instead of the artifact:

1. Querying the licence **vocabulary endpoint** and treating a 404 there as what
   the ingest accepts. It is not; the ingest normalises casing.
2. Reading the **RDM-era field names** (`rights`, `subjects`,
   `creators[].person_or_org`) against an API endpoint that returns the **legacy**
   shape (`metadata.license`, `metadata.keywords`, `creators[].orcid`). Every
   field reported as absent was present throughout.

What remains true is the reason this release exists: the previous tag predated
`.zenodo.json` and `CITATION.cff`, and Zenodo archives the **tag**, not the
default branch. DOI: [10.5281/zenodo.21713871](https://doi.org/10.5281/zenodo.21713871).

### Notes

No functional change. Tools, guards and dependency pins are identical to 0.2.1.

## [0.2.1] — 2026-07-31

### Added

- **Published to the official MCP registry** (`io.github.musharna/phylokit-mcp`)
  via `server.json` and an OIDC workflow, so the server is discoverable from MCP
  clients and directories rather than only from PyPI.

  This needed a release rather than a docs commit. The registry proves PyPI
  ownership by finding an `mcp-name` marker in the package README **as published
  to PyPI**, and PyPI captures `long_description` at release time — so a marker
  sitting on `main` verifies nothing. It is the same mechanism that kept
  `plantcv-mcp`'s "Not published to PyPI" line live on its project page after the
  fix had merged.

- **`tests/test_registry_metadata.py`.** `server.json` states the version in
  three places and nothing else makes them agree with `pyproject.toml`; a stale
  one is rejected by the registry during a release, after the version is spent.
  The README marker is checked against the name `server.json` declares, since
  that exact string is what the registry greps for.

  `server.json` also declares `OMP_NUM_THREADS` and `MKL_NUM_THREADS`, which this
  server pins for reproducibility — a caller overriding them silently loses the
  determinism the tool reports.

### Notes

No functional change. Tools, the mandatory bootstrap and dependency pins are
identical to 0.2.0.

Running the real `mcp-publisher validate` against `server.json` is what caught a
100-character cap on `description` in the sibling `breedsim-mcp` — a constraint
no schema read surfaced, and one that would otherwise have failed the publish
after the version was already on PyPI. The workflow validates before
authenticating for that reason.

## [0.2.0] — 2026-07-30

### Changed

- **Migrated to `mcp` 2.x.** `mcp.server.fastmcp` no longer exists in 2.0.0, but
  `FastMCP` was **renamed, not removed** — it is now
  `mcp.server.mcpserver.MCPServer`, the same class with the same decorator and
  the same `annotations` / `structured_output` kwargs. `ToolError` moved to
  `mcp.server.mcpserver.exceptions`; `mcp.types.ToolAnnotations` did not move.

  The dependency moves to `mcp>=2,<3` rather than widening to `<3`. This package
  imports `mcp.server.mcpserver`, which does not exist in 1.x, so a range
  spanning both majors could resolve to a version that cannot import the server.
  The old 1.28.1 floor was a security floor, not a feature one, and every 2.x
  release is above it.

- **`mcp.types` field names went camelCase → snake_case** (`inputSchema` →
  `input_schema`). This touched a test assertion, not the server.

### Fixed

- **`__version__` is read from installed metadata instead of restated.** It was a
  literal beside a `pyproject.toml` version, with nothing enforcing agreement.
  They agreed here — but `plantcv-mcp` shipped 0.2.0 reporting `"0.1.0"` from
  exactly that arrangement, so this is a latent form of a defect that has already
  shipped elsewhere.

  There was no version test at all; `tests/test_version.py` now compares the
  reported version against what `pyproject.toml` declares. Deliberately not
  against `importlib.metadata`, which is what `__version__` now reads _from_ —
  asserting those agree would compare a value to itself. Confirmed to fail on a
  reintroduced literal before being kept.

## [0.1.0] — 2026-07-30

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
