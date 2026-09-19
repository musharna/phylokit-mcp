# Changelog

All notable changes to `phylokit-mcp` are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] — 2026-09-18

### Added

- **`align_sequences`: the step before `infer_tree`.** The server refused ragged
  input and told the caller to align it somewhere else. It now runs MAFFT
  (`--auto --thread 1 --preservecase --inputorder`, with `--nuc` or `--amino`
  passed from the declared `sequence_type` so MAFFT never sniffs the molecule
  type either) and returns FASTA that `infer_tree` and
  `select_substitution_model` accept unchanged, with the same alignment summary
  and advisories the other tools return.
  - MAFFT is a system binary, not a wheel, so it is found on `PATH` at call
    time. Without it the five tree tools work as before, `capabilities` reports
    `aligner_version: null`, and `align_sequences` returns a refusal naming the
    install command. An installed MAFFT that prints nothing or hangs on
    `--version` is reported in `aligner_error` instead of failing `capabilities`. CI installs it through `scripts/guardrails-setup.sh`, the
    hook the template's guardrails workflows already look for, so `ci.yml` and
    `guardrails.yml` share one install step.
  - Sequences reach MAFFT through a file in a private temporary directory and
    the command is an argument list; nothing caller-supplied is on a command
    line. 2-200 sequences of at most 100,000 residues (the existing caps), a
    600 s wall-clock cap reported as a refusal, names held to the Newick-safe
    set the tree tools require.
  - Two things MAFFT does silently are checked rather than trusted. It exits 0
    on a single sequence and returns it unchanged, so the count check is made
    here. And it strips gap characters and lower-cases nucleotides in its input
    without saying so, so gapped input is refused up front and every output
    row, degapped, must equal its input sequence — with names and order
    unchanged — or the call fails.
  - The known-answer control plants a 12-base deletion and a 9-base insertion
    in a 120-base ancestor, chosen so neither indel can slide a column at equal
    score, and asserts the gaps come back in exactly those rows and columns. End
    to end, sequences cut from a simulated alignment are refused by `infer_tree`,
    aligned, accepted, and recover the true topology. Replacing MAFFT's output
    with right-padding — the cheapest thing that produces equal-length rows —
    fails both: the columns are wrong and the tree comes back at RF 4.
    Thirteen mutants added to `docs/mutate.sh` (11-23), all killed.
  - A two-sequence alignment is returned, flagged `ready_for_infer_tree: false`
    with a `too_few_taxa_for_a_tree` advisory, since the tree tools need four.
  - A MAFFT crash (non-zero exit, output that is not FASTA, altered residues) is
    a `RuntimeError` carrying MAFFT's stderr. Under the refusal boundary that is
    a crash, so it is masked to the model and lands in the server log, as the
    SDK intends for anything that is not an anticipated refusal.
  - Cite MAFFT when this tool produced the alignment: Katoh & Standley 2013,
    doi:10.1093/molbev/mst010 (added to `CITATION.cff` references).

### Changed

- The ragged-alignment refusal now points at `align_sequences` instead of
  saying the server does not align.
- `docs/mutate.sh` sets `PYTHONDONTWRITEBYTECODE=1`: a same-size source swap can
  leave bytecode that still runs as the mutant after the restore.

- **mypy is declared, configured and gated in CI.** It was none of those things,
  which meant `uv run mypy src` picked up an ambient install that could not
  resolve this project's dependencies, and reported eight failures that were
  mostly artefacts of running it wrong — `mcp` "not found" while importing fine
  at runtime, and `ToolAnnotations(read_only_hint=...)` flagged as a bad keyword.
  That last one is correct code: the MCP types set a camelCase `alias_generator`
  with `populate_by_name`, so the field IS snake_case and only the wire format is
  camelCase. Without the pydantic plugin mypy reads only the alias signature, and
  "fixing" the source to match would have pushed a wire spelling into Python to
  satisfy a checker that was mis-modelling it. Now: mypy in the dev group,
  `[tool.mypy]` with the pydantic plugin, `python_version` matching
  `requires-python` (3.12 — checking at 3.11 made mypy fail to *parse* piqtree,
  which uses PEP 695 syntax), and `ignore_missing_imports` scoped to `cogent3.*`
  alone rather than set globally, so it cannot swallow a real resolution failure
  the way it would have swallowed the `mcp` one.

### Fixed

- **Removed a stale `type: ignore[typeddict-item]`** that `warn_unused_ignores`
  surfaced. The four others in the file were each checked by removal and are
  load-bearing — dropping any one produces an error, so this was the only dead
  suppression rather than the tidiest-looking one.

- **Refusals reach the calling agent again under mcp >= 2.1.** mcp 2.1.0
  (python-sdk #3314) treats any exception other than `ToolError` as a crash:
  the model sees only `Error executing tool <name>` and the reason stays in the
  server log. Every refusal this server raises on purpose is a `ValueError` — a
  ragged alignment, an unknown criterion, too few replicates — and its message
  is the product; under 2.0 the text went through regardless, so nothing in the
  server said so. Each tool is now wrapped at registration to re-raise those as
  `ToolError`. The tool functions keep raising their own types for the unit
  tests that import them directly, and a genuine crash stays masked as the SDK
  intends. `test_a_ragged_alignment_surfaces_as_a_tool_error_not_a_crash`
  matches the refusal TEXT at the tool layer, which is what makes it able to
  fail at all: asserting only that the call errored is byte-identical whether
  the reason survived or was masked.

### Added

- **A guard that the refusal boundary covers every tool, not just the one under
  test.** The wrapper is applied by hand at each registration, and the text
  assertion above pins a single tool — it passes unchanged on the day a sixth
  tool is registered without the wrapper, which is the same bug returning with
  no failing test. `tests/test_refusal_boundary.py` reads the tool registry off
  the built server and asserts every registered callable came out of
  `_surfaces_refusals`, identified by code object rather than by `__wrapped__`
  or `__name__` — `functools.wraps` copies the wrapped function's name onto the
  wrapper, so neither of those can tell the boundary from any other decorator.
  Reading the registry rather than parsing the source keeps the guard valid for
  either registration shape. A second test asserts both directions in one place:
  a refusal converts, and a genuine `TypeError` still propagates unmasked, so an
  over-broad `except` cannot pass by destroying the SDK's crash signal.

## [0.4.0] — 2026-08-02

### Changed

- **The server is built by a `build_server()` factory instead of a module-level
  singleton.** `phylokit_mcp.server` no longer constructs an `MCPServer` at import
  time, so importing the module has no side effect and a test that needs a server
  gets its own. This matches `breedsim-mcp` and `plantcv-mcp`, which already did
  it; this package and `ldraw-mcp` were the outliers.

  **The tool functions stay at module level** and are registered inside the
  factory. The sibling servers nest their tool definitions instead, which they
  can do because their tools are thin wrappers over logic living in other
  modules. Here the tool functions *are* the implementation and are unit-tested
  by direct import, so nesting them would put the logic out of reach. Registering
  rather than nesting keeps both properties and collects the whole tool surface
  in one readable block.

  Two tests pin the property, because otherwise nothing stops the singleton
  returning: one asserts two calls yield distinct, **fully populated** servers,
  and one asserts the module has no `mcp` attribute. Both were seen to fail
  against deliberate mutants — a reintroduced module-level singleton, a cached
  factory, and a factory that registers no tools. That third mutant is why the
  first test checks the tool surface: `a is not b` alone passes for two empty
  servers.

  Verified from an installed wheel, not just in-repo: 5 tools registered and the
  `phylokit-mcp` console script still resolves.

- `_READ_ONLY` now uses the snake_case `ToolAnnotations` spellings, matching the
  sibling servers. **This is cosmetic, not a fix** — measured, the camelCase
  aliases build an equal object with the same `read_only_hint` attribute.

### Added

- **Community-health and repo-hygiene files, matching the standard set by
  `data-aggregator-mcp` and `plant-genomics-mcp`.** An earlier parity audit
  compared this repo only against `ldraw-mcp`, which is itself thin on these, so
  the whole tier went unnoticed: `CONTRIBUTING.md`, `SECURITY.md`, issue forms
  (bug report + feature request + a config pointing security reports at private
  advisories), a pull-request template, `.editorconfig`, `.mcp.json`, `glama.json`,
  a CodeQL workflow, and a Dependabot config.

  **Dependabot uses the `uv` ecosystem, not `pip`.** This is a uv-locked project;
  the pip ecosystem would update `pyproject.toml` and leave `uv.lock` stale, which
  CI installs with `--frozen` and would fail on. Dependabot's native uv support
  reads both together.

  `CONTRIBUTING.md` and `SECURITY.md` were added to the sdist allow-list.
  hatchling's allow-list drops anything unlisted **silently** — verified with
  `tar tzf` on a real build rather than assumed, the same way a `NOTICE` was
  previously found missing.

  `SECURITY.md` documents this server's actual trust boundary: inference runs
  in-process through piqtree with no subprocess and no shell, alignments arrive as
  inline data rather than filesystem paths, and the realistic risk is resource
  exhaustion rather than code execution.

- **A "Configure your MCP client" section in the README.** This server had none at
  all — the one section an MCP server's README cannot do without. It notes that
  `uvx` must resolve a Python 3.12+ interpreter, since that is piqtree's wheel floor.

- **README gained a Glama badge**, verified HTTP 200 with a bogus name as a
  negative control.

- **`project.urls` gained `Changelog`.**

- **`server.json` is validated against the registry's own published schema.**
  `breedsim-mcp` v0.4.0 was tagged, uploaded to PyPI and GitHub-released before
  the MCP registry refused it with a 422: its description had grown past a
  100-character cap that nothing local measured. The publish workflow is the only
  thing that checks registry constraints, and it runs on tag push — after the
  version is already burned. This server's description is 83 characters and would
  have passed, but nothing here was measuring it either.

  Rather than copy the one constant, the check validates the whole document
  against the dated `$schema` `server.json` already declares, which is the
  registry's own statement of what it accepts. That covers the four other length
  caps and the required-field list as well. The schema is vendored at
  `tests/server.schema.json` rather than fetched, keeping the suite offline and
  deterministic, and a test asserts the vendored copy's `$id` still matches the
  declared `$schema` so the pin cannot drift silently.

  Verified by mutation rather than assumed: an over-long description was written
  into the real `server.json` and the suite watched to fail on it, with a
  non-length failure (a missing required field) and an in-test positive control
  so a validator that raised on everything could not read as a working guard.

## [0.3.0] — 2026-08-01

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

- **Topology recovery is now checked against a known simulated tree.** Sequences
  are evolved along a tree this server never sees, and inference is asserted to
  recover it. Every prior test compared the server against itself, so a change
  that broke inference consistently would have kept them all green.

  It ships with BUILT-IN NEGATIVE CONTROLS, and they earned their place
  immediately: a blinded Robinson-Foulds comparator — one that had stopped
  looking at the trees at all — was caught by the controls while the primary
  "we recovered the tree" assertion passed happily, because RF = 0 reads as
  perfect recovery whether you compared the right things or nothing at all.

- **Every tree states the engine that built it and the units of its branch
  lengths.** `engine` carries the piqtree/IQ-TREE version and
  `branch_length_units` is `"substitutions per site"`. Newick carries bare
  numbers, and reading those as time or as percent divergence gives a confidently
  wrong answer; the units are not recoverable from the string. The engine version
  was previously reachable only via the separate `capabilities` call, so a saved
  tree was not self-describing — and a test asserts the two agree, because two
  sources that can disagree are worse than one.

### Fixed

- **A signal-free alignment now gets a diagnosis instead of an internal error.**
  Four identical sequences reached IQ-TREE, which cannot fit a likelihood and
  failed inside piqtree with `IQ-TREE output is malformed, likelihood not found`.
  That message is upstream's and describes a PARSING failure, so a caller reads it
  as this server being broken and retries — when the real answer is that the data
  support every topology equally and no retry will change that. The condition is
  knowable before the call, and is now checked there.

  The check is deliberately NOT part of `validate`. That function answers "is this
  a well-formed alignment", and a signal-free alignment is perfectly well formed;
  `select_substitution_model` and the molecule-type checks have legitimate reasons
  to accept one. Putting it in `validate` broke two existing tests whose minimal
  fixtures were correct — the check was in the wrong layer, not the tests.

- **Branch lengths pinned at the optimiser's ceiling are now flagged.** Saturated
  data drive IQ-TREE to its upper bound of 10 substitutions per site, and lengths
  of `9.9999989` were returned unremarked, indistinguishable from a measurement.
  A `saturated_branch_lengths` advisory now names them and states the number is a
  FLOOR rather than an estimate. The support warnings were already good; a
  degenerate FIT had no equivalent.

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
