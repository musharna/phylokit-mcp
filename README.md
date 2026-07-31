# phylokit-mcp

[![ci](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![python](https://img.shields.io/pypi/pyversions/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![license](https://img.shields.io/pypi/l/phylokit-mcp)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21713870.svg)](https://doi.org/10.5281/zenodo.21713870)

<!-- mcp-name: io.github.musharna/phylokit-mcp -->

Phylogenetic inference over MCP, driving IQ-TREE 2 through
[piqtree](https://github.com/cogent3/piqtree).

**A topology without support is not a result.** `infer_tree` always runs a
bootstrap and always returns per-clade support. There is no flag to skip it.

5 tools, 62 tests against real IQ-TREE (no mocked engine), 10 mutation checks,
and a real-process JSON-RPC handshake test.

## Why the rule

A maximum-likelihood tree looks identical whether or not the data support it.
Measured here, on alignments simulated from a _known_ 7-taxon tree so the right
answer is not in doubt:

| sites | informative sites | recovered the true tree? | lowest clade support |
| ----: | ----------------: | ------------------------ | -------------------: |
|   300 |                51 | yes, exactly             |                 1.00 |
|    60 |                11 | **no — RF 2**            |                 0.57 |

At 60 sites the tree contains a clade (`C,D,G`) that does not exist and omits
one that does (`E,F,G`). Both runs return a fully resolved Newick string of the
same shape; nothing about the topology itself distinguishes them. The support
values do — and the false clade is the _lowest-supported_ one in the tree.

That is the entire argument for this server. Returning a bare tree returns a
result the caller cannot evaluate.

## What it reports that a Newick string cannot

- **Conflicting clades** — groupings the data support at ≥0.70 that are _absent_
  from the reported tree. A support-annotated Newick string has nowhere to
  attach these, so the standard format silently drops them.
- **`fraction_resolved`** — the share of clades clearing 0.70. The headline
  number, before any individual grouping is repeated as fact.
- **Model runners-up with ΔAIC** — not just a winner. On the 300-site alignment
  above, simulated under **JC**, the AIC winner is **F81**, with several models
  inside the conventional ±2 indistinguishability margin. A winner without its
  margin is a claim the numbers do not support.
- **Length versus evidence** — `n_parsimony_informative` alongside `n_sites`.
  A 10,000-site alignment of near-identical sequences supports nothing.

## Tools

| tool                        | what it does                                                                       |
| --------------------------- | ---------------------------------------------------------------------------------- |
| `infer_tree`                | ML tree **plus** bootstrap support, per clade. Never one without the other.        |
| `select_substitution_model` | Ranks 100+ models with ΔAIC/AICc/BIC, and says when the criteria disagree.         |
| `compare_trees`             | Robinson–Foulds distance and the clades that differ. Compares splits, not strings. |
| `simulate_alignment`        | Generates sequences along a tree you specify — the positive control.               |
| `capabilities`              | Engine version, 215 substitution models, enforced limits.                          |

## Install

```bash
pip install phylokit-mcp
```

piqtree ships prebuilt wheels, so there is no compiler, no R and no conda step —
but it requires **Python 3.12+**, and so does this package.

## Reproducibility, stated precisely

Measured, not assumed:

- **Across fresh processes: exact.** Three runs of an identical 30-replicate
  bootstrap returned byte-identical support.
- **Within one long-lived process: not bit-exact.** Passing the same `rand_seed`
  does not fully reset IQ-TREE's internal state — building the same tree three
  times gave call 1 == call 2 but call 3 different.

The practical size: over six repeated 50-replicate calls, three of four clades
were bit-identical and one moved **0.02** — a single replicate flipping, well
inside the bootstrap's own sampling error (~0.07 at 50 replicates). The topology
and every conclusion were unchanged. This is reported in every response rather
than papered over, because an MCP server is long-lived by design and that is
exactly the condition which exposes it.

Threads are pinned to 1 before piqtree is imported: likelihood sums accumulate in
thread-completion order, floating-point addition is not associative, and
near-tied topologies can flip on the last bits.

## Limitations

- **Nucleotide alignments only.** Protein and codon models are not exposed.
- **Bootstrap only** — no aLRT, no approximate Bayes, no UFBoot. Support is the
  nonparametric bootstrap (Felsenstein 1985), computed here rather than read back
  from IQ-TREE, because piqtree 0.8.3 runs `bootstrap_replicates` but does not
  expose the resulting values.
- **Cost is linear in replicates.** ~130 ms per replicate at 7 taxa / 300 sites,
  and it grows with taxon count. Capped at 200 taxa and 1000 replicates.
- **It does not align sequences.** Ragged input is refused, not guessed at.
- **Unrooted trees.** No rooting, no dating, no ancestral reconstruction.

## Licence

**GPL-2.0-only.** The "only" is load-bearing: piqtree declares `GPL-2.0-only`,
which is _incompatible_ with GPL-3.0, so the distributed combination cannot be
GPL-3. cogent3 is BSD and imposes nothing.

Unofficial. Not affiliated with, endorsed by, or sponsored by the IQ-TREE authors
or the cogent3 project. **IQ-TREE 2 is academic software and expects to be cited**
— if results from this server appear in published work, cite IQ-TREE 2 as
directed at [iqtree.org](http://www.iqtree.org/), not this wrapper. See
[NOTICE](NOTICE).
