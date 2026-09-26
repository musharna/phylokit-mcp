# phylokit-mcp

[![ci](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![python](https://img.shields.io/pypi/pyversions/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![license](https://img.shields.io/pypi/l/phylokit-mcp)](LICENSE)
[![Glama](https://glama.ai/mcp/servers/musharna/phylokit-mcp/badges/score.svg)](https://glama.ai/mcp/servers/musharna/phylokit-mcp)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21713870.svg)](https://doi.org/10.5281/zenodo.21713870)

<!-- mcp-name: io.github.musharna/phylokit-mcp -->

Phylogenetic inference over MCP, driving IQ-TREE 3 through
[piqtree](https://github.com/cogent3/piqtree) 0.8, which builds IQ-TREE 3 into
its wheel.

**A topology without support is not a result.** `infer_tree` always runs a
bootstrap and always returns per-clade support. There is no flag to skip it.

6 tools. The test suite runs real IQ-TREE and real MAFFT (no mocked engine),
and includes 23 mutation checks and a real-process JSON-RPC handshake test.

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
  above, simulated under **JC**, the AIC winner is **not JC** (TPM2u and TPM2
  tie at the top), and several other models, JC among them at ΔAIC 0.7, sit
  inside the conventional ±2 indistinguishability margin (`seed=1`, piqtree
  0.8.3). A winner without its margin is a claim the numbers do not support.
- **Length versus evidence** — `n_parsimony_informative` alongside `n_sites`.
  A 10,000-site alignment of near-identical sequences supports nothing.

## Tools

| tool                        | what it does                                                                       |
| --------------------------- | ---------------------------------------------------------------------------------- |
| `infer_tree`                | ML tree **plus** bootstrap support, per clade. Never one without the other.        |
| `select_substitution_model` | Ranks 100+ models with ΔAIC/AICc/BIC, and says when the criteria disagree.         |
| `compare_trees`             | Robinson–Foulds distance and the clades that differ. Compares splits, not strings. |
| `simulate_alignment`        | Generates sequences along a tree you specify — the positive control.               |
| `align_sequences`           | Aligns unaligned FASTA with MAFFT; the output goes straight into `infer_tree`.     |
| `capabilities`              | piqtree and IQ-TREE versions, 215 substitution models, MAFFT version, limits.      |

## Install

```bash
pip install phylokit-mcp
```

piqtree ships prebuilt wheels, so there is no compiler, no R and no conda step —
but it requires **Python 3.12+**, and so does this package.

`align_sequences` is the one tool that needs something pip cannot install: the
[MAFFT](https://mafft.cbrc.jp/alignment/software/) binary on `PATH`
(`apt install mafft`, `brew install mafft`, or `conda install -c bioconda mafft`).
The other five tools work without it, `capabilities` reports
`aligner_version: null`, and calling `align_sequences` returns a refusal that
names the install rather than a crash. A MAFFT that is installed but does not
answer `--version` is a different state: `aligner_version` is still `null` and
`aligner_error` says what happened.

## Configure your MCP client

```json
{
  "mcpServers": {
    "phylokit": {
      "command": "uvx",
      "args": ["phylokit-mcp"]
    }
  }
}
```

`uvx` fetches the released package on demand, so this needs no prior install — but
it must resolve a **Python 3.12+** interpreter, since that is piqtree's wheel floor.
If `uvx` picks an older one, pin it with `"args": ["--python", "3.12", "phylokit-mcp"]`.

If you installed it yourself instead, `"command": "phylokit-mcp"` works when the
executable is on your `PATH`; give the absolute path to the entry point in the
environment you installed into if it is not.

The same file ships as [`.mcp.json`](.mcp.json) in this repo, which Claude Code
picks up automatically when the repo is your working directory.

## Reproducibility, stated precisely

Measured, not assumed:

- **Not bit-exact, in any setting.** The same request with the same seed — on
  repeat in one process, or in a fresh process — can return branch lengths and a
  log-likelihood that differ in the trailing digits (eight fresh processes gave
  five distinct log-likelihoods, spread ~2e-6). IQ-TREE reads the wall clock
  during its search: freezing `gettimeofday()` alone made every run
  bit-identical. piqtree exposes no option to take the clock out, so the server
  reports `deterministic_across_processes: false` rather than promise it.
- **Support can move by a replicate flipping.** Over six repeated 50-replicate
  calls, three of four clades were bit-identical and one moved **0.02**, well
  inside the bootstrap's own sampling error (~0.07 at 50 replicates). The
  topology and every conclusion were unchanged. The column resampling itself is
  numpy-seeded and exact.

Compare trees with `compare_trees`, and numbers with a tolerance — never by
string equality.

Threads are pinned to 1 before piqtree is imported: likelihood sums accumulate in
thread-completion order, floating-point addition is not associative, and
near-tied topologies can flip on the last bits. Pinning is necessary, not
sufficient.

## Limitations

- **Nucleotide and protein alignments.** Pass `sequence_type="protein"` and a
  protein model (`LG`, `WAG`, …). Codon models are still not exposed.
  The molecule type is **declared, never sniffed**: an alignment of only A/C/G/T
  is a valid protein alignment too (Ala/Cys/Gly/Thr), so guessing would fit a
  nucleotide model to protein data and return a tree, a likelihood and support
  values that are all wrong and none of which complain.
- **Bootstrap only** — no aLRT, no approximate Bayes, no UFBoot. Support is the
  nonparametric bootstrap (Felsenstein 1985), computed here rather than read back
  from IQ-TREE, because piqtree 0.8.3 runs `bootstrap_replicates` but does not
  expose the resulting values.
- **Cost is linear in replicates.** Each replicate is a full maximum-likelihood
  search on a resampled alignment, so it grows with taxon count, alignment
  length and model complexity. Capped at 200 taxa and 1000 replicates.
- **Alignment is MAFFT `--auto`, single-threaded, and nothing else.** No choice
  of strategy, no profile alignment, no trimming, at most 200 sequences of
  100,000 residues, and a 600 s wall-clock cap. Input that already contains gaps
  is refused rather than silently degapped. The tree tools still refuse ragged
  input; they do not align it for you.
- **Unrooted trees.** No rooting, no dating, no ancestral reconstruction.

## Licence

**GPL-2.0-only.** The "only" is load-bearing: piqtree declares `GPL-2.0-only`,
which is _incompatible_ with GPL-3.0, so the distributed combination cannot be
GPL-3. cogent3 is BSD and imposes nothing.

Unofficial. Not affiliated with, endorsed by, or sponsored by the IQ-TREE authors
or the cogent3 project. **IQ-TREE is academic software and expects to be cited**
— if results from this server appear in published work, cite IQ-TREE 3 as
directed at [iqtree.org](http://www.iqtree.org/) (currently Wong et al. 2026,
[doi:10.1093/molbev/msag117](https://doi.org/10.1093/molbev/msag117), plus the
paper for any method used, such as ModelFinder), and piqtree, McArthur et al.
2026, [doi:10.1093/molbev/msag061](https://doi.org/10.1093/molbev/msag061) —
not this wrapper. The same holds for MAFFT when `align_sequences` produced the
alignment: Katoh & Standley 2013,
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010). MAFFT is
BSD-licensed and is run as a separate program, not linked. See [NOTICE](NOTICE).
