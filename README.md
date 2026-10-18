# phylokit-mcp

phylokit-mcp is an MCP server that lets AI agents build maximum-likelihood
trees with IQ-TREE 3 and align sequences with MAFFT.

[![CI](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/musharna/phylokit-mcp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![Python](https://img.shields.io/pypi/pyversions/phylokit-mcp)](https://pypi.org/project/phylokit-mcp/)
[![License](https://img.shields.io/pypi/l/phylokit-mcp)](https://github.com/musharna/phylokit-mcp/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21713870.svg)](https://doi.org/10.5281/zenodo.21713870)
[![Glama](https://glama.ai/mcp/servers/musharna/phylokit-mcp/badges/score.svg)](https://glama.ai/mcp/servers/musharna/phylokit-mcp)

<!-- mcp-name: io.github.musharna/phylokit-mcp -->

**A topology without support is not a result.** `infer_tree` always runs a
bootstrap and returns per-clade support. There is no flag to skip it.

## Why the rule

A maximum-likelihood tree looks the same whether or not the data support it.
Measured on alignments simulated from a _known_ 7-taxon tree (JC, 100
replicates, seed 1):

| Sites | Informative sites | Recovered the true tree? | Lowest clade support |
| ----: | ----------------: | ------------------------ | -------------------: |
|   300 |                51 | yes, exactly             |                 1.00 |
|    60 |                11 | **no, RF 2**             |                 0.57 |

At 60 sites the tree contains a clade (`C,D,G`) that does not exist and omits
one that does (`E,F,G`). Both runs return a fully resolved Newick string of the
same shape. Only the support values tell them apart, and the false clade is the
_lowest-supported_ one in the tree.

## Does it help?

Not in our benchmark: with this server Claude was no more accurate than with a
Python shell, and it was slower and cost more.

<p align="center">
  <img src="https://raw.githubusercontent.com/musharna/phylokit-mcp/main/benchmarks/agent_qa/agent_qa.png" alt="Bar chart of 114 answers per setup (38 questions, 3 runs each). Claude alone: 114 declined. Read only: 76 right, 24 wrong, 12 declined, 2 invalid. Python, MAFFT and web: 113 right, 1 invalid. phylokit-mcp: 110 right, 4 wrong." width="780">
</p>

Claude Opus 5.5 answered 38 questions 3 times each. Counting answers, as the
chart does, the server got 110 of 114 right and a Python shell with piqtree and
MAFFT got 113; counting questions, 36 of 38 were right in every run against 37
(sign test p = 1). Per question the server took 6.0 times as long (95% CI 4.4
to 8.3) and cost 14.4 times as much (9.9 to 20.2), each a geometric mean of
per-question ratios. A Read-only control got most answers right on half the
kinds of question by reading the file alone.

The benchmark does not score support values; method, controls and other
limitations are in [benchmarks/agent_qa](https://github.com/musharna/phylokit-mcp/blob/main/benchmarks/agent_qa/README.md).

## What it reports that a Newick string cannot

- **Conflicting clades**: groupings the data support at ≥0.70 that are _absent_
  from the reported tree. A support-annotated Newick string has nowhere to put
  them.
- **`fraction_resolved`**: the share of clades clearing 0.70. The headline
  number, before any single grouping is repeated as fact.
- **Model runners-up with ΔAIC**, not just a winner. On the 300-site alignment
  above, simulated under **JC**, the AIC winner is **not JC** (TPM2u and TPM2
  tie), and several models, JC among them at ΔAIC 0.7, sit inside the
  conventional ±2 margin (`seed=1`, piqtree 0.8.3).
- **Length versus evidence**: `n_parsimony_informative` alongside `n_sites`. A
  10,000-site alignment of near-identical sequences supports nothing.

## Tools

| Tool                        | What it does                                                                       |
| --------------------------- | ---------------------------------------------------------------------------------- |
| `infer_tree`                | ML tree **plus** bootstrap support, per clade.                                     |
| `select_substitution_model` | Ranks 100+ models with ΔAIC/AICc/BIC, and says when the criteria disagree.         |
| `compare_trees`             | Robinson–Foulds distance and the clades that differ. Compares splits, not strings. |
| `simulate_alignment`        | Generates sequences along a tree you specify: the positive control.                |
| `align_sequences`           | Aligns unaligned FASTA with MAFFT; the output goes straight into `infer_tree`.     |
| `capabilities`              | piqtree and IQ-TREE versions, 215 substitution models, MAFFT version, limits.      |

The tests run real IQ-TREE and real MAFFT (no mocked engine), with 23 mutation
checks and a real-process JSON-RPC handshake test.

## Install

```bash
pip install phylokit-mcp
```

Requires **Python 3.12+**. IQ-TREE 3 comes through
[piqtree](https://github.com/cogent3/piqtree) 0.8, which builds it into its
wheel. `align_sequences` also needs the
[MAFFT](https://mafft.cbrc.jp/alignment/software/) binary on `PATH`; the other
five tools work without it. Details: [docs/INSTALL.md](https://github.com/musharna/phylokit-mcp/blob/main/docs/INSTALL.md).

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

If `uvx` picks a Python older than 3.12, use
`"args": ["--python", "3.12", "phylokit-mcp"]`.

## Reproducibility

Not bit-exact, in any setting: IQ-TREE reads the wall clock during its search,
so the same request and seed can differ in the trailing digits of branch
lengths and log-likelihood, and a support value can move by a replicate
flipping. The server reports `deterministic_across_processes: false`. Compare
trees with `compare_trees` and numbers with a tolerance, never by string
equality. Measurements: [docs/REPRODUCIBILITY.md](https://github.com/musharna/phylokit-mcp/blob/main/docs/REPRODUCIBILITY.md).

## Limitations

- **Nucleotide and protein, no codon models.** Pass `sequence_type="protein"`
  and a protein model (`LG`, `WAG`, …). The molecule type is **declared, never
  sniffed**: an A/C/G/T-only alignment is also valid protein (Ala/Cys/Gly/Thr),
  and a wrong guess returns a wrong tree, likelihood and support with no
  warning.
- **Bootstrap only**: no aLRT, approximate Bayes or UFBoot. Support is the
  nonparametric bootstrap (Felsenstein 1985), computed here because piqtree
  0.8.3 runs `bootstrap_replicates` but does not expose the values.
- **Cost is linear in replicates**, each a full ML search. Capped at 200 taxa
  and 1000 replicates.
- **Alignment is MAFFT `--auto`, single-threaded, and nothing else**: no
  strategy choice, profile alignment or trimming; at most 200 sequences of
  100,000 residues and 600 s. Gapped input is refused, not degapped. The tree
  tools refuse ragged input rather than align it.
- **Unrooted trees.** No rooting, dating or ancestral reconstruction.

## License

**GPL-2.0-only**, because piqtree is `GPL-2.0-only`, which is _incompatible_
with GPL-3.0. cogent3 is BSD.

If you publish results from this server, cite IQ-TREE 3, piqtree and (when it
made the alignment) MAFFT, not this wrapper. The references are in
[CITATION.cff](https://github.com/musharna/phylokit-mcp/blob/main/CITATION.cff),
and licensing and affiliation are in
[NOTICE](https://github.com/musharna/phylokit-mcp/blob/main/NOTICE).
