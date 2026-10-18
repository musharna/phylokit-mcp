# Install and client setup, in detail

The README has the short version. This page holds the rest.

## The engine

IQ-TREE 3 comes through [piqtree](https://github.com/cogent3/piqtree) 0.8,
which builds it into its wheel. piqtree ships prebuilt wheels for Python 3.12+,
so there is no compiler, R or conda step.

## Without MAFFT

`align_sequences` needs the [MAFFT](https://mafft.cbrc.jp/alignment/software/)
binary on `PATH` (`apt install mafft`, `brew install mafft` or
`conda install -c bioconda mafft`). Without it the other five tools work,
`capabilities` reports `aligner_version: null`, and `align_sequences` returns a
refusal naming the install rather than crashing. If MAFFT is installed but does
not answer `--version`, `aligner_error` says what happened.

## Client configuration

`uvx` fetches the released package on demand but must resolve a Python 3.12+
interpreter; if it picks an older one, use
`"args": ["--python", "3.12", "phylokit-mcp"]`. For your own install, use
`"command": "phylokit-mcp"`, or the entry point's absolute path if it is not on
`PATH`. The same config ships as [`.mcp.json`](../.mcp.json), which Claude Code
picks up when this repo is the working directory.
