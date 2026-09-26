# Security Policy

## Supported versions

`phylokit-mcp` ships fixes against the latest released version only. The current
release is **v0.4.0**. Please reproduce any issue on the latest release
(`uvx phylokit-mcp` always pulls it) before reporting.

| Version         | Supported          |
| --------------- | ------------------ |
| latest (0.4.x) | :white_check_mark: |
| < latest        | :x:                |

## Reporting a vulnerability

**Please do not open a public issue for a security vulnerability.**

Report privately with GitHub's **"Report a vulnerability"** button under the
repo's **Security** tab (private security advisories).

Please include a description of the issue, the affected version, and a minimal
reproduction. You can expect an initial acknowledgement within a few days. Once a
fix ships, you'll be credited in the release notes unless you ask otherwise.

## Security model

This server runs **phylogenetic inference in-process** through `piqtree`, which
embeds IQ-TREE 3 as a library.

- **No shell.** IQ-TREE is invoked through piqtree's Python bindings, not by
  building a command line. The one subprocess is MAFFT, for `align_sequences`:
  it is started from an argument list (never through a shell) with a fixed set
  of options, reads the sequences from a file the server writes in a private
  temporary directory, has stdin closed, and is killed at a timeout.
- **Alignments arrive as data, not as paths.** The tools take sequence data
  inline; the server does not open caller-supplied file paths, so it cannot be
  steered into reading arbitrary files off the host.
- **No network access.** The server performs no outbound requests.
- **Untrusted alignments are parsed by cogent3.** Malformed input is a parser
  concern, and parse failures surface as errors rather than being swallowed.

The realistic risk here is resource exhaustion, not code execution: maximum-likelihood
inference is CPU-bound and superlinear in taxa and sites, and nothing in the server
caps the size of an alignment a caller may submit. A large alignment can occupy a
core for a long time. Run it where that is acceptable, and prefer a bounded
deployment if the model driving it is exposed to untrusted input.
