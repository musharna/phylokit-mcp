# Security Policy

## Supported versions

`phylokit-mcp` ships fixes against the latest release on PyPI only. Please
reproduce any issue on the latest release (`uvx phylokit-mcp` always pulls it)
before reporting.

| Version        | Supported          |
| -------------- | ------------------ |
| latest release | :white_check_mark: |
| older releases | :x:                |

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
- **File paths are opened as given unless you set roots.** Every tool that
  takes FASTA or Newick also takes a path the server opens itself
  (`fasta_path`, `trees_path`), and `align_sequences` can write its result to
  `out_path`. By default a relative path resolves against the server's working
  directory, but nothing keeps it there: `../`, an absolute path, `~` and a
  symlinked directory all reach any file the server's user can read, or any
  directory it can write to.
  - **`--root DIR` confines them.** Start the server as
    `phylokit-mcp --root /data/alignments` (repeat the flag for more
    directories, or set `PHYLOKIT_ROOTS`, separated by `:`; `;` on Windows).
    Every `fasta_path` and `trees_path` must then be a file inside a root, and
    `out_path` must create its file inside one. `..` and symlinks are followed
    before the path is judged, so a symlink inside a root that points outside
    is refused. The refusal names the roots and is the same whether or not the
    outside file exists. The path is judged twice: by name before anything is
    opened, and again on the descriptor the server opened, so a path re-pointed
    between the two is refused before a byte is read or a file is created.
    The second check needs Linux or macOS; elsewhere a server with roots set
    refuses every path. Roots cover the paths a caller passes. MAFFT's private
    temporary directory is the server's own and is not inside them.
  - A file that parses as FASTA or Newick is returned as data: its sequence
    names and residues (`map_positions` echoes both) or its taxon names. A
    file that does not parse is refused, and the refusal can still carry
    details of it, such as sequence names and lengths or the set of
    unexpected characters. A FASTA file that does not start with `>` is refused
    without quoting it.
  - `out_path` only creates a new file. It refuses any name that already
    exists, including a symlink at that name, so it cannot overwrite or follow
    one; it never creates a directory.
  - If the model driving the server may act on untrusted input, set `--root`,
    and run the server as a user that can read and write only what it should,
    or in a container.
- **No network access.** The server performs no outbound requests.
- **Untrusted alignments are parsed by cogent3.** Malformed input is a parser
  concern, and parse failures surface as errors rather than being swallowed.

The other realistic risk is resource exhaustion, not code execution. Inputs are
capped: an alignment at 200 taxa and 100,000 sites, a file read at 32 MiB,
bootstrap replicates at 1,000 (UFBoot2 at 10,000), and MAFFT is killed after 600 s.
Inference itself has no time limit, and maximum-likelihood search is CPU-bound and
superlinear in taxa and sites, so an alignment inside the caps can still occupy a
core for a long time. Run it where that is acceptable, and prefer a bounded
deployment if the model driving it is exposed to untrusted input.
