# Agent benchmark: phylogenetics questions with and without the server

Does Claude answer phylogenetics questions better with phylokit-mcp than
without it? 46 questions, asked 3 times each to Claude Opus 5.5 in Claude Code
under three setups:

- **Claude alone**: no tools.
- **Claude + Python, MAFFT, web**: what a user has without this server. Bash,
  Read, WebSearch and WebFetch, a Python with piqtree, cogent3 and numpy (the
  server's own locked versions, but not the server), and MAFFT on `PATH`.
- **Claude + phylokit-mcp**: this server (0.6.0 from PyPI) only, with MAFFT on
  `PATH`.

TODO: run date and Claude Code version (the full sweep has not run yet).

![Right, wrong or declined per setup](agent_qa.png)

| Setup | Right | Wrong | Declined | Questions right in all 3 runs | Cost / answer, median (mean) | Seconds / answer, median (mean) | Slowest |
|---|---|---|---|---|---|---|---|
| Claude alone | TODO | TODO | TODO | TODO of 46 | TODO | TODO | TODO |
| Claude + Python, MAFFT, web | TODO | TODO | TODO | TODO of 46 | TODO | TODO | TODO |
| Claude + phylokit-mcp | TODO | TODO | TODO | TODO of 46 | TODO | TODO | TODO |

Out of 138 answers per setup. The three runs of a question are not
independent, so the questions column is the fairer count. The numbers come
from [`summary.tsv`](summary.tsv), which `plot.R` writes from
[`results.jsonl`](results.jsonl).

## What it shows

TODO after the sweep. Check before writing: per-question counts, not only
totals; medians beside means; and for any label put on a wrong answer, the
transcript that shows it.

![Cost and time per answer, every run](agent_qa_cost_time.png)

![Right, wrong or declined by kind of question](agent_qa_by_category.png)

## What it does not show

- **The data are simulated.** Trees are drawn at random and sequences evolved
  down them under JC69 or the HKY85 family (indels added for the alignment
  questions). Real alignments have rate variation, compositional shifts and
  alignment error that none of these questions contain.
- **The questions are answerable on purpose.** A draw was kept only when
  neighbour joining and piqtree's ML tree recovered the generating tree
  (sister questions), MAFFT put the asked residue in the true column
  (homology), the generating model won BIC by at least 2 (model choice), and
  both distance estimates fell in the tolerance (distances). Hard cases, where
  the data do not support the answer, are not tested, and those are where the
  server's bootstrap support matters most.
- **Support values are not scored.** No question asks how well supported a
  clade is.
- **One model, one date.**

## Questions

[`questions.json`](questions.json) holds the 46 questions with their answers,
built on 2026-10-16 with seed 20261016 by `build_questions.py`. The data each
question needs is in the prompt and also written as a file
(`alignment.fasta`, `sequences.fasta` or `trees.nwk`) into the run's empty
working directory.

| Kind | Count | Server tool it exercises | Answer key |
|---|---|---|---|
| Sister taxon in the ML tree (8 taxa, 600 sites, JC69) | 8 | `infer_tree` | the asked taxon's cherry partner in the generating tree |
| Robinson-Foulds distance (8-12 taxa, 0-7 NNI moves apart) | 8 | `compare_trees` | splits of the two Newick strings as written, symmetric difference |
| Homologous residue after alignment (6 seqs, 4 DNA + 4 protein) | 8 | `align_sequences` | the true alignment from an indel simulator that tracks every residue's column |
| Model with the lowest BIC among JC69, K80, F81, HKY85 (6 taxa, 1000 sites) | 8 | `select_substitution_model` | the generating model, 2 per model |
| JC69 distance between two taxa (6 taxa, 1000 sites) | 8 | `infer_tree` (branch lengths) | path length between the two tips in the generating tree |
| Parsimony-informative sites (6-10 taxa, 150-300 sites) | 6 | `infer_tree`, `align_sequences` (`n_parsimony_informative`) | counted from the alignment |

The key comes from [`truth.py`](truth.py): numpy and plain Python, no
`phylokit_mcp`, cogent3 or piqtree, so a bug in the server or the libraries
under it cannot also be in the key (`tests/test_agent_qa_score.py` checks the
imports). piqtree and MAFFT are used only to check that a draw is answerable;
1 of the draws was discarded (a homology question where MAFFT disagreed), and
`questions.json` records the count. Rebuilding gives the same questions and
answers (checked by building twice and comparing).

## Scoring

Each reply ends with `ANSWER: <value>` or `ANSWER: UNKNOWN`. `score.py`
decides the verdict from that line alone:

- **right**: the value matches the key.
- **wrong**: the value does not match.
- **declined**: `UNKNOWN`.
- **no answer**: there is no `ANSWER:` line.

Matching rules by kind, and why each is as strict as it is:

- **Sister taxon**: the taxon name, without case. The draw was checked to give
  the generating tree under NJ and ML, so there is one right name.
- **RF distance, informative sites**: the integer, exactly. Both are
  deterministic functions of the text in the question. "0.6" is not read as
  6.
- **Homologous residue**: the integer, exactly. Positions have no gap in any
  sequence within 8 columns, so the homology is not in doubt.
- **Model**: any name of the right model. IQ-TREE writes HKY85 as `HKY+F` and
  F81 as `F81+F`; `K2P+F` is HKY85, not K80.
- **Distance**: within 3 large-sample SDs of the JC69 distance estimate at
  the true distance, sd = sqrt(p(1-p)/n) / (1 - 4p/3). An estimate from data
  differs from the generating value by sampling error, so an exact match would
  score correct work as wrong. True distances run from 0.05 to 1.6; the best
  single guessed number would be right on 3 of 8.

A wrong answer that is still a legitimate option (a taxon in the alignment,
one of the four models) is also flagged `wrong_id`. Each row also records,
read from the transcript, whether any tool call named `phylokit` or the
benchmark's key file (`mentions_phylokit`, `mentions_key`).
Tests: `tests/test_agent_qa_score.py`.

## Reproduce

```bash
export MAFFT=/path/to/mafft                                 # if MAFFT is not on PATH
uv run python benchmarks/agent_qa/build_questions.py        # rebuild the key (writes questions.json)
uv run python benchmarks/agent_qa/run.py --setup-alt        # once: the alt arm's Python, without the server
uv run python benchmarks/agent_qa/run.py --reps 3 --workers 2   # 414 claude -p runs; resumable
uv run python benchmarks/agent_qa/run.py --rescore          # re-score stored replies after a key fix
Rscript benchmarks/agent_qa/plot.R                          # summary.tsv + figures
```

Each run is an isolated `claude -p`:
- `--setting-sources ""` and an empty working directory holding only the
  question's data file.
- `--strict-mcp-config`, so no other MCP servers.
- The same short system prompt in every setup.

Raw transcripts go to `runs/`, which is gitignored. A pilot of 2 questions x
3 setups (q01 sister taxon, q17 homologous residue) ran end to end before the
sweep: every run exited 0 with the expected tools and server, and every reply
scored.
