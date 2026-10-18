# Agent benchmark: phylogenetics questions with and without the server

Does Claude answer phylogenetics questions better with phylokit-mcp than
without it? 38 questions, asked 3 times each to Claude Opus 5.5 in Claude Code
under four setups:

- **Claude alone**: no tools. The data are in a file it cannot open, so the
  right reply is `UNKNOWN`. It is the floor, not a competitor.
- **Claude + Read only (control)**: it can see every byte of the data file and
  has nothing to compute with: no shell, no server, no web. It tests whether
  the questions need computation at all.
- **Claude + Python, MAFFT, web**: what a user has without this server. Bash,
  Read, WebSearch and WebFetch, a Python with piqtree (IQ-TREE 3), cogent3 and
  numpy (the server's own locked versions, but not the server), and MAFFT on
  `PATH`. Its prompt says these are installed.
- **Claude + phylokit-mcp**: this server (0.6.0 from PyPI) plus Read, with
  MAFFT on `PATH`. It needs Read because the server takes sequences as a
  string, not a file path.

The comparison that matters is the third against the fourth: the same engine
(IQ-TREE 3 through piqtree, MAFFT) driven by Claude from a shell, or through
the server.

Run 2026-10-16 (Claude alone, Python, phylokit-mcp) and 2026-10-16/17 (Read
only), Claude Code 2.1.302, model `claude-opus-5-5`.

![Right, wrong or declined per setup](agent_qa.png)

| Setup | Right | Wrong | Declined | Invalid | Questions right in all 3 runs | Cost / answer, median (mean) | Seconds / answer, median (mean) | Slowest |
|---|---|---|---|---|---|---|---|---|
| Claude alone | 0 | 0 | 114 | 0 | 0 of 38 | $0.012 ($0.012) | 5.4 (5.8) | 11 s |
| Claude + Read only (control) | 76 | 24 | 12 | 2 | 21 of 38 | $0.18 ($0.27) | 34.5 (79.2) | 1688 s |
| Claude + Python, MAFFT, web | 113 | 0 | 0 | 1 | 37 of 38 | $0.029 ($0.030) | 13.3 (15.6) | 43 s |
| Claude + phylokit-mcp | 110 | 4 | 0 | 0 | 36 of 38 | $0.55 ($0.55) | 112.2 (133.7) | 663 s |

Out of 114 answers per setup. The three runs of a question are not
independent, so the questions column is the fairer count. An invalid row
counts as not right, whatever its reply scored. The numbers come from
[`summary.tsv`](summary.tsv) and [`paired.tsv`](paired.tsv), which
`analyze.py` writes from [`results.jsonl`](results.jsonl) and
[`results_read.jsonl`](results_read.jsonl); 95% intervals are bootstrap over
questions.

## What it shows

- **The server does not make Claude more accurate than a Python shell with
  the same libraries.** 36 of 38 questions right in every run with the server,
  37 of 38 with Python. They disagree on 3 questions, 2 for Python and 1 for
  the server (exact sign test p = 1). Accuracy over answers: 0.965 (95% CI
  0.90-1.00) with the server, 0.991 (0.97-1.00) with Python.
- **The server is slower and dearer.** Per question it took 6.0x as long (95%
  CI 4.4-8.3) and cost 14.4x as much (9.9-20.2) as the Python shell.
- **The server's 4 wrong answers are all homology, and all counting slips
  after the alignment.** `align_sequences` returns the alignment; the model
  then counted residues along the aligned rows by hand, and said so in each
  reply: 337 and 339 against a key of 338 (q13), and 356 against 366 (q14).
  The server has no tool that maps a residue to its column, so this step is
  left to the model.
- **Tools help over reading alone where the question needs computation.**
  Right in every run: 37 (Python) and 36 (server) of 38, against 21 for Read
  only; discordant 16:0 (p = 3e-5) and 16:1 (p = 3e-4).
- **Python's one invalid row** printed a traceback naming the alt
  environment's path, which then sat inside the benchmark directory (q10,
  run 2); the answer was right but is counted as not right.

![Accuracy beside the no-LLM baselines](agent_qa_accuracy.png)

![Right, wrong or declined by kind of question](agent_qa_by_category.png)

![Cost and time per answer, every run](agent_qa_cost_time.png)

## Limitation: half the kinds of question are solvable by reading

Version 2 was built so that the inputs would be past what can be worked by
eye. The Read-only control shows that held for three kinds of question and
not for the other three:

| Kind | Read only: right / 18 (24 for model) | Questions right in all 3 runs | How it got there (from its replies) |
|---|---|---|---|
| JC69 distance | 18 | 6 of 6 | counted differences in blocks of 10, applied the JC69 formula |
| RF distance | 17 | 5 of 6 | listed every split of both trees by hand |
| Homologous residue | 16 | 4 of 6 | matched the flanking motif in the other sequence, no alignment |
| Sister taxon | 11 (2 invalid) | 3 of 6 | |
| Model with lowest BIC | 11 | 3 of 8 | base composition plus a transition share counted on a sample |
| Parsimony-informative sites | 3 (12 declined) | 0 of 6 | |

On distance, RF and homology, a model with no tools gets almost every
answer, so those 18 questions cannot show a difference between a tool and no
tool. The no-LLM baselines point the same way: matching the residue's +-6
flanking motif (no alignment) gets 5 of 6 homology questions, because the
no-gap margin that makes the key certain also makes the residue findable by
its context. The distance tolerance (3 SDs of the estimate) admits a careful
hand count.

On the other three kinds (20 questions, a subset chosen **after** the Read
control ran, so post hoc):

| Setup | Right / 60 | 95% CI | Questions right in all 3 runs |
|---|---|---|---|
| Claude alone | 0 | 0.00-0.00 | 0 of 20 |
| Claude + Read only | 25 | 0.23-0.60 | 6 of 20 |
| Claude + Python, MAFFT, web | 60 | 1.00-1.00 | 20 of 20 |
| Claude + phylokit-mcp | 60 | 1.00-1.00 | 20 of 20 |

There the tools separate cleanly from reading (14:0 discordant questions
for each, p = 1.2e-4), and the server and Python still tie. On this subset the
server took 11.9x as long (95% CI 8.9-15.6) and cost 32x as much (27-39) as
Python.

The Read control's 2 invalid rows (q01 and q02, run 1) are a harness fault,
not the model's: the first launch put run directories under a path that
names the server and the benchmark, so every Read of the data file tripped
the leak rule. The launch was stopped after those two runs and restarted
with run directories in the system temp dir. By the one rule for every
setup the two rows stay in and count as not right, which costs Read only one
right answer (q02). Over valid rows alone it scores 0.679.

## Why version 2

Version 1 put 6-8-taxon data in the prompt. In a partial sweep (87 of 414
runs), Claude with no tools answered 28 of 30 right, including every sister
taxon, RF and homologous-residue question: the inputs were small enough to
work by eye, so the sweep said nothing about the server. Version 2 scales the
inputs (16-24 taxa and 800-1200 sites, 30-50-tip trees, ten indel-rich
sequences of 300-700 residues) and gives them as files only. A test checks
that no question carries its data in the prompt. The Read control above
shows the scaling was enough for three kinds and not for RF, distance and
homology.

## What it does not show

- **The data are simulated.** Trees are drawn at random and sequences evolved
  down them under JC69 or the HKY85 family (indels added for the alignment
  questions). Real alignments have rate variation, compositional shifts and
  alignment error that none of these questions contain.
- **The questions are answerable on purpose.** A draw was kept only when
  piqtree's ML tree contained the asked cherry (sister questions), MAFFT put
  the asked residue in the true column (homology), the generating model won
  BIC by at least 2 (model choice), and both distance estimates fell in the
  tolerance (distances). Cases where the data do not support the answer are
  not tested, and those are where the server's bootstrap support matters most.
- **Support values are not scored.** No question asks how well supported a
  clade is.
- **Input size is bounded by the server's interface.** The server takes the
  alignment as a string argument, so the model has to read the file and copy
  its contents into the tool call. The largest file here is 20 KB, and copying
  it costs output tokens. Much larger alignments would need a path argument.
- **One model, one date.**

## Questions

[`questions.json`](questions.json) holds the 38 questions with their answers
(version 2, seed 20261016, built 2026-10-17 UTC by `build_questions.py`).
The data each question needs is written as a file (`alignment.fasta`,
`sequences.fasta` or `trees.nwk`) into the run's empty working directory, and
the prompt names the file.

| Kind | Count | Server tool it exercises | Answer key |
|---|---|---|---|
| Sister taxon in the ML tree (24 taxa, 800 sites, JC69, internal edges 0.01-0.05) | 6 | `infer_tree` | the asked taxon's cherry partner in the generating tree |
| Robinson-Foulds distance (30-50 tips, 3-30 NNI moves apart) | 6 | `compare_trees` | splits of the two Newick strings as written, symmetric difference |
| Homologous residue after alignment (10 seqs, 3 DNA + 3 protein, indel-rich) | 6 | `align_sequences` | the true alignment from an indel simulator that tracks every residue's column |
| Model with the lowest BIC among JC69, K80, F81, HKY85 (16 taxa, 1200 sites) | 8 | `select_substitution_model` | the generating model, 2 per model |
| JC69 distance between two taxa (16 taxa, 1000 sites) | 6 | `infer_tree` (branch lengths) | path length between the two tips in the generating tree |
| Parsimony-informative sites (16-24 taxa, 800-1000 sites) | 6 | `infer_tree`, `align_sequences` (`n_parsimony_informative`) | counted from the alignment |

The key comes from [`truth.py`](truth.py): numpy and plain Python, no
`phylokit_mcp`, cogent3 or piqtree, so a bug in the server or the libraries
under it cannot also be in the key (`tests/test_agent_qa_score.py` checks the
imports). piqtree and MAFFT are used only to check that a draw is answerable.
`questions.json` records the draws per category. 2 of 8 homology draws were
discarded because MAFFT disagreed with the true column. No sister, model or
distance draw was discarded. In 2 of the 6 sister questions the taxon nearest
the asked one by JC distance is not its sister, so the answer needs a tree,
not a distance lookup. Rebuilding gives the same questions and answers.

## Scoring

Each reply ends with `ANSWER: <value>` or `ANSWER: UNKNOWN`. `score.py`
decides the verdict from that line alone:

- **right**: the value matches the key.
- **wrong**: the value does not match.
- **declined**: `UNKNOWN`.
- **no answer**: there is no `ANSWER:` line.

Matching rules by kind, and why each is as strict as it is:

- **Sister taxon**: the taxon name, without case. The ML tree on the data was
  checked to contain the cherry.
- **RF distance, informative sites**: the integer, exactly. Both are
  deterministic functions of the file. "0.6" is not read as 6.
- **Homologous residue**: the integer, exactly. Neither sequence has a gap
  within 6 columns of the asked one, and MAFFT agrees with the true alignment.
- **Model**: any name of the right model. IQ-TREE writes HKY85 as `HKY+F` and
  F81 as `F81+F`; `K2P+F` is HKY85, not K80.
- **Distance**: within 3 large-sample SDs of the JC69 distance estimate at
  the true distance, sd = sqrt(p(1-p)/n) / (1 - 4p/3). An estimate from data
  differs from the generating value by sampling error, so an exact match would
  score correct work as wrong. The best single guessed number would be right
  on 2 of 6.

A wrong answer that is still a legitimate option (a taxon in the alignment,
one of the four models) is also flagged `wrong_id`. Each row also records,
read from the transcript, whether any tool call named `phylokit` or the
benchmark's key file (`mentions_phylokit`, `mentions_key`).

One validity rule for every setup (`run.py`'s `validity`): the run exited 0
with a result, its tools were exactly the setup's, a setup with tools made at
least one call, and no tool input or output touched the key, the benchmark's
code, its transcripts or another run's directory. Invalid rows are kept and
reported, never retried.

`controls.py` holds the positive control (the key scores right, known-wrong
answers score wrong, checked before any run) and the no-LLM baselines: the
best constant per kind, chosen with the key (an upper bound), and one
shortcut per kind that a reader without tools could take. They are in
[`baselines.tsv`](baselines.tsv). Tests: `tests/test_agent_qa_score.py`,
`tests/test_agent_qa_controls.py`.

## Reproduce

```bash
export MAFFT=/path/to/mafft                                 # if MAFFT is not on PATH
uv run python benchmarks/agent_qa/build_questions.py        # rebuild the key (writes questions.json)
uv run python benchmarks/agent_qa/run.py --setup-alt        # once: the alt arm's Python, without the server
uv run python benchmarks/agent_qa/run.py --reps 3 --workers 2   # 342 claude -p runs; resumable
uv run python benchmarks/agent_qa/run.py --arms read --reps 3 --workers 1 \
    --results benchmarks/agent_qa/results_read.jsonl        # the Read-only control, 114 runs
uv run python benchmarks/agent_qa/run.py --rescore          # re-score stored replies (add --results for the read file)
R=benchmarks/agent_qa/results.jsonl,benchmarks/agent_qa/results_read.jsonl
uv run python benchmarks/agent_qa/analyze.py $R             # summary.tsv, paired.tsv, baselines.tsv
Rscript benchmarks/agent_qa/plot.R $R                       # figures
```

Each run is an isolated `claude -p`:
- `--setting-sources ""` and an empty working directory in the system temp
  dir, holding only the question's data file.
- `--strict-mcp-config`, so no other MCP servers.
- The same short system prompt in every setup. The user prompt is the same
  apart from the alt arm's one-sentence note on what its shell has.

Raw transcripts go to `~/phylokit-mcp-artifacts/agent_qa/runs/` (or
`$PHYLOQA_ARTIFACTS`), outside the repo.
