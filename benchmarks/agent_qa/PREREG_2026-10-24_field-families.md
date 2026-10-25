# Pre-registration amendment: field families, 2026-10-24

Committed before any model has seen these questions. Anything run that is not
listed here is exploratory and labelled so.

## Why

The 38-question set ties at every tier because each question prescribes its
method (the model, the tree builder), so it sits at ceiling. These three
families take their designs from how phylogenetics scores its own tools
against simulated truth (ModelFinder, UFBoot2, SATé), give the task in prose,
and never name a substitution model.

## Families (10 questions each, `field_2026-10-24/questions.json`)

Builder: `field_2026-10-24/build_field.py`, seed 20261024, each draw from its
own seed, kept in draw order. Every generating tree is drawn by `truth.py`
(numpy only); every sequence is simulated by AliSim (IQ-TREE 3.1.3). The
IQ-TREE CLI and MAFFT v7.526 are used only to check draws and to compute the
reference pipeline and the shortcuts. Nothing in the key comes from
phylokit-mcp (`tests/test_agent_qa_field.py` reads the imports from the AST).
AliSim outputs, MAFFT alignments and generating trees are written to
`~/phylokit-mcp-artifacts/agent_qa/field-families-2026-10-24/key/<qid>/`, and
their sha256 is in the key.

**P1 `lba_sister`, long-branch attraction.**
- 8-20 taxa, internal edges 0.01-0.04, tips 0.02-0.1, two long tips
  (1.0-2.0) in two different cherries, 1000-3000 sites, GTR+F+G4 with
  alpha 0.3-0.5.
- Ask: which taxon forms a cherry with the long-branch taxon.
- Kept only if the ML tree under the generating model (GTR+F+G4) AND the ML
  tree under ModelFinder's BIC pick both have the true cherry, and the JC ML
  tree does not.
- Score: exact.
- Retention: 18 drawn, 10 kept. 8 discarded: the generating-model ML missed
  the cherry on 6, and the JC ML had it (no attraction) on 2.

**P2 `supported_clades`, reportable clades at 95%.**
- 30-100 taxa, 300-1000 sites, GTR+F+G4 with alpha 0.5-1.5. Internal edges
  are 0.02-0.12, except that each is set near zero (0.001-0.01) with
  probability 0.5.
- Ask: every clade (split) you would report as well supported at the 95%
  level. The answer is one Newick tree holding exactly the reported clades.
- Score: precision and recall of the reported splits against the generating
  tree's. Right iff both are at or above the band.
- Band: the reference pipeline's precision and recall, minus 3 binomial SE.
  The reference pipeline is IQ-TREE ModelFinder (BIC) + UFBoot2 with 1000
  replicates, seed 1, and it reports a clade when support ≥ 95. The SE uses
  the Laplace estimate (x+1)/(n+2), so a reference at 100% precision still has
  a band; n is the number reported (precision) or n − 3 (recall). The 95
  threshold is UFBoot2's own, not the server's 0.70.
- Discards: a band that admits the star tree, or a draw where the
  unbootstrapped ML tree reporting every clade lands inside the band. The
  second is a shortcut filter. It is applied in a named regime (near-zero
  edges) and its count is published.
- Retention: 14 drawn, 10 kept. 4 discarded, all because the every-clade ML
  tree was inside the band.
- Reference precision ranges 0.83-1.00 and recall 0.42-0.63. The bands are
  precision ≥ 0.67-0.88 and recall ≥ 0.21-0.45.

**P3 `align_tree_homology`, unaligned indel-rich DNA.**
- AliSim `--indel 0.05,0.1` (default size distributions) on 20-50 taxa, root
  length 500-1000, internal edges 0.02-0.1, tips 0.05-0.2, GTR+F+G4 with
  alpha 0.5-1.5. Only the unaligned FASTA is given.
- Answer `<residue> | <Newick>`. Right iff both parts are right.
- (a) Tree: right iff its nRF to the generating tree is ≤ the reference's mean
  nRF + max(3 sd, one split), where one split is 2/(2n − 6) in nRF for n taxa.
  The mean and sd are taken over 20 seeded runs: MAFFT --auto, then IQ-TREE ModelFinder
  + ML at seed 1, then ModelFinder's pick fitted at seeds 2-20. The first
  build ran ModelFinder at all 20 seeds. It was stopped after one batch,
  because that design projected the build past 2 h; no question from it was
  kept.
- (b) Homology: exact against AliSim's true alignment. The existing homology
  rule applies: the column is in the second half, there is no gap within 6
  columns in either row, the two positions differ by ≥ 5, and MAFFT --auto
  agrees with the true alignment. A draw that fails any of these is
  discarded.
- Retention: 11 drawn, 10 kept. 1 discarded: MAFFT disagreed on the column.
- **Tolerance amendment (2026-10-24, coordinator decision, before any arm
  ran).** The rule was mean + 3 sd. It is now mean + max(3 sd, one split).
  - Why: the 3 sd term was meant to absorb tree-search noise. On every
    question the 20 reference runs give the same nRF (sd 0), so 3 sd alone
    meant "match the reference pipeline exactly". That fails a legitimate
    alternative alignment, such as MAFFT without --auto, by a single split.
  - Now the floor binds on all 10. A tree may have one more wrong split than
    the reference pipeline's.
  - The reference has 0 wrong splits on 4 questions, 1 on 5 and 2 on 1, so
    1-3 wrong splits are allowed.
  - Tolerances run from nRF 0.030 to 0.074.
  - The no-alignment shortcut (nRF 0.81-1.0) is still wrong on 10 of 10.
  - `tests/test_agent_qa_field.py::test_p3_admits_one_split_beyond_the_reference_and_not_two`
    checks this on every question. It was seen to fail on the 3 sd key.

Inputs stay inside the server's caps (≤ 200 taxa, ≤ 100,000 sites). No server
feature was added. Gaps that may cost the server, noted and not fixed:
- `infer_tree` has a standard bootstrap only (20-1000 replicates), and a
  standard bootstrap at 95 is more conservative than UFBoot2 at 95, so P2
  recall may fall under the band.
- `simulate_alignment` has no indel model; it is not needed by any family.

## Controls (`tests/test_agent_qa_field.py`)

- **Oracle:** the 30 key answers score right, and a wrong answer to each
  scores wrong. `run.py` checks this before any run.
- **Positive control:** the reference pipeline re-run at seed 2, never the
  server, is right on 30 of 30. Further checks:
  - with IQ-TREE present: P1 at seed 3;
  - the alt arm's own piqtree with GTR+F+G4 on P1.
- **Shortcuts:**
  - P1, the JC ML tree: wrong on 10 of 10. Each time it pairs the asked long
    branch with the other one. The no-+G GTR tree is wrong on 8 of 10.
  - P2, the ML tree with every clade reported: wrong on 10 of 10. That count
    is after the filter above; 4 draws were discarded for it.
  - P3, no alignment (sequences padded at the end, then ML): wrong on 10 of
    10, with nRF 0.81-1.0. Same relative position plus the reference tree:
    wrong on 9 of 10.
- **Seen to fail.** 8 of 8 mutants were killed under `mutate-run`:
  - the scorer:
    - P2 recall only;
    - P3 ignores the residue;
    - P3 ignores the tree;
    - P3 strict band edge;
  - the builder's SE without the Laplace adjustment;
  - the positive control reading a shortcut's reply;
  - a web tool left in the web-off alt;
  - the key: p1q01's truth set to the JC answer.
- **After the tolerance amendment**, 8 of 8 more mutants were killed:
  - the builder:
    - no floor;
    - min for max;
    - one split as 1/(2n−6);
    - one split as 4/(2n−6);
  - the scorer:
    - a strict edge;
    - no float slack;
    - the tree ignored;
  - the key: p3q02's tolerance set back to 0.
- The controls were re-run on the amended key:
  - the reference at seed 2 is right on 30 of 30;
  - no alignment is wrong on 10 of 10;
  - the same relative position is wrong on 9 of 10;
  - P1 and P2 are unchanged.
- **Caught in the build.** The first key rounded each band to 4 decimals.
  That put the P3 edge below the reference's own nRF, and the seed-2
  reference failed 6 of 10. The bands are now rebuilt from integer counts
  (`exact_tolerances`), unrounded.

## Pilot

- Arm `alt` only, model `claude-sonnet-5-5`, 1 rep, all 30 questions, one
  worker.
- **Web off.** The key's `alt_web: false` gives the alt arm `Bash,Read` only;
  the validity rule requires exactly those tools, and a run whose init event
  lists WebSearch or WebFetch is invalid.
- The alt note naming the installed libraries (piqtree, cogent3, numpy, MAFFT)
  is unchanged (design-note decision). It is the note every earlier sweep
  used, so a difference comes from the questions, not from the prompt.

```bash
export PHYLOQA_ARTIFACTS=~/phylokit-mcp-artifacts/agent_qa/field-families-2026-10-24
Q=benchmarks/agent_qa/field_2026-10-24/questions.json
uv run python benchmarks/agent_qa/run.py --questions $Q --setup-alt   # once, no model
uv run python benchmarks/agent_qa/run.py --questions $Q --arms alt --model claude-sonnet-5-5 \
  --reps 1 --workers 1 --results benchmarks/agent_qa/field_2026-10-24/pilot_alt_sonnet.jsonl
```

`--dry-run` on the same command prints the alt tools and every prompt, and
runs nothing.

## Headroom rule

- A family has **headroom** if alt is not right on ≥ 3 of its 10 questions.
  Wrong, declined, no answer and invalid all count as not right.
- Only headroom families get an mcp arm (`phylokit-mcp==0.7.0`, Sonnet 5.5,
  1 rep), within what is left of the shared $40 pilot budget.
- A family without headroom is reported as at ceiling and gets no further
  runs.
- Any full run is pre-registered separately. A tie on a full run is an
  acceptable result and will be reported as one.

## Amendment, 2026-10-25: P2 and P3 mcp arm re-run on a revised server

The mcp arm on `phylokit-mcp==0.7.0` lost both headroom families: on P2 it got 0
of 10 right (alt 4), and on P3 4 right, 0 wrong and 6 timeouts (alt 5). The
transcripts trace that to three server seams, not to load. P2 failed the same way
10 of 10 at a load of 13-30.

- **No fast bootstrap, and silence.** `infer_tree` had the nonparametric
  bootstrap only, one full ML search per replicate. One replicate takes 25-93 s
  on these alignments, and the call sent nothing until it finished, so Claude
  Code aborted it at 1822 s.
- **Inline-only alignment.** `align_sequences` returned the alignment inline only.
  The arm has no file-writing tool, so the agent retyped ~15 x 962 characters into
  `infer_tree` and mis-copied a row each time.

The user decided (2026-10-25, "go ahead with 1-3") to fix those seams and re-run
the mcp arm on P2 and P3.

**This run is exploratory, and is labelled so in every table.** The server was
changed after the pilot showed what these 20 questions need, so a better score
here is not a held-out result. A confirmatory claim needs fresh draws from the
same builder: new seeds, the same discard rules, built before any arm runs.

- **Server:** `git+https://github.com/musharna/phylokit-mcp@dfd248cd5eab98db9b581cea3d9b20102553c89f`
  (branch `fix/field-seams`, draft PR #80, not released). It adds three things:
  - `infer_tree(support_method="ufboot")`, IQ-TREE's UFBoot2 with clades
    supported at 0.95;
  - a progress heartbeat every 15 s on every tool;
  - `align_sequences(out_path=...)`.
  The default support method is unchanged. Nothing tells the agent which method
  to pick, beyond the tool's own description and server instructions.
- **Fixed, unchanged from the pilot:** questions, key, bands, scorer, prompt,
  arm tools (Read + phylokit), model `claude-sonnet-5-5`, 1 rep, 2 workers, and
  the 1800 s run cap.
- **Reported:** right / wrong / timeout per question, beside the 0.7.0 mcp arm and
  the alt arm, plus which `support_method` and `replicates` each run used. A
  tie or a loss is reported as one.
- **Results** go to `field_2026-10-24/pilot_mcp_seams/results.jsonl` (their own
  dir, so `versions.json` pins one server).
- **Budget:** what is left of the pilot's $40, through the same `LAUNCH_STOP` gate
  (160.23). The run is estimated at about $9.
