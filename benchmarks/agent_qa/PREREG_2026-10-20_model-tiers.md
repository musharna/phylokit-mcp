# Pre-registration: model tiers, 2026-10-20

Committed and pushed before either confirmatory sweep. Anything run after a
sweep that is not listed here is exploratory and is labelled so.

## Question

On Claude Opus 5.5 the server tied the Python shell on accuracy (10-16: 36
against 37 of 38 questions right in every run; 10-18: 36 against 35). Opus
writes the piqtree and MAFFT scripts itself, so accuracy sat at a ceiling.
This asks whether phylokit-mcp lifts a weaker model's accuracy.

## Tiers

Claude Haiku 5.5 (`claude-haiku-5-5`), then Claude Sonnet 5.5
(`claude-sonnet-5-5`), same design. Both run in this session, Sonnet right
after Haiku (brief amendment, 2026-10-20 18:50 EDT). No other tier is
confirmatory.

## Server and setup

- Server: phylokit-mcp 0.7.0 from PyPI, `--server-spec phylokit-mcp==0.7.0`,
  its libraries constrained to the repo's `uv.lock` export (as on 10-18,
  `--constraints alt-env.requirements.txt`, identical to the 10-18 file).
  `check_mcp` reads the installed version back and stops on anything but
  0.7.0. 0.7.0 carries both post-sweep fixes from 10-18 (CHANGELOG 0.7.0,
  Fixed): tool calls run one at a time under a process lock (concurrent calls
  vs piqtree's chdir), and `select_substitution_model` reads IQ-TREE's raw
  output when the winner is a name piqtree cannot parse (bare `K3Pu`).
- Alt: the same locked piqtree 0.8.3, cogent3 2026.7.30a0, numpy 2.4.6, and
  MAFFT v7.526, checked against `questions.json`'s `checked_with` before any
  run.
- `run.py --model <id>`; the init event's `model` is stored in every row and
  in `versions.json` (`model_served`), and a run served by another model
  stops the sweep.
- Claude Code 2.1.305 (10-18 ran 2.1.302).

## Key

`questions.json` version 2 (38 questions, seed 20261016), unchanged. The key
is computed from simulators, with no live data, so nothing it depends on can
have moved; its answerability checks depend on piqtree 0.8.3 and MAFFT
v7.526, which both arms run. The oracle check runs before each sweep and
must pass (38 right, 38 wrong).

## Arms, reps, order

`bare`, `alt` and `mcp`, exactly as `arm_setup` defines them (the `read`
control is not re-run). 3 reps. Order: rep, then question, then arm (bare,
alt, mcp), so arms are interleaved by question. One worker. No Opus row is
reused in the primary test.

```bash
export PHYLOQA_ARTIFACTS=~/phylokit-mcp-artifacts/agent_qa/model-tiers-2026-10-20
D=benchmarks/agent_qa/tiers_2026-10-20
for M in haiku sonnet; do
  flock ~/.claude/jobs/9cc8114d/tmp/sweep.lock env HEAVY_RUN_MEM=6G heavy-run \
    uv run python benchmarks/agent_qa/run.py --model claude-$M-5-5 \
    --arms bare,alt,mcp --reps 3 --workers 1 --server-spec phylokit-mcp==0.7.0 \
    --results $D/$M/results.jsonl
done
```

Transcripts: `$PHYLOQA_ARTIFACTS/runs-claude-<tier>-5-5/`.

## Primary endpoint and rule

Per tier: questions right in all 3 runs (an invalid row is not right), mcp
against alt. Test: McNemar's exact test on the discordant questions
(`analyze.sign_test`, as on 10-16 and 10-18); interval: bootstrap 95% CI of
the difference mcp − alt in questions right in all 3 runs (10,000 resamples
of questions, seed 20261016). Holm across the 2 tiers.

- **LIFT** at a tier: Holm-adjusted p < 0.05, mcp right-only questions
  outnumber alt's, and the CI lies wholly above 0. **DROP**: the mirror image.
- **WIN**: a lift at one tier or more and no drop at the other. **LOSS**: a
  drop and no lift. **MIXED**: both. **TIE**: neither.
- Implemented as `analyze.tier_rule` (tested in
  `tests/test_agent_qa_controls.py`).

Secondary, reported, not deciding:

- **Floor:** the bare arm's right answers per tier (expected 0).
- **Substitute claim:** the tier's mcp arm against the 10-18 Opus alt rows
  (`rerun_2026-10-18/results.jsonl`), same statistics
  (`analyze.py --substitute=`). Not concurrent (2 days apart, Claude Code
  2.1.302 vs 2.1.305); labelled as such wherever reported.
- **Mechanism:** server tool calls per mcp run, and mcp runs with no server
  call ("the mechanism did not fire"); `diagnose.py`'s path/inline,
  `map_positions` and `models` counts.
- Cost and seconds per answer (the 10-18 ratios), and accuracy by kind.

## Controls

The oracle check, the no-LLM baselines (unchanged from 10-16), one validity
rule for every arm, the leak scan, and invalid rows kept and counted as not
right. Rate limits: `run.py` reads each transcript's last `rate_limit_event`;
a rejected request or five-hour use ≥ 0.85 holds later runs until the window
resets, and seven-day use ≥ 0.90 stops the sweep (rows kept, reported). After
each sweep, `run.py --drop-api-failures` moves every row the API ended (a
rejected request, an API error, or a timeout with no result; read from the
transcript only) to `results.api_failures.jsonl`, and the sweep is resumed to
run those triples once. Rows the API ends again stay, and count as invalid.

No other retries. An interrupted sweep is resumed, never restarted, and is
not stopped early on interim results.

## Analysis

```bash
for M in haiku sonnet; do
  uv run python benchmarks/agent_qa/analyze.py $D/$M/results.jsonl \
    --substitute=benchmarks/agent_qa/rerun_2026-10-18/results.jsonl
  uv run python benchmarks/agent_qa/diagnose.py $PHYLOQA_ARTIFACTS/runs-claude-$M-5-5
done
```

then the verdict, `analyze.py --tiers=haiku=$D/haiku/results.jsonl,sonnet=$D/sonnet/results.jsonl`
(`tier_rule` on full-precision paired rows; it prints the Holm p it decided on).

## Reporting

Every result is reported, whatever it is, in a dated README section. The
Opus results stay published.

## Pilot

No run has been made on either tier before this file. The sweep lock is
shared with five other benchmarks, so the pilot runs in the same hold of the
lock as the sweeps, immediately before Haiku
(`~/phylokit-mcp-artifacts/agent_qa/model-tiers-2026-10-20/scripts/sweep.sh`):
q13 and q23, all three arms, 1 rep, on each tier (12 runs), written to
`model-tiers-2026-10-20/pilot/`, never analysed. It stops the sweeps only on a
harness failure (non-zero exit, no result, wrong tool set or server not
connected, or a run served by another model), never on a score. Before
running, the oracle (38/38 right, 38/38 wrong), the alt env, MAFFT v7.526 and
`phylokit-mcp==0.7.0` (installed 0.7.0, piqtree 0.8.3) were checked.

The usage reading before this file: seven_day 0.69, five_hour 0.13 (18:35
EDT, newest transcript on the account).

## Deviations, recorded during the sweep (2026-10-21, 00:18 EDT, before any analysis)

- **Usage stop.** The Haiku sweep started at 23:03 EDT and the usage gate
  stopped it at 00:17 EDT, when a transcript read seven-day use 0.90 (five-hour
  0.73). 87 of 342 runs are in `tiers_2026-10-20/haiku/results.jsonl`: rep 1,
  q01 to q29, all three arms. They are kept. Nothing after the stop has run;
  Sonnet has not started. The account's seven-day window resets Sun 10-25,
  12:00 EDT. The sweep is to be resumed, not restarted.
- **Claude Code** auto-updated to 2.1.306 before the sweep (this file says
  2.1.305); `versions.json` records 2.1.306.
- **CPU load.** Other work shared the machine (load average about 15 of 22
  cores; from 23:50 a second benchmark lane). The arms are interleaved, so only
  seconds per run is affected.
- **A harness artifact in the validity rule** (as first recorded at 00:18;
  its description of the cause was wrong in one detail, corrected below).
  11 of the 29 alt rows were invalid, all of them right answers, and all for
  one reason: the leak scan matched `agent_qa` in a tool output that printed
  the harness's own tool install path, which sits under
  `~/phylokit-mcp-artifacts/agent_qa/model-tiers-2026-10-20/`. The mcp arm has
  no shell, so it cannot trip this, and an invalid row counts as not right.

### Decision: narrowed leak rule (2026-10-21, 00:22 EDT, coordinator's decision, before any more rows and before any analysis)

Every scan hit in the 87 rows (all arms) was listed and classified by hand.
There are 21 hits in 11 rows, all alt. None in bare or mcp, none in a tool
input, no other run's directory, no "named phylokit":

| Row | Hits | What printed the path |
|---|---|---|
| alt q16 | 1 | `which mafft`: `.../bin/mafft` |
| alt q03, q06, q21, q24, q27 | 1 each | the alt env's `site-packages` path (a `print` of a module's location, or `sys.path`) |
| alt q05, q20, q23, q25, q26, q27 | 2-6 each | Python tracebacks and `<module 'piqtree.iqtree' from '.../alt-env/lib/python3.14/site-packages/...'>` |

All 21 are paths inside the harness's own tool installs (`alt-env/`, `bin/`).
None touches the key, the benchmark's code, transcripts, results or another
run. There is no genuine leak among them.

Rule from here on, for ALL rows of ALL arms, as the PRIMARY result:
before the leak scan and the "named phylokit" check, any path inside the
alt env or the MAFFT link directory (absolute or `~/` form) is masked
(`run.mask_tool_paths`). A path that climbs out of an install (`/..`) is not
masked, and every other pattern applies unchanged: the key file, the
benchmark's code, `results.jsonl`, the artifacts and transcripts dirs, the
repo's worktrees, its GitHub directory and another run's directory still
invalidate (`tests/test_agent_qa_controls.py::test_the_narrowed_leak_rule_*`,
5/5 mutants killed). The artifacts are NOT moved, so the remaining rows run in
the same environment.

The rule as written is kept as a SECONDARY result: every row records
`valid_as_written`, `invalid_reasons_as_written` and `leaks_as_written`, and
`analyze.py --as-written` reports the result under it. `run.py --rescore`
applied both to the 87 rows: under the narrowed rule 0 are invalid (the 11 alt
rows become valid), and under the rule as written the same 11 are invalid as
before (the as-written column reproduces the original judgement on every row;
no verdict changed).

Why: the rule as written measured where the harness installed the alt arm's
tools, not leakage, and only an arm with a shell can print that path. The
narrowed rule counts more alt answers right, so it works against the server:
it is the conservative choice for the primary. Disclosure: when the rule was
chosen, these 11 rows' verdicts (all right) had been seen; no accuracy
comparison between arms had been run.

### Resume (2026-10-21, 20:20 EDT, coordinator's decision)

The Haiku sweep resumes from q30 rep 1 before the Sun 10-25 reset: the
seven-day window rolls, and a live probe at 18:12 EDT read seven_day 0.57 and
five_hour 0.49. The same command, environment and stop rule (seven_day ≥
0.90) apply, with the narrowed leak rule as primary. Claude Code is now
2.1.307 (the first 87 rows ran 2.1.306); the version is accepted, recorded in
every row (`claude_code`), and reported. Sonnet does not start.

### Second stop (2026-10-21, about 22:35 EDT, by the coordinator)

Stopped at 129 rows on a seven_day 0.98 reading from a second weekly window.
That reading came from another benchmark's transcripts; this harness's
transcripts showed only the other weekly window (0.62 to 0.65), so its gate
could not see it. The coordinator stopped the sweep with SIGTERM after run
42 of 255 (mcp q05 rep 2). All 129 rows are complete (43 per arm, rep 1
q01 to q38 and rep 2 q01 to q05). The run in flight (bare q06 rep 2) left no
row and no transcript, and it will run again on resume.
