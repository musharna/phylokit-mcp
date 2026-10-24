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

### Gate change and resume (2026-10-21, committed 22:45:35 EDT, coordinator's decision, before any more rows)

The user's /usage reads 66% weekly, so the Sat 23:00 window is the binding
limit and the 0.98 readings (Sun 12:00 window) are judged spurious. From
here on, a transcript's seven_day >= 0.90 triggers a fresh probe
(`claude -p "say ok" --model claude-haiku-5-5 --output-format stream-json --verbose`,
the newest event's `unifiedWindows.seven_day`; `run.probe_seven_day`). The
sweep stops only if the probe also reads >= 0.90. A probe that gives no
reading stops it. Otherwise both readings are logged to `versions.json`
(`usage_gate`) and the sweep continues. The five-hour pause, the narrowed leak rule
(primary) and Claude Code 2.1.307 are unchanged. Test:
`test_a_high_seven_day_reading_stops_only_if_a_fresh_probe_agrees` (0.98
reading + 0.66 probe continues; 0.98 + 0.95 stops). It fails on the old gate
(mutate-run, 3/3 mutants killed). A real probe just before the commit read 0.67. Resume
from bare q06 rep 2.

### Load and outage (2026-10-22, 04:50 EDT, coordinator's decision for load; arm-neutral)

- **Load.** From the next resume on, every row records the machine's
  1-minute load average at the run's start and end (`load1_start`,
  `load1_end`), timeouts included. The once-only re-runs at the end of the
  sweep start only when the 1-minute load is below 22 (the core count). The
  load is checked every 10 minutes for up to 2 hours, then the re-runs start
  anyway. Why: mcp q01 rep 3 hit the 1,800 s timeout at a load of 32 to 37
  from other sessions' jobs. The server's tree search runs on this machine,
  so contention falls on the mcp arm. The 234 rows before this point carry no
  load reading.
- **DNS outage.** From about 04:15 EDT the machine could not resolve
  `api.anthropic.com` (`EAI_AGAIN`). bare, alt and mcp q02 rep 3 ended in
  `api_error` after about 450 s of retries. The sweep was stopped at 04:47
  EDT (234 rows, all complete; the run in flight, bare q03 rep 3, left no
  row). The cause was the Windows-side VPN's tunnel mode breaking WSL's DNS
  proxy. The sweep resumes on the coordinator's message, not when the name
  resolves again. The 3 outage rows (bare, alt and mcp q02 rep 3; each
  transcript ends `API Error: Can't reach the API server ... (EAI_AGAIN)`)
  are outage failures, NOT the registered single re-run. At 05:15 EDT they
  were moved to `results.outage_2026-10-22.jsonl` (marked `outage`). Their
  transcripts were copied to
  `~/phylokit-mcp-artifacts/agent_qa/model-tiers-2026-10-20/outage-2026-10-22/`,
  and the resumed sweep runs those triples fresh. The registered rule
  (`--drop-api-failures`, then one re-run) still applies to them afterwards,
  as to every row. The timed-out row (mcp q01 rep 3, load 32 to 37) is not an
  outage row; it gets the registered single re-run at the end. 231 rows
  remain in `results.jsonl`.

### Second DNS outage and DNS precheck (2026-10-22, recorded 15:00 EDT, coordinator's decision; arm-neutral, before any analysis)

- **Second outage.** The sweep resumed at 05:35 EDT and wrote 9 rows (240 in
  all). At 06:4x it entered a registered five-hour pause (five_hour 0.88,
  seven_day 0.87; until 08:20 EDT). DNS failed again at about 06:48 and the
  coordinator stopped the sweep in that pause at 06:58 (SIGTERM; no run in
  flight, no half-written row; the pause entry is in the sweep log, not in
  `versions.json`). The timed-out mcp q04 rep 3 row (06:15:15 to 06:45:38
  EDT, load 17 to 18) is NOT an outage row: its window ends before 06:48, its
  transcript has no network error, and its last API call answered at 06:45:34.
  It timed out because Haiku asked `infer_tree` for 100 bootstrap replicates,
  which ran past Claude Code's 1,815 s silent-tool limit. It stays a
  registered timeout and gets the single re-run at the end.
- **DNS precheck.** From the next resume on, `getent hosts api.anthropic.com`
  runs before each run (`run.py --dns-precheck`, inside the usage gate, after
  any pause) and before each pass of the resume script. On failure the sweep
  stops cleanly: no new run starts and finished rows stay. A run already in
  flight when DNS fails can still end in an outage error; such rows are set
  aside by time window as outage rows, as with q02 rep 3.
- **Third stop (2026-10-22, 15:34 EDT, the registered gate).** After resuming at 14:55, the sweep wrote 3 rows (243 in all). Then the gate stopped it at seven_day 0.90, confirmed by its fresh probe (0.90). The last row, mcp q05 rep 3, is a registered timeout: Haiku asked `infer_tree` for 100 replicates, past the 1,827 s silent-tool limit, and there was no network error. The sweep resumes on the coordinator's word after the weekly reset (Sat 10-24, 23:00 EDT).

## Amendment, 2026-10-23 ~19:20 EDT: API billing (before any further run)

- **Billing.** The user's decision: every remaining run in this design (both
  tiers) bills an Anthropic API key, not the Claude subscription. Same model
  IDs, prompts, arms, reps, harness and controls. Claude Code is 2.1.309 now
  (recorded per row).
- **How.** A benchmark-only `claude` wrapper, first on `PATH` for the sweep,
  loads the key and then `exec`s the real CLI (so the harness's timeout still
  kills Claude itself). It records each run's cost (Claude Code's
  `total_cost_usd`; for a killed run, estimated from its token counts) and
  refuses to start a run once the total reaches $171 of the user's $175 cap
  (exit 75, no transcript). A run whose `init` event does not show
  `apiKeySource: ANTHROPIC_API_KEY` stops every later launch.
- **Gate.** API calls report no subscription windows, so the seven-day stop
  and the five-hour pause no longer fire. `run.py` instead asks the spend gate
  (`BENCH_SPEND_GATE`) before every run and, on a refusal, stops cleanly as
  at the seven-day stop: finished rows stay, nothing new starts. Test added,
  seen to fail with the check disabled (`mutate-run`).
- **Set aside by exit code, not outcome.** A run the wrapper refuses (exit
  75) is a harness failure, not a result: it is set aside and re-run once.
- **Concurrency.** Up to 4 servers' sweeps run at once (2 before). The
  workers per sweep are unchanged.
- **Harness fault, 20:38-20:52 EDT, and set-aside (recorded before any analysis).**
  The first API-billed launch hid the key from the agent's tools with
  `CLAUDE_CODE_SUBPROCESS_ENV_SCRUB=1`. That setting also forces Claude
  Code's permission mode to `default` (checked live: with it, a run asking for
  `bypassPermissions` starts in `default`; without it, on either billing
  source, it starts in `bypassPermissions`). The data-aggregator arms ask for
  `bypassPermissions`, so their tools were denied (91 of 91 Haiku runs had
  denials; median 7 turns against 21.5 before). The coordinator stopped every
  sweep at 20:52. Rule, by time window and not by outcome: **every run
  started under that wrapper, on every server and tier, is set aside and
  re-run** (`*.scrub_2026-10-23*`, `_scrub/` for data-aggregator), including
  servers whose arms ask for `default` mode and showed no change. The key now
  reaches Claude Code through an `apiKeyHelper` flag setting (`--settings`),
  never the environment; the run's `init` must show `apiKeySource:
  apiKeyHelper` and the permission mode the harness asked for, or every later
  launch stops.
