# Confirmatory P2/P3 run, 2026-10-25/26

Pre-registered in `PREREG_2026-10-24_field-families.md`, "Amendment, 2026-10-25:
confirmatory P2 and P3 run on held-out draws" (9199023, 444e979; pushed before the
draws were built). Key and controls: b6cdb0b, built from draws 2000+, with 10 of
10 controls passing and 4 of 4 planted faults killed before any arm ran.

**Result: WIN** for `phylokit-mcp@dfcaed8` (draft PR #80) over the alt arm.

| | mcp right | alt right | alt timeouts (1800 s) |
|---|---|---|---|
| P2 `supported_clades` | 10/10 | 2/10 | 8 (7 no answer, 1 declined) |
| P3 `align_tree_homology` | 10/10 | 9/10 | 0 (1 wrong) |
| pooled | 20/20 | 11/20 | 8 |

- **Primary rule** (exact two-sided sign test, `analyze.sign_test`):
  - **as written:** 10 mcp-only to 0 alt-only questions, p = 0.002. Invalid
    counts as not right, and alt p3q07 is invalid on a false-positive leak flag
    (below).
  - **with that flag corrected:** 9 to 0, p = 0.004.
  - Either way it is a WIN.
- **mcp arm behaviour:**
  - all 20 runs called `infer_tree(support_method="ufboot")`;
  - in 19 of 20, the answer tree appears verbatim in a tool output (copied from
    `newick_*`);
  - 0 timeouts; runs took 71-232 s.
- **Cost:** $1.32 by the run's own records. The spend ledger moved
  $135.34 → $137.45, which includes estimated costs for the timeouts.

## What the win rests on (reported, not deciding)

- **Alt's timeouts.** All 8 of alt's P2 misses are 1800 s timeouts, with 9-16
  tool calls each. In all 8 the agent ran UFBoot itself, through piqtree
  (`bootstrap_replicates=1000` or `-B 1000`), on 30-90 taxa. It relaunched
  with changed options and polled background jobs, and did not finish. That is
  the same method the server runs in one call. Drop the 8
  alt-timeout questions and the remaining 12 split 1:0 (p3q02), p = 1, a tie.
  So on P2 the server wins on getting a supported tree inside the run budget,
  not on reaching a different answer by another route.
- **Load.** The laptop was shared with other sessions' jobs: 1-minute load
  20-77 on 22 cores during the run. Alt's timeouts started at loads of 24-67.
  Both arms compute on this machine (the server runs IQ-TREE locally too), and
  the mcp arm finished every run under the same load. Timing under load was
  not part of the pre-registered rule.
- **P2 by construction (stated in the amendment).** The P2 key is UFBoot2 1000
  replicates at 0.95, which is the server's `ufboot` method. P3's key (MAFFT +
  ModelFinder + ML on the IQ-TREE CLI) has no such overlap: there mcp is 10/10
  and alt 9/10.

## Harness defect found (follow-up, scoring unchanged)

The leak scan flags alt p2q03 and p3q07. `RUN_DIR` (`phyloqa-[a-z]+-[\w]+`)
stops at the first `-`, and Claude Code writes background-task output under the
cwd with `_` turned into `-` (`-tmp-phyloqa-alt-6alf-1m7-cwd`). So a run's own
directory reads as a foreign one (`phyloqa-alt-6alf`). No run read another run's
files or the key (`mentions_key` is false on all 40). The pre-registered rule
was applied as written.

## Files

- `results.jsonl`: 40 rows.
- `versions.json`.
- Transcripts:
  `~/phylokit-mcp-artifacts/agent_qa/field-families-2026-10-25-confirm/runs-claude-sonnet-5-5/`.
