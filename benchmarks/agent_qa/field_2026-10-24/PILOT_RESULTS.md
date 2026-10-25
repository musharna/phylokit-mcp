# Field-family pilot results (Sonnet 5.5, 1 rep)

Number right out of 10. Headroom means alt is not right on at least 3 of 10.

| Family | alt (Bash+Read, web off) | mcp 0.7.0 | mcp dfd248c, **exploratory** |
|---|---|---|---|
| P1 `lba_sister` | 10 (no headroom) | not run | not run |
| P2 `supported_clades` | 4 (4 timeouts) | 0 (10 timeouts at 1822 s) | 4 |
| P3 `align_tree_homology` | 5 | 4 (6 timeouts) | 6 |

- **Files:** alt is `pilot_alt_sonnet.jsonl`, 0.7.0 is `pilot_mcp/`, and dfd248c is `pilot_mcp_seams/`.
- **Why the dfd248c column is exploratory:** the server was revised after the 0.7.0 arm showed what these questions need (PREREG amendment, 2026-10-25). It is not a held-out result.
- **dfd248c run:**
  - Every run used `support_method="ufboot"`, and every P3 run wrote its alignment with `out_path`.
  - There were 0 timeouts; runs took 68-481 s and $1.87 in all.
- **dfd248c misses:** all 10 are trees the agent re-typed from the tool output and broke (`parsed: False`, unbalanced parentheses). Every answer that parsed was right.
- **Counterfactual re-score, not a run:** scoring the tree the server returned in each of those transcripts gives P2 10/10 and P3 10/10. P3 uses the agent's own residue and the server's ML topology. P2 uses the server's clades at ≥0.95, built with `splits_to_newick`.
  - The P2 key is UFBoot2 ≥95, the method the server now offers, so the P2 match is partly by construction.
- **What a confirmatory claim needs:** fresh draws from `build_field.py` (new seeds), run on the server version that returns `newick_supported`/`newick_topology` (PR #80, a820b74 or later), against alt.
