# Reproducibility, stated precisely

Measured, not assumed:

- **Not bit-exact, in any setting.** The same request with the same seed — on
  repeat in one process, or in a fresh process — can return branch lengths and a
  log-likelihood that differ in the trailing digits (eight fresh processes gave
  five distinct log-likelihoods, spread ~2e-6). IQ-TREE reads the wall clock
  during its search: freezing `gettimeofday()` alone made every run
  bit-identical. piqtree exposes no option to take the clock out, so the server
  reports `deterministic_across_processes: false` rather than promise it.
- **Support can move by a replicate flipping.** Over six repeated 50-replicate
  calls, three of four clades were bit-identical and one moved **0.02**, well
  inside the bootstrap's own sampling error (~0.07 at 50 replicates). The
  topology and every conclusion were unchanged. The column resampling itself is
  numpy-seeded and exact.
- **UFBoot2 support moves further, and the seed does not pin it.** IQ-TREE
  draws the UFBoot2 replicates inside its clock-dependent search. On the 60-site
  test fixture at one seed, 48 fresh processes gave one clade list but put the
  weakest clade anywhere from 0.22 to 0.34; 12 repeats in one process ranged
  0.07. With the clock frozen every run was identical. The server reports 0.15
  for `support_method="ufboot"` and `bootstrap_resampling_exact: false`. A clade
  near the 0.95 line can be called supported on one call and not on the next.

Compare trees with `compare_trees`, and numbers with a tolerance — never by
string equality.

Threads are pinned to 1 before piqtree is imported: likelihood sums accumulate in
thread-completion order, floating-point addition is not associative, and
near-tied topologies can flip on the last bits. Pinning is necessary, not
sufficient.
