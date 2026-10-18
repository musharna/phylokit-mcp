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

Compare trees with `compare_trees`, and numbers with a tolerance — never by
string equality.

Threads are pinned to 1 before piqtree is imported: likelihood sums accumulate in
thread-completion order, floating-point addition is not associative, and
near-tied topologies can flip on the last bits. Pinning is necessary, not
sufficient.
