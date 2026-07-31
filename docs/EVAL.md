# Eval — recovering a tree that was known before inference ran

Every other test in this repo checks that the code does what the code intends. A
regression test can tell you a number stopped changing; it cannot tell you the
tree is wrong. This file records the check that compares inference to a ground
truth fixed in advance.

## Design

Sequences are evolved down a **known** topology `((a,b),(c,d))` under
Jukes-Cantor, inference is run on the result, and the inferred tree is scored by
**Robinson-Foulds distance** to the generating tree. RF compares splits, so it
is invariant to the many Newick spellings of one topology.

**The simulator is written inside the test on purpose.** Using cogent3 to evolve
the sequences would share a library with the inference path; forty lines of JC69
share nothing, and the generating process stays fully visible:

    p(substitution on a branch of length t) = 3/4 · (1 − e^(−4t/3))

## Measured

Recovery rate over 10 seeds, before any of it was written into an assertion:

| internal edge | tip length | 1000 sites | 100 sites | 30 sites  |
| ------------- | ---------- | ---------- | --------- | --------- |
| 0.25          | 0.05       | 10/10      | 10/10     | **10/10** |
| 0.02          | 0.30       | 9/10       | 5/10      | 3/10      |

**The discriminating variable is the internal edge, not the alignment length.**
A 30-site alignment across a 0.25 internal edge recovers every time; a 1000-site
alignment across a 0.02 edge with long tips does not. That is this server's own
doctrine — _alignment length is not evidence_ — arrived at from the opposite
direction, and it is why the eval asserts on the edge regime rather than on
sequence length.

## The two controls

A benchmark that always scores full marks measures nothing.

**Signal destroyed.** Shuffling each site independently across taxa preserves
every base composition and the exact alignment dimensions while destroying
shared ancestry. Recovery must fail. Inference still returns a tree, confidently
— which is the whole point: an unsupported topology looks exactly like a
supported one from the outside.

**Resolution.** The hard regime must recover _less often_ than the easy one,
otherwise the eval could not notice inference getting worse.

## Seen to fail

| mutant                                                    | result        |
| --------------------------------------------------------- | ------------- |
| truth compared against the wrong topology `((a,c),(b,d))` | RED (7 tests) |
| simulator puts `a` and `c` in the same clade              | RED (6 tests) |
| RF replaced by a constant zero — everything "recovered"   | RED (2 tests) |

The third is the instructive one. Faking RF as zero leaves **five of the seven
tests passing**, because the recovery assertions are then trivially satisfied.
Only the shuffle control and the resolution test catch it. A recovery benchmark
without those two would have reported a perfect score for an inference engine
that had been replaced by a constant.

## Not covered

Branch LENGTHS are not scored, only topology. Support values are checked
elsewhere for calibration in the ordinary test-suite sense, not against a
simulated truth.
