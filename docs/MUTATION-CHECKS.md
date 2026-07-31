# Mutation checks

A test that has never been seen to fail is not evidence. Each mutant below was
applied to the source, the suite was run, and the result recorded. The gate
asserts the suite is **green before mutating** — without that, an already-red
suite reports every mutant as killed.

Reproduce: `bash docs/mutate.sh` (restores every file afterwards and re-checks
green).

## Round 1 — 2026-07-30, at 0.1.0

| #   | mutant                                                    | file             | result                    |
| --- | --------------------------------------------------------- | ---------------- | ------------------------- |
| 1   | every clade's support forced to 1.0                       | `bootstrap.py`   | killed                    |
| 2   | splits not canonicalised — no reference-taxon anchor      | `splits.py`      | killed                    |
| 3   | bootstrap resamples **without** replacement (permutation) | `bootstrap.py`   | killed                    |
| 4   | trivial splits kept (terminal edges counted as clades)    | `splits.py`      | **survived, then killed** |
| 5   | parsimony-informative counts singleton sites              | `alignment.py`   | killed                    |
| 6   | ragged alignment accepted instead of refused              | `alignment.py`   | killed                    |
| 7   | AIC drops its parameter penalty                           | `inference.py`   | killed                    |
| 8   | model ranking reports the winner with no runners-up       | `inference.py`   | killed                    |
| 9   | the unsupported-clade advisory never fires                | `diagnostics.py` | killed                    |
| 10  | server claims bit-exact in-process reproducibility        | `engine.py`      | killed                    |

### Mutant 4 is the one worth reading

Removing the trivial-split guard changed **nothing** on the 7-taxon fixture the
whole suite was built on, because that tree shape never produces such a split.
The guard was untested and looked tested.

The discriminating input is a **ladder rooted at a tip**:
`(A,(B,(C,(D,E))))`. There the edge below `A` subtends `{B,C,D,E}` — every other
taxon. Every tree on those taxa contains that split, so counting it would add a
guaranteed-1.00 clade to every result and inflate `fraction_resolved` for free.

Two tests now cover it: one asserting the trivial split is dropped, one asserting
that three different rootings of the same topology yield one identical split set.
The second also pins a related hazard — a rooted tree emits the _same_ bipartition
from both root children, and only returning a `set` stops it being double-counted.

Re-run after the fix: **killed**.

## Not mutated, deliberately

**"Support below a threshold" in place of the topology check.** It is tempting to
test correctness by asserting `fraction_resolved > 0.9` on the easy fixture and
calling it a day. That would pass against a `support()` returning 1.0
unconditionally — mutant 1 — which is why the suite instead asserts that the
false clade scores _lower than every true clade_ on a case where the tree is
genuinely wrong. The comparison is between clades within one run, so no absolute
threshold can satisfy it by accident.

**The `hard_alignment` fixture itself.** `test_a_short_alignment_produces_a_WRONG_tree`
is a guard on the guard: if that assertion ever starts failing, the fixture has
stopped being hard and every support test below it is only measuring the easy
regime, where support is 1.00 everywhere and nothing discriminates.


## Protein alignments

| mutant                                              | result |
| --------------------------------------------------- | ------ |
| moltype never reaches `build_tree` (always dna)      | RED    |
| moltype never reaches `model_finder`                 | RED    |
| parsimony signal always counted in the DNA alphabet  | RED    |
| protein alphabet accepted under the `dna` default    | RED    |

The first two are the ones that matter: they are what "the server accepted
`sequence_type='protein'` and then quietly ran a nucleotide analysis" looks like
from the inside. `test_a_protein_tree_is_built_with_a_protein_model` is only
evidence because of its negative control — a nucleotide model on protein data
RAISES, so building a tree at all proves the molecule type reached the engine.
Without that control the test would pass for any model that happened to work.
