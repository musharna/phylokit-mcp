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


## align_sequences — 2026-09-18, at 0.5.0

Against MAFFT v7.526. Each was also run once with the failure printed, to confirm
the test failed for the reason in the right-hand column and not some other one.

| #   | mutant                                                  | killed by                                                        |
| --- | ------------------------------------------------------- | ---------------------------------------------------------------- |
| 11  | a single sequence is sent to MAFFT                      | `DID NOT RAISE AlignmentError` — MAFFT exits 0 on one sequence   |
| 12  | alphabet check off                                      | `DID NOT RAISE AlignmentError`                                   |
| 13  | gapped input accepted                                   | the degap check fires instead: MAFFT stripped the gaps silently  |
| 14  | MAFFT always told the input is nucleotide               | `'--nuc' == '--amino'`                                           |
| 15  | output rows not checked against the input               | `DID NOT RAISE` on a stand-in that drops a residue               |
| 16  | `--preservecase` dropped                                | 8 of 9 tests: MAFFT lower-cases, the degap check rejects it      |
| 17  | non-zero exit ignored                                   | error text no longer carries `exited 3` and the stderr           |
| 18  | missing MAFFT raises `RuntimeError`                     | the refusal type, and the text no longer survives `call_tool`    |
| 19  | output replaced by right-padded input                   | planted-indel columns wrong; end-to-end tree at RF 4             |
| 20  | output name/order not checked                           | `DID NOT RAISE` on a stand-in that reorders                      |
| 21  | two aligned sequences reported ready for `infer_tree`   | `ready_for_infer_tree is False`                                  |
| 22  | no wall-clock cap                                       | the hung stand-in runs to completion; wrong error                |
| 23  | per-sequence length cap off                             | `DID NOT RAISE AlignmentError`                                   |

Mutant 19 is the one worth reading. Equal-length rows are all `infer_tree` checks
for, and right-padding produces them. Only a control whose homology was planted
can tell that from an alignment.

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
