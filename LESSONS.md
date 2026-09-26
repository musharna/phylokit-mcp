# Lessons

One line per miss: the class, and the mechanism that now catches it.

- 2026-09-22 — A reproducibility claim was tested on the coarsest output (support) and hid drift in branch lengths and lnL; claims are now tested on every number a caller would compare, across fresh processes, with a frozen-clock positive control (`tests/test_support.py`).
- 2026-09-22 — The refusal wrapper converts only `ValueError`, and cogent3/piqtree report bad input as other types, so user errors surfaced as opaque crashes; errors are now classified at the input-parsing boundary, and `tests/test_input_boundary.py` asserts every known library input error carries its reason through a real MCP session.
- 2026-09-22 — `validate()` kept its own copy of the alphabet and it drifted from cogent3's; it is now derived from cogent3, and a test compares the two on every printable character.
- 2026-09-22 — A test asserted `ready_for_infer_tree is True` for output `infer_tree` refuses, encoding the bug; "ready" is now computed by `infer_tree`'s own checks.
- 2026-09-25 — `rank_models` trusted that a dict iterated once per model, but piqtree re-adds the best models under non-string keys, so the winner was ranked as its own runner-up; names are collapsed at the boundary with a disagreement check, and `tests/test_models.py::test_each_model_is_ranked_once_and_never_its_own_runner_up` checks uniqueness on real engine output.
