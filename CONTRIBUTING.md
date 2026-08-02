# Contributing

Thanks for helping improve `phylokit-mcp`. It is an MCP server exposing IQ-TREE 2 phylogenetic inference via piqtree, and it never returns a topology without its support.

## Dev setup

Requires Python >=3.12 and [`uv`](https://docs.astral.sh/uv/).

Python 3.12+ is a hard floor: piqtree publishes wheels for 3.12 and up, and it
is the inference engine. No system packages are needed — piqtree ships IQ-TREE 2
as a wheel, so there is no separate binary to install and no version skew between
a CLI and the library.

```bash
uv sync
```

## Running the tests

```bash
uv run pytest -q
uv run ruff check .
uv run ruff format --check .
```

The suite runs real inference. `docs/mutate.sh` drives the mutation checks.

## Testing rules

These are the rules this repo is actually held to. They share one idea: **a test is
worth only what it can fail on.**

### 1. Never trust a test you have not seen fail

Before a test is trusted, run it against the broken state — for a bug fix, the
pre-fix code you still have — and confirm it fails *for the stated reason*, not
merely that it fails. A test asserting "some exception was raised" passes against
broken code too.

### 2. A negative result needs a positive control

A test asserting that something is refused must also assert, **in the same test**,
that the legitimate path still succeeds. Otherwise a harness that raises on
everything reads as "the guard works". `tests/test_registry_metadata.py` is the
worked example: it rejects an over-long `server.json` description *and* validates
the real document in the same function.

### 3. Prefer an external oracle to self-comparison

A test that compares the server against itself will keep passing when the server is
consistently wrong. Where a ground truth exists — a known geometry, a simulated
tree, an analytic expectation — assert against that instead. See `docs/EVAL.md`.

### 4. Mutation checks

`docs/MUTATION-CHECKS.md` records deliberate mutations introduced to confirm the
suite catches them. **A surviving mutant is the coverage report.** If you change
behaviour in a guarded area, add the mutation you used to prove the guard works.

## Pull requests

- Update `CHANGELOG.md` for any user-facing change (Keep a Changelog format).
- Update the README if tool signatures or configuration change.
- Bump nothing else: version lives in `pyproject.toml`, `server.json` (x3) and
  `CITATION.cff`, and a release PR moves them together. `tests/test_registry_metadata.py`
  enforces that they agree, and validates `server.json` against the MCP registry's
  published schema so a constraint violation fails here rather than after a tag is cut.
- Fail loud. Do not add silent fallbacks or swallow errors; surface enough context
  (inputs, what was attempted) that the failure is debuggable.
