"""The IQ-TREE binding, and what can and cannot be made reproducible about it.

Threads are pinned before piqtree is imported. IQ-TREE parallelises the
likelihood evaluation, and floating-point addition is not associative, so the
order in which partial likelihoods are summed changes the last bits of the
log-likelihood — which changes which of two near-tied topologies wins. A run
that is not thread-pinned is not reproducible even with a fixed seed. Pinning is
necessary but NOT sufficient -- see the measurement below.
"""

from __future__ import annotations

import os
from typing import Any

# Set BEFORE the import below. IQ-TREE reads these when its thread pool
# initialises, which happens at import, so setting them afterwards is a no-op
# that looks like it worked.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

_INSTALL_HINT = (
    "phylokit-mcp needs piqtree, the official IQ-TREE Python bindings (IQ-TREE 3 in piqtree 0.8).\n"
    "Install it with:  pip install 'piqtree>=0.8,<0.9'\n"
    "piqtree publishes prebuilt wheels for Python 3.12+ on Linux and macOS; if "
    "pip is trying to build from source, check that your Python is 3.12 or newer."
)

_piqtree: Any = None


def piqtree() -> Any:
    """The piqtree module, with an install error that names the fix.

    Imported through a function rather than at module scope so that an import
    sorter cannot hoist it above the thread pinning above, which would silently
    remove the reproducibility guarantee this module exists to provide.
    """
    global _piqtree
    if _piqtree is None:
        try:
            import piqtree as _p
        except ImportError as exc:  # pragma: no cover - environment-dependent
            raise RuntimeError(_INSTALL_HINT) from exc
        _piqtree = _p
    return _piqtree


def engine_version() -> str:
    """piqtree's version. Not IQ-TREE's: see `iqtree_version`."""
    return str(getattr(piqtree(), "__version__", "unknown"))


def iqtree_version() -> str:
    """The version of the IQ-TREE build compiled into piqtree (3.x for 0.8)."""
    return str(getattr(piqtree(), "__iqtree_version__", "unknown"))


def threads_pinned() -> bool:
    return os.environ.get("OMP_NUM_THREADS") == "1"


# Measured on piqtree 0.8.3, not assumed. The same request with the same
# `rand_seed` does NOT reproduce bit-exactly -- neither on repeat within one
# process NOR across fresh processes. Eight fresh processes building the same
# ML tree gave five distinct log-likelihoods (spread ~2e-6) and branch lengths
# differing from the 7th significant digit.
#
# The cause is IQ-TREE reading the wall clock during the search. Freezing
# gettimeofday() alone (an LD_PRELOAD shim; freezing time(), clock() or
# getrusage() instead changed nothing) made every run bit-identical, in-process
# and across processes. The thread count was already 1 (piqtree's default, and
# OMP_NUM_THREADS above), and the seed does reach the engine, so neither of
# those is the lever. piqtree exposes no option to take the clock out of the
# search, so this server cannot make the result exact and does not claim to.
# An earlier version reported deterministic_across_processes=True on the
# strength of a test that compared support values only -- the one output
# coarse enough to hide the drift.
#
# The practical size: branch lengths and the log-likelihood move in the
# trailing digits, and the bootstrap replicates are subject to the same drift,
# so a support value can move by a replicate flipping. Over six repeated
# 50-replicate calls, three of four clades were bit-identical and one moved
# 0.02 -- well inside the bootstrap's own sampling error (about 0.07 at 50
# replicates). The column resampling itself is numpy-seeded and exact.
SUPPORT_DRIFT = 0.05
# Kept under its old name because the test suite imports it. The within-process
# bound is unchanged; it is now known to apply across processes too.
IN_PROCESS_DRIFT = SUPPORT_DRIFT

# UFBoot2 support moves much more, and not by a replicate flipping. IQ-TREE
# draws the replicates itself and scores them on the trees its clock-dependent
# search happens to visit, so the support is a different estimate on every run.
# The 0.05 above was measured on the nonparametric bootstrap and was reported
# for UFBoot calls too. Measured on the 60-site test fixture, seed 7, 1000
# replicates: 48 fresh processes gave one clade list and 48 distinct results;
# the weakest clade ranged 0.22 to 0.34 (sd 0.025 per run), the others 0.08,
# 0.04 and 0.02; 12 repeats in one process ranged 0.07. With the clock frozen
# (the gettimeofday shim) every run was identical, so the clock is the whole
# cause. Two calls differ with sd about 0.035, which puts 0.15 about four sd out.
UFBOOT_SUPPORT_DRIFT = 0.15

_RESAMPLING = {
    "bootstrap": (
        "That is a replicate flipping, smaller than the bootstrap's own "
        "sampling error. The seed does fix the column resampling exactly."
    ),
    "ufboot": (
        "IQ-TREE draws the UFBoot2 replicates inside its search, so the seed "
        "does not pin the support either. A clade near the 0.95 line can be "
        "called supported on one call and not on the next."
    ),
}


def support_drift(support_method: str) -> float:
    """The support drift this server claims for `support_method`."""
    if support_method == "bootstrap":
        return SUPPORT_DRIFT
    if support_method == "ufboot":
        return UFBOOT_SUPPORT_DRIFT
    raise ValueError(f"no measured support drift for support_method {support_method!r}")


def reproducibility_note(support_method: str) -> str:
    drift = support_drift(support_method)
    return (
        "Not bit-exact, in any setting: IQ-TREE reads the wall clock during its "
        "search, so the same request with the same seed -- repeated in this "
        "process or sent to a fresh one -- can return branch lengths and a "
        "log-likelihood that differ in the trailing digits, and a support value "
        f"that differs by up to about {drift:.2f}. "
        f"{_RESAMPLING[support_method]} Compare trees with compare_trees, and "
        "numbers with a tolerance, never by string equality."
    )


SEED_MIN = 0
# IQ-TREE takes the seed as a C int; numpy's generator needs it non-negative.
SEED_MAX = 2**31 - 1


def validate_seed(seed: int) -> None:
    """Refuse a seed the engine cannot take, BEFORE any work is done.

    Outside this range the failure was a pybind TypeError (masked as a crash)
    or, for a negative seed, numpy's ValueError -- raised by the bootstrap only
    after the full maximum-likelihood search had already run.
    """
    if isinstance(seed, bool) or not SEED_MIN <= seed <= SEED_MAX:
        raise ValueError(
            f"seed must be an integer between {SEED_MIN} and {SEED_MAX} "
            f"(IQ-TREE takes a 32-bit signed seed), got {seed!r}."
        )


def reproducibility(support_method: str) -> dict:
    """What a repeat of this call can and cannot reproduce, for its method."""
    drift = support_drift(support_method)
    return {
        "deterministic_across_processes": False,
        "bit_exact_on_repeat_within_process": False,
        "max_support_drift_within_process": drift,
        "max_support_drift_across_processes": drift,
        "bootstrap_resampling_exact": support_method == "bootstrap",
        "threads_pinned": threads_pinned(),
        "note": reproducibility_note(support_method),
    }
