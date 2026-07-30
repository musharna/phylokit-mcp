"""The IQ-TREE binding, and the two things that make its answers reproducible.

Threads are pinned before piqtree is imported. IQ-TREE parallelises the
likelihood evaluation, and floating-point addition is not associative, so the
order in which partial likelihoods are summed changes the last bits of the
log-likelihood — which changes which of two near-tied topologies wins. A run
that is not thread-pinned is not reproducible even with a fixed seed.
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
    "phylokit-mcp needs piqtree, the official IQ-TREE 2 Python bindings.\n"
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
    return str(getattr(piqtree(), "__version__", "unknown"))


def threads_pinned() -> bool:
    return os.environ.get("OMP_NUM_THREADS") == "1"


# Measured on piqtree 0.8.3, not assumed. Passing the same `rand_seed` does NOT
# fully reset IQ-TREE's internal state within a process: building the same tree
# three times in a row gave call 1 == call 2 but call 3 different. Across FRESH
# PROCESSES the same call sequence reproduces exactly — three runs of an
# identical 30-replicate bootstrap returned byte-identical support.
#
# The practical size of the drift, measured over six repeated 50-replicate calls
# in one process: three of four clades bit-identical, one moved 0.02 — a single
# replicate flipping. That is well inside the bootstrap's own sampling error
# (about 0.07 at 50 replicates), so the topology and every conclusion drawn from
# it were unchanged. This is reported rather than papered over because the server
# is long-lived by design, which is exactly the condition that exposes it.
IN_PROCESS_DRIFT = 0.05

REPRODUCIBILITY_NOTE = (
    "Identical across fresh processes: the same request in a new server process "
    "reproduces exactly. Within ONE long-lived process, a repeated identical "
    f"request may move a support value by up to about {IN_PROCESS_DRIFT:.2f}, "
    "because IQ-TREE's internal state is not fully reset by rand_seed. That is "
    "smaller than the bootstrap's own sampling error, so it does not change the "
    "topology or which clades are supported."
)


def reproducibility() -> dict:
    return {
        "deterministic_across_processes": True,
        "bit_exact_on_repeat_within_process": False,
        "max_support_drift_within_process": IN_PROCESS_DRIFT,
        "threads_pinned": threads_pinned(),
        "note": REPRODUCIBILITY_NOTE,
    }
