"""Build the 2026-10-24 headline (field) families, with the key fixed before any run.

Three families, 10 questions each (PREREG_2026-10-24_field-families.md):

- ``lba_sister`` (P1): long-branch attraction. 8-20 taxa, two long tip
  branches that are not sisters, short internal edges, simulated by AliSim
  under GTR+G4 with alpha 0.3-0.5. A draw is kept only if the ML tree under
  the generating model AND the ML tree under ModelFinder's BIC pick both
  contain the true cherry, and the JC ML tree does not. The asked taxon is a
  long-branch one; the answer is its cherry partner in the generating tree.
- ``supported_clades`` (P2): 30-100 taxa, some internal edges near zero
  (0.001-0.02), 300-1000 sites, AliSim GTR+G4. The answer is the set of
  clades the agent reports at the 95% level, scored by precision and recall
  against the generating tree. The pass band is the reference pipeline's
  precision and recall minus 3 binomial SE; the reference pipeline is
  IQ-TREE ModelFinder + UFBoot2 (1000 replicates), clades with support >= 95.
- ``align_tree_homology`` (P3): AliSim ``--indel 0.05,0.1`` on 20-50 taxa;
  the agent gets only the unaligned FASTA. (a) the tree, scored by normalised
  RF to the generating tree, passing if <= mean + 3 sd of the reference's
  nRF over 20 seeded runs (MAFFT --auto + ModelFinder + ML); (b) the residue
  of taxon B in the same column as residue i of taxon A, scored exactly
  against AliSim's true alignment (the existing homology rule: no gap in
  either row within 6 columns, the two positions differ by >= 5, and MAFFT
  --auto agrees with the true alignment).

Where the truth comes from: every generating tree is drawn here (``truth.py``,
numpy only) and every sequence is simulated by AliSim (``iqtree3 --alisim``).
IQ-TREE (the CLI, never phylokit-mcp) and MAFFT are used only to CHECK a
draw and to compute the reference pipeline and the shortcut controls, never
to set an answer. This script never imports ``phylokit_mcp``, cogent3 or
piqtree.

No question names a substitution model.

Outputs:
- ``questions.json`` beside this script: the prompts, the data files, the
  truth, the tolerances, the reference and shortcut answers, and the draw
  bookkeeping (drawn / kept / discarded, with the reason).
- ``$PHYLOQA_ARTIFACTS/key/`` (outside the repo): every AliSim output
  (aligned, with gaps), every reference run's tree, by draw; their sha256s
  are written into ``questions.json``.

    IQTREE3=/path/to/iqtree3 MAFFT=/path/to/mafft \\
      HEAVY_RUN_MEM=8G heavy-run uv run python \\
      benchmarks/agent_qa/field_2026-10-24/build_field.py --workers 4
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import truth as T

SEED = 20261024
FAMILY_SEED = {"lba_sister": 1, "supported_clades": 2, "align_tree_homology": 3}
N_PER_FAMILY = 10
ARTIFACTS = Path(
    os.environ.get(
        "PHYLOQA_ARTIFACTS",
        Path.home()
        / "phylokit-mcp-artifacts"
        / "agent_qa"
        / "field-families-2026-10-24",
    )
)
KEY_DIR = ARTIFACTS / "key"
FASTA = "alignment.fasta"
UNALIGNED = "sequences.fasta"

# P2: a clade is reported by the reference pipeline when its UFBoot2 support
# is at least this. 95 is UFBoot2's own recommended threshold for a clade to
# be called supported (Hoang 2018), set here, not taken from the server.
UFBOOT_THRESHOLD = 95.0
UFBOOT_REPLICATES = 1000
# P3: the reference pipeline is re-run with this many seeds; the tree
# tolerance is mean + 3 sd of their nRF to the generating tree.
P3_REF_SEEDS = 20
# P3 homology rule (as in build_questions.py): no gap within this many
# columns of the asked one in either row; the two positions differ by >= 5.
HOMOLOGY_MARGIN = 6
HOMOLOGY_MIN_OFFSET = 5


def tool(env: str, *names: str) -> str:
    """The path of an outside tool: ``$<env>``, else the first name on PATH.
    Fails loud: a key built without its checks is not this key."""
    path = os.environ.get(env) or next(
        (p for p in (shutil.which(n) for n in names) if p), None
    )
    if not path:
        raise SystemExit(f"{names[0]} not found: put it on PATH or set {env}")
    return path


def iqtree_bin() -> str:
    return tool("IQTREE3", "iqtree3")


def mafft_bin() -> str:
    return tool("MAFFT", "mafft")


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# --- outside tools ---------------------------------------------------------------


def read_fasta(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            rows[name] = ""
        elif name is not None:
            rows[name] += line.strip()
    return rows


def _run(cmd: list[str], cwd: Path, timeout: int = 3600) -> str:
    done = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, check=False
    )
    if done.returncode != 0:
        raise RuntimeError(
            f"{' '.join(cmd)} exited {done.returncode}\n{done.stdout[-2000:]}\n"
            f"{done.stderr[-2000:]}"
        )
    return done.stdout


def alisim(
    newick: str,
    model: str,
    length: int,
    seed: int,
    indel: tuple[float, float] | None = None,
) -> dict[str, str]:
    """AliSim's alignment of ``length`` root sites down ``newick``, as rows
    (with ``-`` gaps when there are indels)."""
    with tempfile.TemporaryDirectory(prefix="alisim-") as d:
        tmp = Path(d)
        (tmp / "tree.nwk").write_text(newick + "\n")
        cmd = [
            iqtree_bin(), "--alisim", "sim", "-t", "tree.nwk", "-m", model,
            "--length", str(length), "--seed", str(seed), "-af", "fasta",
            "-T", "1", "-redo",
        ]  # fmt: skip
        if indel is not None:
            cmd += ["--indel", f"{indel[0]},{indel[1]}", "--no-unaligned"]
        _run(cmd, tmp)
        return read_fasta((tmp / "sim.fa").read_text())


def iqtree_ml(
    seqs: dict[str, str],
    model: str,
    seed: int,
    ufboot: int | None = None,
) -> dict[str, Any]:
    """IQ-TREE's ML tree for ``seqs`` under ``model`` ("MFP" runs ModelFinder,
    BIC by default): its Newick (UFBoot supports as internal labels when
    ``ufboot``) and the model it fitted."""
    with tempfile.TemporaryDirectory(prefix="iqtree-") as d:
        tmp = Path(d)
        (tmp / "in.fa").write_text(T.fasta(seqs))
        cmd = [
            iqtree_bin(), "-s", "in.fa", "-m", model, "--seed", str(seed),
            "-T", "1", "--prefix", "run", "-quiet", "-redo",
        ]  # fmt: skip
        if ufboot:
            cmd += ["-B", str(ufboot)]
        _run(cmd, tmp)
        report = (tmp / "run.iqtree").read_text()
        m = re.search(r"Best-fit model according to \w+: (\S+)", report)
        fitted = m.group(1) if m else model
        return {"newick": (tmp / "run.treefile").read_text().strip(), "model": fitted}


def mafft_auto(seqs: dict[str, str]) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="mafft-") as d:
        f = Path(d) / "in.fasta"
        f.write_text(T.fasta(seqs))
        out = _run(
            [mafft_bin(), "--auto", "--quiet", "--thread", "1", "--preservecase",
             "--nuc", str(f)],
            Path(d), timeout=900,
        )  # fmt: skip
    return read_fasta(out)


def mafft_version() -> str:
    done = subprocess.run(
        [mafft_bin(), "--version"], capture_output=True, text=True, check=False
    )
    return (done.stderr + done.stdout).strip().splitlines()[0]


def iqtree_version() -> str:
    out = subprocess.run(
        [iqtree_bin(), "--version"], capture_output=True, text=True, check=False
    ).stdout
    return out.strip().splitlines()[0]


# --- trees and models ------------------------------------------------------------


def names(n: int, rng: np.random.Generator) -> list[str]:
    """Tip labels t1..tn in a random order, so label order says nothing."""
    labels = [f"t{i}" for i in range(1, n + 1)]
    return [labels[i] for i in rng.permutation(n)]


def gtr_model(rng: np.random.Generator, alpha: tuple[float, float]) -> dict[str, Any]:
    """A GTR+F+G4 draw: exchangeabilities around typical empirical values
    (transitions faster), base frequencies from a Dirichlet, alpha uniform."""
    rates = [
        float(rng.uniform(lo, hi))
        for lo, hi in ((0.8, 1.6), (2.0, 5.0), (0.5, 1.2), (0.6, 1.4), (2.5, 6.0))
    ]
    freqs = rng.dirichlet([12, 8, 8, 12])
    a = float(rng.uniform(*alpha))
    spec = (
        "GTR{" + "/".join(f"{r:.3f}" for r in rates) + "}"
        "+F{"
        + "/".join(f"{x:.4f}" for x in freqs[:3])
        + f"/{1 - sum(round(x, 4) for x in freqs[:3]):.4f}"
        "}" + f"+G4{{{a:.3f}}}"
    )
    return {"alisim": spec, "rates": rates, "freqs": [float(x) for x in freqs],
            "alpha": a}  # fmt: skip


def splits_of(newick: str) -> set[T.Split]:
    return T.splits(T.parse_newick(newick))


def has_clade(adj: T.Adj, clade: set[str]) -> bool:
    taxa = set(T.tips(adj))
    want = {frozenset(clade), frozenset(taxa - clade)}
    return bool(want & T.splits(adj))


def cherry_partner(adj: T.Adj, taxon: str) -> str | None:
    for a, b in T.cherries(adj):
        if taxon in (a, b):
            return b if a == taxon else a
    return None


def nrf(a: set[T.Split], b: set[T.Split], n_taxa: int) -> float:
    return T.robinson_foulds(a, b) / (2 * (n_taxa - 3))


def collapse(newick: str, keep: set[T.Split]) -> str:
    """``newick`` (unrooted) rewritten with only the splits in ``keep``: every
    other internal branch contracted to a polytomy. No branch lengths."""
    adj = T.parse_newick(newick)
    taxa = frozenset(T.tips(adj))
    ref = min(taxa)

    def name_of(side: frozenset[str]) -> frozenset[str]:
        return taxa - side if ref in side else side

    # Contract each internal edge whose split is not kept.
    changed = True
    while changed:
        changed = False
        for u, v in sorted(T._edges(adj), key=str):
            if not (isinstance(u, int) and isinstance(v, int)):
                continue
            side = frozenset(T._side(adj, v, u))
            if name_of(side) in keep:
                continue
            for w, ln in list(adj[v].items()):
                if w != u:
                    adj[u][w] = adj[w][u] = ln
                    del adj[w][v]
            del adj[u][v]
            del adj[v]
            changed = True
            break
    root = next(iter(adj[ref]))

    def fmt(node: object, prev: object) -> str:
        if isinstance(node, str):
            return node
        kids = sorted((k for k in adj[node] if k != prev), key=str)
        return "(" + ",".join(fmt(k, node) for k in kids) + ")"

    kids = sorted(adj[root], key=str)
    return "(" + ",".join(fmt(k, root) for k in kids) + ");"


def binomial_se(hits: int, n: int) -> float:
    """SE of a proportion, on the Laplace-adjusted estimate (x+1)/(n+2), so a
    reference at 100% precision still has a band (SE at p=1 is 0)."""
    p = (hits + 1) / (n + 2)
    return math.sqrt(p * (1 - p) / max(n, 1))


# --- draw bookkeeping ------------------------------------------------------------


CACHE = ARTIFACTS / "draw-cache"


def _fingerprint(family: str) -> str:
    """What a cached draw must match: the family's settings and the source of
    its draw functions. A change to a shared helper is not seen here; rebuild
    with --no-cache after one."""
    import inspect

    fns = {"lba_sister": [p1_draw], "supported_clades": [p2_draw, p2_checked],
           "align_tree_homology": [p3_draw]}[family]  # fmt: skip
    params = {"lba_sister": P1, "supported_clades": P2, "align_tree_homology": P3}[
        family
    ]
    src = "".join(inspect.getsource(fn) for fn in fns)
    return sha256(json.dumps(params, sort_keys=True) + src + str(SEED))


def cached_draw(args: tuple[str, int, bool]) -> dict[str, Any]:
    """``FAMILIES[family]``'s draw ``i``, read back from the cache when
    ``use_cache`` and the fingerprint matches (a draw is a pure function of its
    seed). Passed in, not a global: pool workers re-import this module."""
    family, i, use_cache = args
    path = CACHE / family / f"{i}.json"
    fp = _fingerprint(family)
    if use_cache and path.exists():
        saved = json.loads(path.read_text())
        if saved["fingerprint"] == fp:
            return saved["result"]
    res = FAMILIES[family][1](i)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"fingerprint": fp, "result": res}))
    return res


def run_draws(
    family: str,
    draw: Callable[[tuple[str, int, bool]], dict[str, Any]],
    workers: int,
    use_cache: bool = True,
    n_keep: int = N_PER_FAMILY,
    max_draws: int = 400,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Evaluate draws 0, 1, 2, ... (each from its own seed, so in parallel or
    not, the same draws are kept) and keep the first ``n_keep`` that pass, in
    draw order. Returns the kept questions and the bookkeeping."""
    kept: list[dict[str, Any]] = []
    reasons: dict[str, int] = {}
    examined = 0
    i = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        while len(kept) < n_keep:
            if i >= max_draws:
                raise SystemExit(f"{family}: {max_draws} draws, only {len(kept)} kept")
            batch = list(range(i, i + workers))
            i += workers
            for res in pool.map(draw, [(family, j, use_cache) for j in batch]):
                if len(kept) >= n_keep:
                    break  # later draws in the batch are not examined
                examined += 1
                if res["kept"]:
                    kept.append(res["question"])
                else:
                    reasons[res["reason"]] = reasons.get(res["reason"], 0) + 1
                print(f"{family} draw {res['draw']}: "
                      f"{'kept' if res['kept'] else res['reason']}", flush=True)  # fmt: skip
    book = {
        "drawn": examined,
        "kept": len(kept),
        "discarded": examined - len(kept),
        "discard_reasons": reasons,
    }
    return kept, book


def rng_for(family: str, i: int) -> np.random.Generator:
    return np.random.default_rng([SEED, FAMILY_SEED[family], i])


def file_note(name: str, what: str) -> str:
    return f"The {what} is in the file `{name}` in the current working directory."


# --- P1: long-branch attraction --------------------------------------------------

# Fixed after a pilot on draws 1000+ (never used in the key), recorded in the
# PREREG: long tips 0.6-1.2 kept 1 of 8 (the generating model missed 6);
# 1.0-2.0 kept 4 of 8 at 1000-3000 sites; internal 0.01-0.04 kept 7 of 12.
P1 = {
    "taxa": (8, 20),
    "internal": (0.01, 0.04),
    "tip": (0.02, 0.1),
    "long": (1.0, 2.0),
    "sites": (1000, 3000),
    "alpha": (0.3, 0.5),
}
GENERATING_FIT = "GTR+F+G4"


def p1_draw(i: int) -> dict[str, Any]:
    rng = rng_for("lba_sister", i)
    n = int(rng.integers(P1["taxa"][0], P1["taxa"][1] + 1))
    taxa = names(n, rng)
    adj = T.random_tree(taxa, rng, internal=P1["internal"], tip=P1["tip"])
    cherries = T.cherries(adj)
    ia, ib = (int(x) for x in rng.choice(len(cherries), 2, replace=False))
    pa, pb = cherries[ia], cherries[ib]
    k = int(rng.integers(2))
    a, c = pa[k], pa[1 - k]
    k = int(rng.integers(2))
    b = pb[k]
    for x in (a, b):
        (nb,) = adj[x]
        length = round(float(rng.uniform(*P1["long"])), 4)
        adj[x][nb] = adj[nb][x] = length
    model = gtr_model(rng, P1["alpha"])
    sites = int(rng.integers(P1["sites"][0], P1["sites"][1] + 1))
    sim_seed = int(rng.integers(1, 2**31 - 1))
    newick = T.to_newick(adj, rng)
    seqs = alisim(newick, model["alisim"], sites, sim_seed)
    fits = {
        "generating": iqtree_ml(seqs, GENERATING_FIT, 1),
        "bic": iqtree_ml(seqs, "MFP", 1),
        "jc": iqtree_ml(seqs, "JC", 1),
        "no_gamma": iqtree_ml(seqs, "GTR+F", 1),
    }
    has = {k: has_clade(T.parse_newick(v["newick"]), {a, c}) for k, v in fits.items()}
    base = {"draw": i, "has_cherry": has}
    if not has["generating"]:
        return base | {"kept": False, "reason": "generating-model ML misses the cherry"}
    if not has["bic"]:
        return base | {"kept": False, "reason": "BIC-pick ML misses the cherry"}
    if has["jc"]:
        return base | {"kept": False, "reason": "JC ML has the cherry (no LBA)"}
    ref2 = iqtree_ml(seqs, "MFP", 2)
    partner = {
        k: cherry_partner(T.parse_newick(v["newick"]), a)
        for k, v in (fits | {"reference_seed2": ref2}).items()
    }
    aln_text = T.fasta(seqs)
    question = {
        "category": "lba_sister",
        "family": "P1",
        "question": (
            f"These {n} DNA sequences ({sites} sites) evolved along an unknown tree. "
            "Infer the phylogeny and treat it as unrooted. Which taxon forms a "
            f"two-taxon clade (a cherry) with {a}? " + file_note(FASTA, "alignment")
        ),
        "answer_format": "a taxon name, e.g. t3",
        "files": {FASTA: aln_text},
        "truth": {"taxon": c, "taxa": sorted(taxa)},
        "source": (
            f"AliSim {model['alisim']} on a {n}-taxon generating tree with long tips "
            f"{a} and {b} (not sisters); {a}'s cherry partner in that tree"
        ),
        "check": {
            "has_true_cherry": has,
            "partner_of_asked": partner,
            "jc_pairs_the_long_branches": partner["jc"] == b,
            "bic_model": fits["bic"]["model"],
            "reference_seed2_model": ref2["model"],
            "reference_reply": partner["reference_seed2"],
            "shortcut_replies": {
                "JC ML tree": partner["jc"],
                "GTR ML tree, no +G": partner["no_gamma"],
            },
        },
        "generator": {
            "newick": newick,
            "model": model,
            "sites": sites,
            "alisim_seed": sim_seed,
            "long_tips": [a, b],
            "alignment_sha256": sha256(aln_text),
        },
    }
    return base | {"kept": True, "question": question}


# --- P2: reportable clades at 95% ------------------------------------------------

P2 = {
    "taxa": (30, 100),
    "sites": (300, 1000),
    "internal": (0.02, 0.12),
    "tip": (0.03, 0.25),
    # Pilot (never in the key): 0.001-0.02 on 35% of internal edges left the
    # every-clade ML tree's precision only 0.006 and 0.023 under the band's
    # floor (draws 1000-1001); 0.001-0.01 on 50% put it outside the band on 5
    # of 6 draws (1000-1005).
    "near_zero": (0.001, 0.01),
    "near_zero_share": 0.5,
    "alpha": (0.5, 1.5),
}


def labelled_splits(newick: str) -> dict[T.Split, float | None]:
    """Each non-trivial split of ``newick`` (named as in ``truth.splits``) with
    the numeric internal label on its branch (IQ-TREE writes UFBoot support
    there), or None when the branch has no label."""
    toks = re.findall(r"\s*([(),:;]|[^(),:;\s]+)", newick.strip())
    pos = 0
    found: list[tuple[frozenset[str], float | None]] = []
    every: set[str] = set()

    def node() -> frozenset[str]:
        nonlocal pos
        if toks[pos] == "(":
            pos += 1
            tips: set[str] = set()
            while True:
                tips |= node()
                if toks[pos] == ":":
                    pos += 2
                if toks[pos] == ",":
                    pos += 1
                    continue
                pos += 1  # ")"
                break
            label = None
            if pos < len(toks) and toks[pos] not in "(),:;":
                try:
                    label = float(toks[pos].split("/")[0])
                except ValueError:
                    label = None
                pos += 1
            found.append((frozenset(tips), label))
            return frozenset(tips)
        name = toks[pos]
        pos += 1
        every.add(name)
        return frozenset({name})

    node()
    taxa = frozenset(every)
    ref = min(taxa)
    out: dict[T.Split, float | None] = {}
    for side, label in found:
        if ref in side:
            side = taxa - side
        if 2 <= len(side) <= len(taxa) - 2:
            out[side] = label
    return out


def precision_recall(reported: set[T.Split], true: set[T.Split]) -> tuple[float, float]:
    hits = len(reported & true)
    precision = hits / len(reported) if reported else 1.0
    return precision, hits / len(true)


def p2_draw(i: int) -> dict[str, Any]:
    rng = rng_for("supported_clades", i)
    n = int(rng.integers(P2["taxa"][0], P2["taxa"][1] + 1))
    taxa = names(n, rng)
    adj = T.random_tree(taxa, rng, internal=P2["internal"], tip=P2["tip"])
    near = 0
    for u, v in sorted(T._edges(adj), key=str):
        internal = isinstance(u, int) and isinstance(v, int)
        if internal and rng.random() < P2["near_zero_share"]:
            adj[u][v] = adj[v][u] = round(float(rng.uniform(*P2["near_zero"])), 4)
            near += 1
    model = gtr_model(rng, P2["alpha"])
    sites = int(rng.integers(P2["sites"][0], P2["sites"][1] + 1))
    sim_seed = int(rng.integers(1, 2**31 - 1))
    newick = T.to_newick(adj, rng)
    seqs = alisim(newick, model["alisim"], sites, sim_seed)
    true = T.splits(adj)
    runs = {s: iqtree_ml(seqs, "MFP", s, ufboot=UFBOOT_REPLICATES) for s in (1, 2)}
    lab = {s: labelled_splits(r["newick"]) for s, r in runs.items()}
    rep = {
        s: {k for k, v in lab[s].items() if v is not None and v >= UFBOOT_THRESHOLD}
        for s in lab
    }
    p_ref, r_ref = precision_recall(rep[1], true)
    hits = len(rep[1] & true)
    band = {
        "precision_min": round(p_ref - 3 * binomial_se(hits, len(rep[1])), 4),
        "recall_min": round(r_ref - 3 * binomial_se(hits, len(true)), 4),
    }
    ml_all = set(lab[1])
    p_ml, r_ml = precision_recall(ml_all, true)
    p2_, r2_ = precision_recall(rep[2], true)
    aln_text = T.fasta(seqs)
    true_newick = T.to_newick(adj, rng, lengths=False)
    question = {
        "category": "supported_clades",
        "family": "P2",
        "question": (
            f"These {n} DNA sequences ({sites} sites) evolved along an unknown tree. "
            "Treat the tree as unrooted, so a clade here is a split of the taxa into "
            "two groups. Which clades would you report as well supported at the 95% "
            "level? Report every such clade and no others. "
            + file_note(FASTA, "alignment")
        ),
        "answer_format": (
            f"one Newick tree of all {n} taxa, on one line, that contains exactly the "
            "clades you report: collapse every other branch into a polytomy (no "
            "support values or branch lengths needed), e.g. ((t1,t2),t3,(t4,t5,t6),t7);"
        ),
        "files": {FASTA: aln_text},
        "truth": {
            "newick": true_newick,
            "taxa": sorted(taxa),
            "n_true_splits": len(true),
            **band,
        },
        "source": (
            f"AliSim {model['alisim']} on a {n}-taxon generating tree, {near} of its "
            f"{n - 3} internal edges near zero; band = IQ-TREE ModelFinder + UFBoot2 "
            f"({UFBOOT_REPLICATES}) clades >= {UFBOOT_THRESHOLD:g}, precision and "
            "recall minus 3 binomial SE"
        ),
        "check": {
            "reference": {
                "precision": round(p_ref, 4),
                "recall": round(r_ref, 4),
                "n_reported": len(rep[1]),
                "model": runs[1]["model"],
            },
            "reference_seed2": {
                "precision": round(p2_, 4),
                "recall": round(r2_, 4),
                "n_reported": len(rep[2]),
                "model": runs[2]["model"],
            },
            "ml_all_clades": {
                "precision": round(p_ml, 4),
                "recall": round(r_ml, 4),
                "n_reported": len(ml_all),
            },
            "reference_reply": collapse(runs[2]["newick"], rep[2]),
            "shortcut_replies": {
                "ML tree, every clade reported": collapse(runs[1]["newick"], ml_all),
            },
        },
        "generator": {
            "newick": newick,
            "model": model,
            "sites": sites,
            "alisim_seed": sim_seed,
            "near_zero_edges": near,
            "alignment_sha256": sha256(aln_text),
        },
    }
    return {"draw": i, "kept": True, "question": question, "reason": None}


# --- P3: unaligned indel-rich DNA -> tree + homology ------------------------------

P3 = {
    "taxa": (20, 50),
    "root_length": (500, 1000),
    "internal": (0.02, 0.1),
    "tip": (0.05, 0.2),
    "alpha": (0.5, 1.5),
    "indel": (0.05, 0.1),
}


def _drop_gap_columns(rows: dict[str, str]) -> dict[str, str]:
    names_ = list(rows)
    cols = [c for c in zip(*rows.values(), strict=True) if set(c) != {"-"}]
    return {n: "".join(c[k] for c in cols) for k, n in enumerate(names_)}


def pad_unaligned(seqs: dict[str, str]) -> dict[str, str]:
    """The skip-alignment shortcut: every sequence as given, gaps added at the
    end to the longest length."""
    width = max(len(s) for s in seqs.values())
    return {n: s + "-" * (width - len(s)) for n, s in seqs.items()}


def p3_draw(i: int) -> dict[str, Any]:
    rng = rng_for("align_tree_homology", i)
    n = int(rng.integers(P3["taxa"][0], P3["taxa"][1] + 1))
    taxa = names(n, rng)
    adj = T.random_tree(taxa, rng, internal=P3["internal"], tip=P3["tip"])
    model = gtr_model(rng, P3["alpha"])
    length = int(rng.integers(P3["root_length"][0], P3["root_length"][1] + 1))
    sim_seed = int(rng.integers(1, 2**31 - 1))
    newick = T.to_newick(adj, rng)
    aligned = _drop_gap_columns(
        alisim(newick, model["alisim"], length, sim_seed, indel=P3["indel"])
    )
    unaligned = {k: "".join(ch for ch in v if ch != "-") for k, v in aligned.items()}
    true_cols = T.gapped_rows(aligned)
    src, dst = (str(x) for x in rng.choice(sorted(taxa), 2, replace=False))
    n_cols = len(aligned[src])
    m = HOMOLOGY_MARGIN
    ok = [
        c
        for c in range(n_cols // 2, n_cols - m)
        if all(true_cols[t][j] is not None for t in (src, dst) for j in range(c - m, c + m + 1))
        and abs((true_cols[src][c] or 0) - (true_cols[dst][c] or 0)) >= HOMOLOGY_MIN_OFFSET
    ]  # fmt: skip
    base = {"draw": i}
    if not ok:
        return base | {"kept": False, "reason": "no column meets the homology rule"}
    col = ok[int(rng.integers(len(ok)))]
    pos, answer = true_cols[src][col], true_cols[dst][col]
    assert pos is not None and answer is not None
    mafft = mafft_auto(unaligned)
    got = T.homologous_position(T.gapped_rows(mafft), src, pos, dst)
    if got != answer:
        return base | {"kept": False, "reason": "MAFFT --auto disagrees on the column"}
    true = T.splits(adj)
    # ModelFinder once (seed 1); its pick is then fitted with seeds 2-20, so
    # the 20 runs vary the tree search. Restructured after the first build
    # batch: 20 ModelFinder runs per draw projected the build past 2 h.
    first = iqtree_ml(mafft, "MFP", 1)
    refs = [first] + [
        iqtree_ml(mafft, first["model"], s) for s in range(2, P3_REF_SEEDS + 1)
    ]
    ref_nrf = [nrf(splits_of(r["newick"]), true, n) for r in refs]
    tol = float(np.mean(ref_nrf) + 3 * np.std(ref_nrf, ddof=1))
    padded = iqtree_ml(pad_unaligned(unaligned), "MFP", 1)
    pad_nrf = nrf(splits_of(padded["newick"]), true, n)
    text = T.fasta(unaligned)
    true_aln_text = T.fasta(aligned)
    question = {
        "category": "align_tree_homology",
        "family": "P3",
        "question": (
            f"The file `{UNALIGNED}` in the current working directory holds {n} "
            "homologous DNA sequences, unaligned; they differ by substitutions and "
            "by many insertions and deletions. (a) Infer their phylogeny, treated as "
            "unrooted. (b) In a multiple sequence alignment of them, which residue of "
            f"{dst} is in the same column as residue {pos} of {src}? Count residues "
            "from 1, ignoring gaps."
        ),
        "answer_format": (
            "the residue number for (b), then ' | ', then the tree for (a) in Newick "
            f"with all {n} taxa on one line, e.g. 57 | ((t1,t2),t3,(t4,t5));"
        ),
        "files": {UNALIGNED: text},
        "truth": {
            "position": answer,
            "newick": T.to_newick(adj, rng, lengths=False),
            "taxa": sorted(taxa),
            "nrf_max": round(tol, 4),
        },
        "source": (
            f"AliSim {model['alisim']} --indel {P3['indel'][0]},{P3['indel'][1]} on a "
            f"{n}-taxon generating tree; (b) from AliSim's true alignment; (a) nRF "
            f"tolerance = mean + 3 sd over {P3_REF_SEEDS} seeded MAFFT --auto + "
            "ML runs (ModelFinder's pick at seed 1, then fitted at seeds 2-20)"
        ),
        "check": {
            "mafft_auto_agrees": True,
            "offset": answer - pos,
            "reference_nrf": [round(x, 4) for x in ref_nrf],
            "reference_models": sorted({r["model"] for r in refs}),
            "reference_reply": f"{got} | {refs[1]['newick']}",
            "padded_nrf": round(pad_nrf, 4),
            "shortcut_replies": {
                "no alignment (pad at end) + ML": f"{pos} | {padded['newick']}",
                "same relative position + reference tree": (
                    f"{round(pos * len(unaligned[dst]) / len(unaligned[src]))} | "
                    f"{refs[1]['newick']}"
                ),
            },
        },
        "generator": {
            "newick": newick,
            "model": model,
            "root_length": length,
            "alisim_seed": sim_seed,
            "indel": list(P3["indel"]),
            "true_alignment_sha256": sha256(true_aln_text),
            "true_alignment_columns": len(aligned[src]),
        },
        "_artifacts": {
            "true_alignment.fasta": true_aln_text,
            "mafft_auto.fasta": T.fasta(mafft),
        },
    }
    return base | {"kept": True, "question": question}


# --- main ------------------------------------------------------------------------


def p2_checked(i: int) -> dict[str, Any]:
    """``p2_draw`` with P2's two discard rules: the band must not admit the
    star tree (recall_min <= 0, or a reference reporting nothing), and the
    named shortcut, the unbootstrapped ML tree reporting every clade, must
    fall outside it. The second is a shortcut-discrimination filter (the
    PREREG's design note): it is applied in a named regime, near-zero
    internal edges, and its discard count is published with the key. Pilot
    (draws 1000-1005, never in the key): the shortcut passed 1 of 6."""
    res = p2_draw(i)
    q = res["question"]
    t, ml = q["truth"], q["check"]["ml_all_clades"]
    if q["check"]["reference"]["n_reported"] == 0 or t["recall_min"] <= 0:
        return {"draw": i, "kept": False, "reason": "band admits the star tree"}
    if ml["precision"] >= t["precision_min"] and ml["recall"] >= t["recall_min"]:
        return {
            "draw": i,
            "kept": False,
            "reason": "every-clade ML tree is inside the band",
        }
    return res


FAMILIES: dict[str, tuple[str, Callable[[int], dict[str, Any]]]] = {
    "lba_sister": ("p1", p1_draw),
    "supported_clades": ("p2", p2_checked),
    "align_tree_homology": ("p3", p3_draw),
}


def _count(x: float, denom: int) -> int:
    """The integer k of a stored 4-decimal value of k / denom; fails loud if
    the value is not one."""
    k = round(x * denom)
    if abs(k / denom - x) > 1e-4:
        raise SystemExit(f"{x} is not a multiple of 1/{denom}")
    return k


def exact_tolerances(q: dict[str, Any]) -> None:
    """Recompute each tolerance at full precision from integer counts.

    The draws store precision, recall and nRF rounded to 4 decimals, and the
    first build derived the tolerances from those. When all 20 P3 reference
    runs share one nRF, k / 2(n-3), the band's edge rounded below k / 2(n-3),
    and the reference pipeline at seed 2 failed its own key (the positive
    control caught it on 6 of 10). Every stored value is k / m for a known m
    (``_count``); the bands are rebuilt from the k's, unrounded."""
    t, c = q["truth"], q["check"]
    if q["category"] == "align_tree_homology":
        denom = 2 * (len(t["taxa"]) - 3)
        rf = [_count(x, denom) for x in c["reference_nrf"]]
        exact = [k / denom for k in rf]
        c["reference_rf"] = rf
        t["nrf_max"] = float(np.mean(exact) + 3 * np.std(exact, ddof=1))
    elif q["category"] == "supported_clades":
        ref, n_true = c["reference"], t["n_true_splits"]
        hits = _count(ref["recall"], n_true)
        if _count(ref["precision"], ref["n_reported"]) != hits:
            raise SystemExit(f"{q['id']}: precision and recall disagree on the hits")
        precision = hits / ref["n_reported"]
        recall = hits / n_true
        c["reference"]["hits"] = hits
        t["precision_min"] = precision - 3 * binomial_se(hits, ref["n_reported"])
        t["recall_min"] = recall - 3 * binomial_se(hits, n_true)


def write_artifacts(q: dict[str, Any]) -> None:
    """The question's data and AliSim's own outputs, outside the repo; their
    sha256s go into the key."""
    d = KEY_DIR / q["id"]
    d.mkdir(parents=True, exist_ok=True)
    extra = q.pop("_artifacts", {})
    files = {
        **q["files"],
        **extra,
        "generating_tree.nwk": q["generator"]["newick"] + "\n",
    }
    q["generator"]["artifacts_sha256"] = {}
    for name, text in files.items():
        (d / name).write_text(text)
        q["generator"]["artifacts_sha256"][name] = sha256(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--families", default=",".join(FAMILIES))
    ap.add_argument("--out", default=str(HERE / "questions.json"))
    ap.add_argument("--no-cache", action="store_true", help="recompute every draw")
    a = ap.parse_args()
    if not 1 <= a.workers <= 4:
        raise SystemExit("--workers must be 1-4 (shared machine)")
    from importlib.metadata import version

    checked = {
        "iqtree3": iqtree_version(),
        "mafft": mafft_version(),
        # The engine both arms run (repo uv.lock); read from metadata, not imported.
        "piqtree": version("piqtree"),
    }
    questions: list[dict[str, Any]] = []
    draws: dict[str, Any] = {}
    timing: dict[str, float] = {}
    for family in a.families.split(","):
        short = FAMILIES[family][0]
        t0 = time.monotonic()
        kept, book = run_draws(family, cached_draw, a.workers, not a.no_cache)
        timing[family] = round(time.monotonic() - t0, 1)
        for k, q in enumerate(kept, 1):
            q["id"] = f"{short}q{k:02d}"
            exact_tolerances(q)
            write_artifacts(q)
            questions.append(q)
        draws[family] = book
    out = {
        "version": "field-2026-10-24",
        "built": datetime.now(UTC).date().isoformat(),
        "seed": SEED,
        "alt_web": False,
        "checked_with": checked,
        "params": {
            "P1": P1,
            "P1_generating_fit": GENERATING_FIT,
            "P2": P2,
            "P2_ufboot": {
                "replicates": UFBOOT_REPLICATES,
                "threshold": UFBOOT_THRESHOLD,
            },
            "P3": P3,
            "P3_reference_seeds": P3_REF_SEEDS,
            "homology": {"margin": HOMOLOGY_MARGIN, "min_offset": HOMOLOGY_MIN_OFFSET},
        },
        "draws": draws,
        "artifacts": str(KEY_DIR).replace(str(Path.home()), "~"),
        "questions": questions,
    }
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"draws": draws, "seconds": timing}, indent=1))


if __name__ == "__main__":
    main()
