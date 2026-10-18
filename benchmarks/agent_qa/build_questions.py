"""Generate the agent benchmark's questions, with answers fixed before any run.

Every answer comes from a known generator in ``truth.py``: a tree drawn at
random, sequences evolved down it, a count or a distance computed on the data
or the tree. This script never imports ``phylokit_mcp``, so a bug in the
server cannot also be in the key.

The inputs are sized past what can be worked by eye (16-24 taxa, 800-1200
sites; 30-50-tip trees; ten 300-700-residue sequences with many indels) and
are given to the agent as FILES only, never in the prompt. Version 1 inlined
6-8-taxon data in the prompt, and Claude with no tools answered 28 of 30
right: the sweep measured nothing about the server.

Two outside tools are used, and only to CHECK that a question is answerable
from its data, never to set its answer:

- piqtree (the upstream engine): the maximum-likelihood JC69 tree must contain
  the asked cherry (sister questions), its tip-to-tip distance must fall in the
  tolerance (distance questions), and the generating model must beat the other
  three candidates by at least 2 BIC units (model questions).
- MAFFT: its alignment must put the asked residue in the same column as the
  true alignment does (homology questions).

A draw that fails its check is discarded and redrawn; the number of draws and
discards per category is written to questions.json, so the selection is
visible.

    uv run python benchmarks/agent_qa/build_questions.py   # writes questions.json
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import truth as T

SEED = 20261016
HERE = Path(__file__).resolve().parent
FASTA = "alignment.fasta"

# Model questions: the four nested HKY-family candidates. Free parameters
# beyond the branch lengths, which all four share and so cancel out of a BIC
# comparison: kappa (1) and base frequencies (3). The departures from JC69 are
# mild (kappa 2.5; frequencies 0.32/0.18) so that the answer takes a fit, not
# a glance at the base counts.
MODELS: dict[str, dict[str, Any]] = {
    "JC69": {"kappa": 1.0, "freqs": (0.25, 0.25, 0.25, 0.25), "k": 0, "iqtree": "JC"},
    "K80": {"kappa": 2.5, "freqs": (0.25, 0.25, 0.25, 0.25), "k": 1, "iqtree": "K2P"},
    "F81": {"kappa": 1.0, "freqs": (0.32, 0.18, 0.18, 0.32), "k": 3, "iqtree": "F81"},
    "HKY85": {"kappa": 2.5, "freqs": (0.32, 0.18, 0.18, 0.32), "k": 4, "iqtree": "HKY"},
}
# Every name a correct answer might use for each model. IQ-TREE spells F81 and
# HKY with "+F", and JC+F / K2P+F are its names for the same two models.
MODEL_ALIASES = {
    "JC69": ["JC69", "JC", "Jukes-Cantor", "Jukes Cantor"],
    "K80": ["K80", "K2P", "Kimura 2-parameter", "Kimura two-parameter", "Kimura 1980"],
    "F81": ["F81", "F81+F", "JC+F", "JC69+F", "Felsenstein 1981", "Felsenstein 81"],
    "HKY85": ["HKY85", "HKY", "HKY+F", "HKY85+F", "K2P+F", "K80+F"],
}

draws: dict[str, dict[str, int]] = {}


def _draw(category: str, kept: bool) -> None:
    d = draws.setdefault(category, {"drawn": 0, "discarded": 0})
    d["drawn"] += 1
    d["discarded"] += 0 if kept else 1


def names(n: int, rng: np.random.Generator) -> list[str]:
    """Tip labels t1..tn in a random order, so label order says nothing about
    the tree."""
    labels = [f"t{i}" for i in range(1, n + 1)]
    return [labels[i] for i in rng.permutation(n)]


def _piqtree():
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import piqtree
    from cogent3 import make_aligned_seqs

    return piqtree, make_aligned_seqs


def ml_tree(seqs: dict[str, str], model: str) -> tuple[T.Adj, float]:
    """piqtree's ML tree, read back through truth.parse_newick, and its lnL."""
    piqtree, make_aligned_seqs = _piqtree()
    tree = piqtree.build_tree(
        make_aligned_seqs(seqs, moltype="dna"), model, rand_seed=1
    )
    adj = T.parse_newick(tree.get_newick(with_distances=True))
    return adj, float(tree.params["lnL"])


def has_clade(adj: T.Adj, clade: set[str]) -> bool:
    """Whether the unrooted tree has the split separating ``clade``."""
    taxa = set(T.tips(adj))
    want = {frozenset(clade), frozenset(taxa - clade)}
    return bool(want & T.splits(adj))


def file_note(name: str, what: str) -> str:
    return f"The {what} is in the file `{name}` in the current working directory."


# --- categories ----------------------------------------------------------------


def sister_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    """24 taxa, 800 sites, short internal edges (0.01-0.05: the hard side of
    docs/EVAL.md's measured regimes) under long tips. A draw is kept only if
    piqtree's ML tree contains the asked cherry."""
    out: list[dict[str, Any]] = []
    n_taxa, n_sites = 24, 800
    while len(out) < n:
        taxa = names(n_taxa, rng)
        adj = T.random_tree(taxa, rng, internal=(0.01, 0.05), tip=(0.05, 0.3))
        seqs = T.simulate(adj, n_sites, rng)
        cherries = T.cherries(adj)
        pair = cherries[int(rng.integers(len(cherries)))]
        ask, answer = (pair[0], pair[1]) if rng.random() < 0.5 else (pair[1], pair[0])
        ml = ml_tree(seqs, "JC")[0]
        kept = has_clade(ml, {ask, answer})
        _draw("sister_taxon", kept)
        if not kept:
            continue
        nm, d = T.jc_distance_matrix(seqs)
        nj = T.neighbor_joining(nm, d)
        i = nm.index(ask)
        nearest = nm[min((j for j in range(len(nm)) if j != i), key=lambda j: d[i, j])]
        out.append(
            {
                "category": "sister_taxon",
                "question": (
                    f"These {n_taxa} DNA sequences ({n_sites} sites) evolved along an "
                    "unknown tree. Infer the maximum-likelihood phylogeny under the JC69 "
                    "model and treat it as unrooted. Which taxon forms a two-taxon clade "
                    f"(a cherry) with {ask}? " + file_note(FASTA, "alignment")
                ),
                "answer_format": "a taxon name, e.g. t3",
                "files": {FASTA: T.fasta(seqs)},
                "truth": {"taxon": answer, "taxa": sorted(taxa)},
                "source": (
                    f"JC69 simulation (truth.py) on a random {n_taxa}-taxon tree; the "
                    "asked taxon's cherry partner in the generating tree"
                ),
                "check": {
                    "piqtree_ml_has_cherry": True,
                    "piqtree_ml_rf_to_truth": T.robinson_foulds(
                        T.splits(ml), T.splits(adj)
                    ),
                    "nj_has_cherry": has_clade(nj, {ask, answer}),
                    "nearest_by_jc_distance_is_sister": nearest == answer,
                },
                "generator": {"newick": T.to_newick(adj, rng), "n_sites": n_sites},
            }
        )
    return out


def rf_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    moves = [3, 6, 10, 15, 20, 30][:n]
    out: list[dict[str, Any]] = []
    for k in moves:
        n_taxa = int(rng.integers(30, 51))
        taxa = names(n_taxa, rng)
        a = T.random_tree(taxa, rng, internal=(0.02, 0.2), tip=(0.05, 0.4))
        b = {u: dict(v) for u, v in a.items()}
        for _ in range(k):
            nni(b, rng)
        nw_a = T.to_newick(a, rng, lengths=True, root_on_edge=bool(rng.random() < 0.5))
        nw_b = T.to_newick(b, rng, lengths=bool(rng.random() < 0.5), root_on_edge=True)
        # The key is computed from the strings in the file, not the objects
        # they were written from.
        rf = T.robinson_foulds(
            T.splits(T.parse_newick(nw_a)), T.splits(T.parse_newick(nw_b))
        )
        _draw("rf_distance", True)
        out.append(
            {
                "category": "rf_distance",
                "question": (
                    f"The file `trees.nwk` in the current working directory holds two "
                    f"phylogenies of the same {n_taxa} taxa in Newick, one per line (tree "
                    "A, then tree B). Treat both as unrooted and ignore branch lengths. "
                    "What is the Robinson-Foulds distance between them, counted as the "
                    "number of non-trivial bipartitions (splits) found in one tree but not "
                    "the other, summed over both trees (not halved, not normalised)?"
                ),
                "answer_format": "an integer, e.g. 6",
                "files": {"trees.nwk": f"{nw_a}\n{nw_b}\n"},
                "truth": {"rf": rf, "max_rf": 2 * (n_taxa - 3)},
                "source": "truth.splits on both Newick strings as written",
                "generator": {"nni_moves": k, "n_taxa": n_taxa},
            }
        )
    return out


def nni(adj: T.Adj, rng: np.random.Generator) -> None:
    """One nearest-neighbour interchange across a random internal edge."""
    internal = sorted(
        e
        for e in (tuple(sorted(x, key=str)) for x in _pairs(adj))
        if all(isinstance(z, int) for z in e)
    )
    u, v = internal[int(rng.integers(len(internal)))]
    a = sorted((x for x in adj[u] if x != v), key=str)[int(rng.integers(2))]
    b = sorted((x for x in adj[v] if x != u), key=str)[int(rng.integers(2))]
    la, lb = adj[u].pop(a), adj[v].pop(b)
    del adj[a][u], adj[b][v]
    adj[u][b] = adj[b][u] = lb
    adj[v][a] = adj[a][v] = la


def _pairs(adj: T.Adj) -> set[tuple[object, object]]:
    return {(u, v) for u in adj for v in adj[u]}


def mafft_rows(seqs: dict[str, str], moltype: str) -> dict[str, str]:
    exe = shutil.which("mafft")
    if exe is None:
        raise SystemExit("MAFFT is not on PATH; the homology check needs it.")
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "in.fasta"
        f.write_text(T.fasta(seqs))
        flag = "--nuc" if moltype == "dna" else "--amino"
        done = subprocess.run(
            [exe, "--auto", "--quiet", "--thread", "1", "--preservecase", flag, str(f)],
            capture_output=True,
            text=True,
            timeout=300,
            check=True,
        )
    rows: dict[str, str] = {}
    name = None
    for line in done.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            rows[name] = ""
        elif name is not None:
            rows[name] += line.strip()
    return rows


def homology_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    """Ten sequences, 500-700 nt or 300-450 aa, indel-rich. The asked column
    is in the second half (so many indels lie upstream of it), the two
    positions differ by at least 5, neither sequence has a gap within 6
    columns of it, and MAFFT agrees with the true alignment."""
    out: list[dict[str, Any]] = []
    n_seqs = 10
    kinds = ["dna", "protein"] * (n // 2) + ["dna"] * (n % 2)
    for moltype in kinds:
        while True:
            taxa = names(n_seqs, rng)
            adj = T.random_tree(taxa, rng, internal=(0.03, 0.1), tip=(0.05, 0.15))
            alphabet = T.DNA if moltype == "dna" else T.AMINO
            length = int(
                rng.integers(500, 701) if moltype == "dna" else rng.integers(300, 451)
            )
            unaligned, aln = T.simulate_with_indels(
                adj, length, rng, alphabet=alphabet, indel_rate=0.08
            )
            src, dst = (str(x) for x in rng.choice(sorted(unaligned), 2, replace=False))
            n_cols = len(aln[src])
            ok = [
                c
                for c in range(n_cols // 2, n_cols - 6)
                if all(
                    aln[t][j] is not None
                    for t in (src, dst)
                    for j in range(c - 6, c + 7)
                )
                and abs((aln[src][c] or 0) - (aln[dst][c] or 0)) >= 5
            ]
            if not ok:
                _draw("homologous_position", False)
                continue
            col = ok[int(rng.integers(len(ok)))]
            pos, answer = aln[src][col], aln[dst][col]
            assert pos is not None and answer is not None
            got = T.homologous_position(
                T.gapped_rows(mafft_rows(unaligned, moltype)), src, pos, dst
            )
            kept = got == answer
            _draw("homologous_position", kept)
            if kept:
                break
        kind = "DNA" if moltype == "dna" else "protein"
        out.append(
            {
                "category": "homologous_position",
                "question": (
                    f"The file `sequences.fasta` in the current working directory holds "
                    f"{n_seqs} homologous {kind} sequences, unaligned; they differ by "
                    "substitutions and by many insertions and deletions. In a multiple "
                    f"sequence alignment of them, which residue of {dst} is in the same "
                    f"column as residue {pos} of {src}? Count residues from 1, ignoring "
                    "gaps."
                ),
                "answer_format": f"an integer: a residue position in {dst}, e.g. 57",
                "files": {"sequences.fasta": T.fasta(unaligned)},
                "truth": {"position": answer, "moltype": moltype},
                "source": (
                    "substitution + indel simulation (truth.py) that tracks every "
                    "residue's column; no gap in either sequence within 6 columns"
                ),
                "check": {"mafft_auto_agrees": True, "offset": answer - pos},
            }
        )
    return out


def model_questions(rng: np.random.Generator, per_model: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    n_taxa, n_sites = 16, 1200
    for model, spec in MODELS.items():
        got = 0
        while got < per_model:
            taxa = names(n_taxa, rng)
            adj = T.random_tree(taxa, rng, internal=(0.02, 0.08), tip=(0.03, 0.12))
            seqs = T.simulate(
                adj, n_sites, rng, kappa=spec["kappa"], freqs=spec["freqs"]
            )
            bic = {
                m: s["k"] * float(np.log(n_sites)) - 2 * ml_tree(seqs, s["iqtree"])[1]
                for m, s in MODELS.items()
            }
            ranked = sorted(bic, key=bic.__getitem__)
            margin = bic[ranked[1]] - bic[ranked[0]]
            kept = ranked[0] == model and margin >= 2
            _draw("model_selection", kept)
            if not kept:
                continue
            got += 1
            out.append(
                {
                    "category": "model_selection",
                    "question": (
                        "Consider four DNA substitution models, all without rate "
                        "heterogeneity or invariant sites: JC69, K80 (K2P), F81 and "
                        "HKY85. With the tree topology and branch lengths estimated by "
                        "maximum likelihood, which of the four has the lowest BIC on this "
                        f"alignment of {n_taxa} sequences ({n_sites} sites)? "
                        + file_note(FASTA, "alignment")
                    ),
                    "answer_format": "one of JC69, K80, F81, HKY85",
                    "files": {FASTA: T.fasta(seqs)},
                    "truth": {"model": model, "aliases": MODEL_ALIASES},
                    "source": (
                        f"{model} simulation (truth.py; kappa={spec['kappa']}, "
                        f"freqs={list(spec['freqs'])}) on a random {n_taxa}-taxon tree"
                    ),
                    "check": {
                        "piqtree_bic_order": ranked,
                        "bic_margin_over_runner_up": round(margin, 2),
                    },
                }
            )
    return out


def distance_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    n_taxa, n_sites = 16, 1000
    # Each question's tree gets its own scale, spaced evenly in log from 0.3
    # to 4 and shuffled, so the true distances spread widely and no single
    # guessed number lands near most of them.
    scales = list(np.geomspace(0.3, 4.0, n)[rng.permutation(n)])
    while len(out) < n:
        taxa = names(n_taxa, rng)
        scale = float(scales[len(out)])
        adj = T.random_tree(
            taxa,
            rng,
            internal=(0.01 * scale, 0.04 * scale),
            tip=(0.03 * scale, 0.12 * scale),
        )
        seqs = T.simulate(adj, n_sites, rng)
        pairs = [(a, b) for a in sorted(taxa) for b in sorted(taxa) if a < b]
        cherry = {tuple(sorted(c)) for c in T.cherries(adj)}
        far = [p for p in pairs if p not in cherry]
        x, y = far[int(rng.integers(len(far)))]
        d = T.path_length(adj, x, y)
        # Tolerance: 3 large-sample SDs of the JC69 distance estimate at the
        # true distance. The pairwise estimate and the ML tree's tip-to-tip
        # distance both have to land inside it, or the draw is discarded.
        tol = 3 * T.jc_distance_sd(d, n_sites)
        pairwise = T.jc_distance(T.p_distance(seqs[x], seqs[y]))
        ml = T.path_length(ml_tree(seqs, "JC")[0], x, y)
        kept = abs(pairwise - d) <= tol and abs(ml - d) <= tol
        _draw("pairwise_distance", kept)
        if not kept:
            continue
        out.append(
            {
                "category": "pairwise_distance",
                "question": (
                    f"These {n_taxa} DNA sequences ({n_sites} sites) evolved along an "
                    f"unknown tree under the JC69 model. Estimate the evolutionary "
                    f"distance between {x} and {y}, in expected substitutions per site. "
                    + file_note(FASTA, "alignment")
                ),
                "answer_format": "a number of substitutions per site, e.g. 0.123",
                "files": {FASTA: T.fasta(seqs)},
                "truth": {"distance": round(d, 4), "tolerance": round(tol, 4)},
                "source": (
                    "path length between the two tips in the generating tree "
                    "(JC69 simulation, truth.py); tolerance 3 SD of the JC69 estimate"
                ),
                "check": {
                    "pairwise_jc_estimate": round(pairwise, 4),
                    "piqtree_ml_patristic": round(ml, 4),
                },
            }
        )
    return out


def informative_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for _ in range(n):
        n_taxa = int(rng.integers(16, 25))
        n_sites = int(rng.integers(800, 1001))
        taxa = names(n_taxa, rng)
        adj = T.random_tree(taxa, rng, internal=(0.01, 0.05), tip=(0.01, 0.08))
        seqs = T.simulate(adj, n_sites, rng)
        _draw("informative_sites", True)
        out.append(
            {
                "category": "informative_sites",
                "question": (
                    f"How many parsimony-informative sites does this alignment of "
                    f"{n_taxa} DNA sequences ({n_sites} sites) have? A site is "
                    "parsimony-informative when at least two different bases each "
                    "occur in at least two of the sequences. "
                    + file_note(FASTA, "alignment")
                ),
                "answer_format": "an integer, e.g. 42",
                "files": {FASTA: T.fasta(seqs)},
                "truth": {"count": T.parsimony_informative(seqs)},
                "source": "truth.parsimony_informative on the alignment as written",
            }
        )
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    questions = [
        *sister_questions(rng, 6),
        *rf_questions(rng, 6),
        *homology_questions(rng, 6),
        *model_questions(rng, 2),
        *distance_questions(rng, 6),
        *informative_questions(rng, 6),
    ]
    piqtree, _ = _piqtree()
    for i, q in enumerate(questions, 1):
        q["id"] = f"q{i:02d}"
    done = subprocess.run(
        [shutil.which("mafft") or "mafft", "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    out = {
        "version": 2,
        "built": datetime.now(UTC).date().isoformat(),
        "seed": SEED,
        "checked_with": {
            "piqtree": piqtree.__version__,
            "mafft": (done.stderr + done.stdout).strip().splitlines()[0],
        },
        "draws": draws,
        "questions": questions,
    }
    (HERE / "questions.json").write_text(json.dumps(out, indent=1) + "\n")
    counts: dict[str, int] = {}
    for q in questions:
        counts[q["category"]] = counts.get(q["category"], 0) + 1
    print(len(questions), counts)
    print("draws:", draws)
    sizes = sorted(len(t) for q in questions for t in q["files"].values())
    print(f"file sizes: {sizes[0]}-{sizes[-1]} bytes")
    # What a constant guess would score on the one tolerance-scored category:
    # the best single number, tried at every true value.
    dist = [q["truth"] for q in questions if q["category"] == "pairwise_distance"]
    best = max(
        sum(abs(g["distance"] - t["distance"]) <= t["tolerance"] for t in dist)
        for g in dist
    )
    print(
        f"pairwise_distance: the best constant guess is right on {best} of {len(dist)}"
    )


if __name__ == "__main__":
    main()
