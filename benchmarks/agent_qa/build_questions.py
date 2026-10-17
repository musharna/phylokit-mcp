"""Generate the agent benchmark's questions, with answers fixed before any run.

Every answer comes from a known generator in ``truth.py``: a tree drawn at
random, sequences evolved down it, a count or a distance computed on the data
or the tree. This script never imports ``phylokit_mcp``, so a bug in the
server cannot also be in the key.

Two outside tools are used, and only to CHECK that a question is answerable,
never to set its answer:

- piqtree (the upstream engine): the maximum-likelihood JC69 tree must have the
  generating topology (sister questions), its tip-to-tip distance must fall in
  the tolerance (distance questions), and the generating model must beat the
  other three candidates by at least 2 BIC units (model questions).
- MAFFT: its alignment must put the asked residue in the same column as the
  true alignment does (homology questions).

A draw that fails its check is discarded and redrawn; the number of discards
is written to questions.json, so the selection is visible.

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
# comparison: kappa (1) and base frequencies (3).
MODELS: dict[str, dict[str, Any]] = {
    "JC69": {"kappa": 1.0, "freqs": (0.25, 0.25, 0.25, 0.25), "k": 0, "iqtree": "JC"},
    "K80": {"kappa": 6.0, "freqs": (0.25, 0.25, 0.25, 0.25), "k": 1, "iqtree": "K2P"},
    "F81": {"kappa": 1.0, "freqs": (0.4, 0.1, 0.1, 0.4), "k": 3, "iqtree": "F81"},
    "HKY85": {"kappa": 6.0, "freqs": (0.4, 0.1, 0.1, 0.4), "k": 4, "iqtree": "HKY"},
}
# Every name a correct answer might use for each model. IQ-TREE spells F81 and
# HKY with "+F", and JC+F / K2P+F are its names for the same two models.
MODEL_ALIASES = {
    "JC69": ["JC69", "JC", "Jukes-Cantor", "Jukes Cantor"],
    "K80": ["K80", "K2P", "Kimura 2-parameter", "Kimura two-parameter", "Kimura 1980"],
    "F81": ["F81", "F81+F", "JC+F", "JC69+F", "Felsenstein 1981", "Felsenstein 81"],
    "HKY85": ["HKY85", "HKY", "HKY+F", "HKY85+F", "K2P+F", "K80+F"],
}

discards: dict[str, int] = {}


def _discard(category: str) -> None:
    discards[category] = discards.get(category, 0) + 1


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


def data_block(label: str, text: str) -> str:
    return f"{label} (also saved as `{FASTA}` in the current directory):\n\n{text}"


# --- categories ----------------------------------------------------------------


def sister_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    while len(out) < n:
        taxa = names(8, rng)
        adj = T.random_tree(taxa, rng, internal=(0.06, 0.15), tip=(0.05, 0.25))
        seqs = T.simulate(adj, 600, rng)
        nm, d = T.jc_distance_matrix(seqs)
        nj_rf = T.robinson_foulds(T.splits(T.neighbor_joining(nm, d)), T.splits(adj))
        ml_rf = T.robinson_foulds(T.splits(ml_tree(seqs, "JC")[0]), T.splits(adj))
        if nj_rf or ml_rf:
            _discard("sister_taxon")
            continue
        pair = T.cherries(adj)[int(rng.integers(len(T.cherries(adj))))]
        ask, answer = (pair[0], pair[1]) if rng.random() < 0.5 else (pair[1], pair[0])
        out.append(
            {
                "category": "sister_taxon",
                "question": (
                    "These 8 DNA sequences evolved along an unknown tree. Infer the "
                    "maximum-likelihood phylogeny under the JC69 model and treat it as "
                    f"unrooted. Which taxon forms a two-taxon clade (a cherry) with {ask}?"
                    "\n\n" + data_block("Alignment", T.fasta(seqs))
                ),
                "answer_format": "a taxon name, e.g. t3",
                "files": {FASTA: T.fasta(seqs)},
                "truth": {"taxon": answer, "taxa": sorted(taxa)},
                "source": (
                    "JC69 simulation (truth.py) on a random 8-taxon tree; the asked "
                    "taxon's cherry partner in the generating tree"
                ),
                "check": {"nj_rf_to_truth": nj_rf, "piqtree_ml_rf_to_truth": ml_rf},
                "generator": {"newick": T.to_newick(adj, rng), "n_sites": 600},
            }
        )
    return out


def rf_questions(rng: np.random.Generator, n: int) -> list[dict[str, Any]]:
    moves = [0, 1, 2, 2, 3, 4, 5, 7][:n]
    out: list[dict[str, Any]] = []
    for k in moves:
        n_taxa = int(rng.integers(8, 13))
        taxa = names(n_taxa, rng)
        a = T.random_tree(taxa, rng, internal=(0.02, 0.2), tip=(0.05, 0.4))
        b = {u: dict(v) for u, v in a.items()}
        for _ in range(k):
            nni(b, rng)
        nw_a = T.to_newick(a, rng, lengths=True, root_on_edge=bool(rng.random() < 0.5))
        nw_b = T.to_newick(b, rng, lengths=bool(rng.random() < 0.5), root_on_edge=True)
        # The key is computed from the strings in the question, not the
        # objects they were written from.
        rf = T.robinson_foulds(
            T.splits(T.parse_newick(nw_a)), T.splits(T.parse_newick(nw_b))
        )
        out.append(
            {
                "category": "rf_distance",
                "question": (
                    f"Here are two phylogenies of the same {n_taxa} taxa, in Newick "
                    "(also saved as `trees.nwk`, one tree per line). Treat both as "
                    "unrooted and ignore branch lengths. What is the Robinson-Foulds "
                    "distance between them, counted as the number of non-trivial "
                    "bipartitions (splits) found in one tree but not the other, summed "
                    "over both trees (not halved, not normalised)?\n\n"
                    f"Tree A: {nw_a}\nTree B: {nw_b}"
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
    out: list[dict[str, Any]] = []
    kinds = ["dna", "protein"] * (n // 2) + ["dna"] * (n % 2)
    for moltype in kinds:
        while True:
            taxa = names(6, rng)
            adj = T.random_tree(taxa, rng, internal=(0.03, 0.1), tip=(0.03, 0.12))
            alphabet = T.DNA if moltype == "dna" else T.AMINO
            length = int(
                rng.integers(150, 221) if moltype == "dna" else rng.integers(100, 161)
            )
            unaligned, aln = T.simulate_with_indels(
                adj, length, rng, alphabet=alphabet, indel_rate=0.04
            )
            src, dst = (str(x) for x in rng.choice(sorted(unaligned), 2, replace=False))
            n_cols = len(aln[src])
            # A column is askable when the true alignment has no gap in ANY
            # sequence within 8 columns of it, and the two positions differ
            # (so the answer is not just the question's number back).
            ok = [
                c
                for c in range(8, n_cols - 8)
                if all(aln[t][j] is not None for t in aln for j in range(c - 8, c + 9))
                and aln[src][c] != aln[dst][c]
            ]
            if not ok:
                _discard("homologous_position")
                continue
            col = ok[int(rng.integers(len(ok)))]
            pos, answer = aln[src][col], aln[dst][col]
            assert pos is not None and answer is not None
            got = T.homologous_position(
                T.gapped_rows(mafft_rows(unaligned, moltype)), src, pos, dst
            )
            if got != answer:
                _discard("homologous_position")
                continue
            break
        kind = "DNA" if moltype == "dna" else "protein"
        out.append(
            {
                "category": "homologous_position",
                "question": (
                    f"These 6 homologous {kind} sequences are unaligned; they differ by "
                    "substitutions and by insertions and deletions. In a multiple sequence "
                    f"alignment of them, which residue of {dst} is in the same column as "
                    f"residue {pos} of {src}? Count residues from 1, ignoring gaps."
                    "\n\n"
                    + "Sequences (also saved as `sequences.fasta` in the current directory):"
                    + f"\n\n{T.fasta(unaligned)}"
                ),
                "answer_format": f"an integer: a residue position in {dst}, e.g. 57",
                "files": {"sequences.fasta": T.fasta(unaligned)},
                "truth": {"position": answer, "moltype": moltype},
                "source": (
                    "substitution + indel simulation (truth.py) that tracks every "
                    "residue's column; no gap in any sequence within 8 columns"
                ),
                "check": {"mafft_auto_agrees": True},
            }
        )
    return out


def model_questions(rng: np.random.Generator, per_model: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    n_sites = 1000
    for model, spec in MODELS.items():
        got = 0
        while got < per_model:
            taxa = names(6, rng)
            adj = T.random_tree(taxa, rng, internal=(0.05, 0.12), tip=(0.05, 0.15))
            seqs = T.simulate(
                adj, n_sites, rng, kappa=spec["kappa"], freqs=spec["freqs"]
            )
            bic = {
                m: s["k"] * float(np.log(n_sites)) - 2 * ml_tree(seqs, s["iqtree"])[1]
                for m, s in MODELS.items()
            }
            ranked = sorted(bic, key=bic.__getitem__)
            margin = bic[ranked[1]] - bic[ranked[0]]
            if ranked[0] != model or margin < 2:
                _discard("model_selection")
                continue
            got += 1
            out.append(
                {
                    "category": "model_selection",
                    "question": (
                        "Consider four DNA substitution models, all without rate "
                        "heterogeneity or invariant sites: JC69, K80 (K2P), F81 and "
                        "HKY85. With the tree topology and branch lengths estimated by "
                        "maximum likelihood, which of the four has the lowest BIC on "
                        "this alignment of 6 sequences?\n\n"
                        + data_block("Alignment", T.fasta(seqs))
                    ),
                    "answer_format": "one of JC69, K80, F81, HKY85",
                    "files": {FASTA: T.fasta(seqs)},
                    "truth": {"model": model, "aliases": MODEL_ALIASES},
                    "source": (
                        f"{model} simulation (truth.py; kappa={spec['kappa']}, "
                        f"freqs={list(spec['freqs'])}) on a random 6-taxon tree"
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
    n_sites = 1000
    # Each question's tree gets its own scale, spaced evenly in log from 0.3
    # to 4 and shuffled, so the true distances spread from about 0.1 to 1.5
    # and no single guessed number lands near most of them.
    scales = list(np.geomspace(0.3, 4.0, n)[rng.permutation(n)])
    while len(out) < n:
        taxa = names(6, rng)
        scale = float(scales[len(out)])
        adj = T.random_tree(
            taxa,
            rng,
            internal=(0.02 * scale, 0.08 * scale),
            tip=(0.03 * scale, 0.15 * scale),
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
        if abs(pairwise - d) > tol or abs(ml - d) > tol:
            _discard("pairwise_distance")
            continue
        out.append(
            {
                "category": "pairwise_distance",
                "question": (
                    "These 6 DNA sequences evolved along an unknown tree under the JC69 "
                    f"model. Estimate the evolutionary distance between {x} and {y}, in "
                    "expected substitutions per site.\n\n"
                    + data_block("Alignment", T.fasta(seqs))
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
        n_taxa = int(rng.integers(6, 11))
        n_sites = int(rng.integers(150, 301))
        taxa = names(n_taxa, rng)
        adj = T.random_tree(taxa, rng, internal=(0.02, 0.1), tip=(0.02, 0.15))
        seqs = T.simulate(adj, n_sites, rng)
        out.append(
            {
                "category": "informative_sites",
                "question": (
                    f"How many parsimony-informative sites does this alignment of "
                    f"{n_taxa} DNA sequences ({n_sites} sites) have? A site is "
                    "parsimony-informative when at least two different bases each "
                    "occur in at least two of the sequences.\n\n"
                    + data_block("Alignment", T.fasta(seqs))
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
        *sister_questions(rng, 8),
        *rf_questions(rng, 8),
        *homology_questions(rng, 8),
        *model_questions(rng, 2),
        *distance_questions(rng, 8),
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
        "built": datetime.now(UTC).date().isoformat(),
        "seed": SEED,
        "checked_with": {
            "piqtree": piqtree.__version__,
            "mafft": (done.stderr + done.stdout).strip().splitlines()[0],
        },
        "discarded_draws": discards,
        "questions": questions,
    }
    (HERE / "questions.json").write_text(json.dumps(out, indent=1) + "\n")
    counts: dict[str, int] = {}
    for q in questions:
        counts[q["category"]] = counts.get(q["category"], 0) + 1
    print(len(questions), counts, "discarded:", discards)
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
