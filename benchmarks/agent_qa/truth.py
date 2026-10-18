"""Ground truth for the agent benchmark, written without the server's libraries.

Nothing here imports ``phylokit_mcp``, cogent3 or piqtree: the trees, the
simulators, the split sets and the counts are numpy and plain Python, so a bug
in the server (or in the libraries it is built on) cannot also be in the key.
``build_questions.py`` uses piqtree and MAFFT only to CHECK that a question is
answerable, never to compute its answer.

Trees are unrooted and held as an adjacency map ``{node: {neighbour: length}}``.
Tips are their names (str); internal nodes are ints.

Substitution models are the HKY85 family, P(t) = expm(Q t) through the
eigendecomposition of the symmetrised rate matrix, with Q scaled to one expected
substitution per unit time:

    Q_ij = pi_j * (kappa if i<->j is a transition else 1),  i != j

JC69 is kappa=1 with equal frequencies, K80 is equal frequencies, F81 is
kappa=1. Protein sequences (the homology questions) evolve under the Poisson
model: every amino acid equally likely, every change equally likely.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable

import numpy as np

DNA = "ACGT"
AMINO = "ACDEFGHIKLMNPQRSTVWY"
_TRANSITIONS = {frozenset("AG"), frozenset("CT")}

Adj = dict[object, dict[object, float]]
Split = frozenset[str]


# --- trees -------------------------------------------------------------------


def random_tree(
    taxa: list[str],
    rng: np.random.Generator,
    internal: tuple[float, float],
    tip: tuple[float, float],
) -> Adj:
    """A random unrooted binary tree by stepwise addition onto a random edge.

    Internal edge lengths are uniform on ``internal``, tip edges on ``tip``.
    """
    if len(taxa) < 3:
        raise ValueError("need at least 3 taxa")
    adj: Adj = {0: {}}
    for t in taxa[:3]:
        adj[t] = {0: 0.0}
        adj[0][t] = 0.0
    nxt = 1
    for t in taxa[3:]:
        edges = sorted(
            {tuple(sorted((u, v), key=str)) for u in adj for v in adj[u]}, key=str
        )
        u, v = edges[int(rng.integers(len(edges)))]
        w = nxt
        nxt += 1
        del adj[u][v], adj[v][u]
        adj[w] = {u: 0.0, v: 0.0, t: 0.0}
        adj[u][w] = adj[v][w] = 0.0
        adj[t] = {w: 0.0}
    # Assign lengths once per undirected edge, in a fixed order.
    for u, v in sorted(_edges(adj), key=str):
        lo, hi = internal if _is_internal(u) and _is_internal(v) else tip
        length = round(float(rng.uniform(lo, hi)), 4)
        adj[u][v] = adj[v][u] = length
    return adj


def _is_internal(node: object) -> bool:
    return isinstance(node, int)


def _edges(adj: Adj) -> set[tuple[object, object]]:
    out = set()
    for u in adj:
        for v in adj[u]:
            out.add(tuple(sorted((u, v), key=str)))
    return out  # type: ignore[return-value]


def tips(adj: Adj) -> list[str]:
    return sorted(n for n in adj if isinstance(n, str))


def _side(adj: Adj, start: object, came_from: object) -> set[str]:
    """Tips reachable from ``start`` without crossing back to ``came_from``."""
    out: set[str] = set()
    stack = [(start, came_from)]
    while stack:
        node, prev = stack.pop()
        if isinstance(node, str):
            out.add(node)
        for nb in adj[node]:
            if nb != prev:
                stack.append((nb, node))
    return out


def splits(adj: Adj) -> set[Split]:
    """Non-trivial unrooted bipartitions, each named by the side WITHOUT the
    alphabetically first taxon (so a split has one name however it is rooted)."""
    taxa = frozenset(tips(adj))
    ref = min(taxa)
    out: set[Split] = set()
    for u, v in _edges(adj):
        side = frozenset(_side(adj, v, u))
        if ref in side:
            side = taxa - side
        if 2 <= len(side) <= len(taxa) - 2:
            out.add(side)
    return out


def robinson_foulds(a: set[Split], b: set[Split]) -> int:
    """Splits in one tree but not the other, summed over both (not halved)."""
    return len(a ^ b)


def path_length(adj: Adj, x: str, y: str) -> float:
    """Sum of branch lengths on the path between two tips."""
    stack: list[tuple[object, object, float]] = [(x, None, 0.0)]
    while stack:
        node, prev, dist = stack.pop()
        if node == y:
            return dist
        for nb, length in adj[node].items():
            if nb != prev:
                stack.append((nb, node, dist + length))
    raise KeyError(f"no path {x} -> {y}")


def cherries(adj: Adj) -> list[tuple[str, str]]:
    """Pairs of tips that share their one neighbour."""
    out = []
    for n in adj:
        if _is_internal(n):
            leaves = sorted(nb for nb in adj[n] if isinstance(nb, str))
            if len(leaves) == 2:
                out.append((leaves[0], leaves[1]))
    return out


def internal_edge_above(adj: Adj, cherry: tuple[str, str]) -> float:
    """Length of the internal edge that separates a cherry from the rest."""
    (parent,) = set(adj[cherry[0]]) & set(adj[cherry[1]])
    (other,) = [nb for nb in adj[parent] if nb not in cherry]
    return adj[parent][other]


def to_newick(
    adj: Adj,
    rng: np.random.Generator,
    *,
    lengths: bool = True,
    root_on_edge: bool = False,
) -> str:
    """One of the many Newick spellings of ``adj``: rooted at a random internal
    node (a trifurcation) or on a random edge, children in random order."""
    internal = sorted(n for n in adj if _is_internal(n))

    def fmt(node: object, prev: object, length: float | None) -> str:
        kids = [nb for nb in adj[node] if nb != prev]
        order = rng.permutation(len(kids))
        if kids:
            body = (
                "("
                + ",".join(fmt(kids[i], node, adj[node][kids[i]]) for i in order)
                + ")"
            )
        else:
            body = str(node)
        if lengths and length is not None:
            body += f":{length:g}"
        return body

    if root_on_edge:
        edges = sorted(_edges(adj), key=str)
        u, v = edges[int(rng.integers(len(edges)))]
        half = adj[u][v] / 2
        parts = [fmt(u, v, round(half, 5)), fmt(v, u, round(adj[u][v] - half, 5))]
        if rng.random() < 0.5:
            parts.reverse()
        return "(" + ",".join(parts) + ");"
    root = internal[int(rng.integers(len(internal)))]
    return fmt(root, None, None) + ";"


_TOKEN = re.compile(r"\s*([(),:;]|[^(),:;\s]+)")


def parse_newick(text: str) -> Adj:
    """A small Newick reader for the strings this module writes (and other
    plain ones): names on tips only, optional ``:length``. A degree-2 root is
    suppressed, so the result is the unrooted tree."""
    toks = _TOKEN.findall(text.strip())
    pos = 0
    adj: Adj = {}
    counter = [0]

    def new_internal() -> int:
        counter[0] += 1
        adj[counter[0]] = {}
        return counter[0]

    def length() -> float:
        nonlocal pos
        if pos < len(toks) and toks[pos] == ":":
            value = float(toks[pos + 1])
            pos += 2
            return value
        return 0.0

    def node() -> object:
        nonlocal pos
        if toks[pos] == "(":
            pos += 1
            me = new_internal()
            while True:
                child = node()
                ln = length()
                adj[me][child] = adj[child][me] = ln
                if toks[pos] == ",":
                    pos += 1
                    continue
                if toks[pos] != ")":
                    raise ValueError(f"expected ) at token {pos}: {toks[pos]!r}")
                pos += 1
                break
            # An internal label (support value) is ignored.
            if pos < len(toks) and toks[pos] not in "(),:;":
                pos += 1
            return me
        name = toks[pos]
        pos += 1
        if name in adj:
            raise ValueError(f"duplicate tip {name!r}")
        adj[name] = {}
        return name

    root = node()
    length()
    if pos >= len(toks) or toks[pos] != ";":
        raise ValueError("Newick must end with ;")
    if _is_internal(root) and len(adj[root]) == 2:
        a, b = adj[root]
        joined = adj[root][a] + adj[root][b]
        del adj[a][root], adj[b][root], adj[root]
        adj[a][b] = adj[b][a] = joined
    return adj


# --- substitution ------------------------------------------------------------


def rate_matrix(kappa: float, freqs: Iterable[float]) -> np.ndarray:
    pi = np.asarray(list(freqs), dtype=float)
    q = np.zeros((4, 4))
    for i, a in enumerate(DNA):
        for j, b in enumerate(DNA):
            if i != j:
                q[i, j] = pi[j] * (kappa if frozenset((a, b)) in _TRANSITIONS else 1.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    q /= -float(np.dot(pi, np.diag(q)))  # one expected substitution per unit time
    return q


def transition_matrix(q: np.ndarray, pi: np.ndarray, t: float) -> np.ndarray:
    """expm(Q t) for a reversible Q, via the symmetric form D^1/2 Q D^-1/2."""
    d = np.sqrt(pi)
    s = (d[:, None] * q) / d[None, :]
    s = (s + s.T) / 2
    w, v = np.linalg.eigh(s)
    p = (v * np.exp(w * t)) @ v.T
    p = (p / d[:, None]) * d[None, :]
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=1, keepdims=True)


def poisson_matrix(n_states: int, t: float) -> np.ndarray:
    same = 1 / n_states + (1 - 1 / n_states) * math.exp(-n_states * t / (n_states - 1))
    p = np.full((n_states, n_states), (1 - same) / (n_states - 1))
    np.fill_diagonal(p, same)
    return p


def _evolve_states(
    states: np.ndarray, p: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    cum = np.cumsum(p, axis=1)
    u = rng.random(len(states))
    out = (u[:, None] > cum[states]).sum(axis=1)
    return np.minimum(out, p.shape[0] - 1)


def simulate(
    adj: Adj,
    n_sites: int,
    rng: np.random.Generator,
    *,
    kappa: float = 1.0,
    freqs: tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25),
) -> dict[str, str]:
    """An HKY-family alignment evolved down ``adj`` (rooted at its first
    internal node; the models are reversible, so the root does not matter)."""
    pi = np.asarray(freqs, dtype=float)
    q = rate_matrix(kappa, pi)
    root = min(n for n in adj if _is_internal(n))
    seqs: dict[object, np.ndarray] = {root: rng.choice(4, size=n_sites, p=pi)}
    order = [(root, None)]
    while order:
        node, prev = order.pop(0)
        for nb in sorted(adj[node], key=str):
            if nb == prev:
                continue
            p = transition_matrix(q, pi, adj[node][nb])
            seqs[nb] = _evolve_states(seqs[node], p, rng)
            order.append((nb, node))
    return {
        t: "".join(DNA[i] for i in seqs[t])
        for t in sorted(n for n in seqs if isinstance(n, str))
    }


def simulate_with_indels(
    adj: Adj,
    root_length: int,
    rng: np.random.Generator,
    *,
    alphabet: str,
    indel_rate: float,
    max_indel: int = 6,
) -> tuple[dict[str, str], dict[str, list[int | None]]]:
    """Sequences with substitutions AND indels, plus the true alignment.

    Every residue carries the id of its column; insertions create new ids, so
    homology is known exactly. Returns ``(unaligned, true_alignment)`` where
    the true alignment maps each tip to a list over columns of either the
    residue's 1-based position in that tip's sequence or None (a gap).
    Substitutions are Poisson over ``alphabet``; indel events per branch are
    Poisson(indel_rate * length * sequence length), half insertions and half
    deletions, of geometric length capped at ``max_indel``.
    """
    k = len(alphabet)
    root = min(n for n in adj if _is_internal(n))
    order_ids: list[int] = list(range(root_length))
    next_id = root_length
    seqs: dict[object, list[tuple[int, int]]] = {
        root: [(i, int(s)) for i, s in enumerate(rng.integers(0, k, root_length))]
    }
    queue = [(root, None)]
    while queue:
        node, prev = queue.pop(0)
        for nb in sorted(adj[node], key=str):
            if nb == prev:
                continue
            t = adj[node][nb]
            parent = seqs[node]
            states = _evolve_states(
                np.array([s for _, s in parent]), poisson_matrix(k, t), rng
            )
            child = [(cid, int(s)) for (cid, _), s in zip(parent, states, strict=True)]
            for _ in range(int(rng.poisson(indel_rate * t * len(child)))):
                size = min(int(rng.geometric(0.4)), max_indel)
                if rng.random() < 0.5 and len(child) > size + 10:
                    start = int(rng.integers(0, len(child) - size + 1))
                    del child[start : start + size]
                else:
                    at = int(rng.integers(0, len(child) + 1))
                    new = [(next_id + j, int(rng.integers(0, k))) for j in range(size)]
                    next_id += size
                    # Place the new columns right after the left neighbour's
                    # column in the global order (or before the first residue).
                    if at == 0:
                        idx = order_ids.index(child[0][0]) if child else 0
                    else:
                        idx = order_ids.index(child[at - 1][0]) + 1
                    order_ids[idx:idx] = [cid for cid, _ in new]
                    child[at:at] = new
            seqs[nb] = child
            queue.append((nb, node))
    tip_names = sorted(n for n in seqs if isinstance(n, str))
    present = {cid for t in tip_names for cid, _ in seqs[t]}
    columns = [cid for cid in order_ids if cid in present]
    unaligned = {t: "".join(alphabet[s] for _, s in seqs[t]) for t in tip_names}
    true_aln: dict[str, list[int | None]] = {}
    for t in tip_names:
        pos = {cid: i + 1 for i, (cid, _) in enumerate(seqs[t])}
        true_aln[t] = [pos.get(cid) for cid in columns]
    return unaligned, true_aln


def homologous_position(
    aln: dict[str, list[int | None]], src: str, pos: int, dst: str
) -> int | None:
    """The residue of ``dst`` in the same column as residue ``pos`` of ``src``."""
    col = aln[src].index(pos)
    return aln[dst][col]


def gapped_rows(aln_rows: dict[str, str]) -> dict[str, list[int | None]]:
    """A gapped alignment (``-`` for gaps) as residue positions per column."""
    out: dict[str, list[int | None]] = {}
    for name, row in aln_rows.items():
        k = 0
        cols: list[int | None] = []
        for ch in row:
            if ch == "-":
                cols.append(None)
            else:
                k += 1
                cols.append(k)
        out[name] = cols
    return out


# --- statistics ----------------------------------------------------------------


def parsimony_informative(seqs: dict[str, str]) -> int:
    """Sites where at least two states each occur in at least two sequences."""
    rows = list(seqs.values())
    count = 0
    for col in zip(*rows, strict=True):
        tally: dict[str, int] = {}
        for ch in col:
            tally[ch] = tally.get(ch, 0) + 1
        if sum(1 for v in tally.values() if v >= 2) >= 2:
            count += 1
    return count


def p_distance(a: str, b: str) -> float:
    return sum(x != y for x, y in zip(a, b, strict=True)) / len(a)


def jc_distance(p: float) -> float:
    return -0.75 * math.log(1 - 4 * p / 3)


def jc_distance_sd(d: float, n_sites: int) -> float:
    """Large-sample SD of the JC69 distance estimate at true distance ``d``:
    Var = p(1-p) / (n (1 - 4p/3)^2), with p the expected p-distance."""
    p = 0.75 * (1 - math.exp(-4 * d / 3))
    return math.sqrt(p * (1 - p) / n_sites) / (1 - 4 * p / 3)


def neighbor_joining(names: list[str], dist: np.ndarray) -> Adj:
    """Saitou & Nei neighbour joining; branch lengths are not needed, only the
    topology, so negative lengths are kept as computed."""
    nodes: list[object] = list(names)
    d = dist.astype(float).copy()
    adj: Adj = {n: {} for n in names}
    nxt = 0
    while len(nodes) > 3:
        n = len(nodes)
        r = d.sum(axis=1)
        q = (n - 2) * d - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        i, j = np.unravel_index(int(np.argmin(q)), q.shape)
        new = nxt
        nxt += 1
        di = 0.5 * d[i, j] + (r[i] - r[j]) / (2 * (n - 2))
        adj[new] = {nodes[i]: di, nodes[j]: d[i, j] - di}
        adj[nodes[i]][new] = di
        adj[nodes[j]][new] = d[i, j] - di
        dk = 0.5 * (d[i] + d[j] - d[i, j])
        keep = [k for k in range(n) if k not in (i, j)]
        d = np.vstack(
            [np.hstack([d[np.ix_(keep, keep)], dk[keep][:, None]]), [*dk[keep], 0.0]]
        )
        nodes = [nodes[k] for k in keep] + [new]
    centre = nxt
    adj[centre] = {}
    for k, node in enumerate(nodes):
        a, b = (x for x in range(3) if x != k)
        length = 0.5 * (d[k, a] + d[k, b] - d[a, b])
        adj[centre][node] = length
        adj[node][centre] = length
    return adj


def jc_distance_matrix(seqs: dict[str, str]) -> tuple[list[str], np.ndarray]:
    names = sorted(seqs)
    m = np.zeros((len(names), len(names)))
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            p = min(p_distance(seqs[names[a]], seqs[names[b]]), 0.74)
            m[a, b] = m[b, a] = jc_distance(p)
    return names, m


def fasta(seqs: dict[str, str]) -> str:
    return "".join(f">{n}\n{s}\n" for n, s in seqs.items())
