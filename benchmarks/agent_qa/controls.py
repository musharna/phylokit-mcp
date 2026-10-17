"""Controls for the agent benchmark, scored through the real scorer, no LLM.

- Positive control (oracle): the key's own answers, written as replies, must
  all score right; a known-wrong reply to every question must all score wrong.
  ``run.py`` checks both before it starts a sweep, and a test checks them in CI.
- Negative controls (baselines): what a fixed guess or a named shortcut gets
  from the same data files the arms are given.
  - ``best constant``: per kind of question, the one answer that is right most
    often. It is picked WITH the key, so it is an upper bound for any fixed
    guess, not a strategy anyone could follow blind.
  - named shortcuts, each a cheap rule a reader might try instead of the
    intended method: the nearest taxon by JC distance, and the cherry partner
    in a neighbour-joining tree, for the sister-taxon questions; the same
    residue number, or the same relative position, for homology; the
    uncorrected p-distance for JC69 distances; variable sites for
    parsimony-informative sites; a base-composition test plus a transition
    share rule for model choice. Robinson-Foulds questions have no shortcut
    here (two 30-50-tip Newick strings), only the constant.

Like ``truth.py``, this imports no ``phylokit_mcp``, cogent3 or piqtree.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any

import numpy as np
import truth as T
from score import score

MODELS = ["JC69", "K80", "F81", "HKY85"]


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


def key_value(q: dict[str, Any]) -> str:
    """The key's answer, as a reply would write it."""
    t, cat = q["truth"], q["category"]
    if cat == "sister_taxon":
        return t["taxon"]
    if cat == "rf_distance":
        return str(t["rf"])
    if cat == "informative_sites":
        return str(t["count"])
    if cat == "homologous_position":
        return str(t["position"])
    if cat == "model_selection":
        return t["model"]
    if cat == "pairwise_distance":
        return f"{t['distance']:.4f}"
    raise ValueError(f"unknown category {cat!r}")


def wrong_value(q: dict[str, Any]) -> str:
    """A well-formed answer that is not the key's."""
    t, cat = q["truth"], q["category"]
    if cat == "sister_taxon":
        return next(x for x in sorted(t["taxa"]) if x != t["taxon"])
    if cat == "rf_distance":
        return str(t["rf"] + 2)
    if cat == "informative_sites":
        return str(t["count"] + 1)
    if cat == "homologous_position":
        return str(t["position"] + 1)
    if cat == "model_selection":
        return MODELS[(MODELS.index(t["model"]) + 1) % len(MODELS)]
    if cat == "pairwise_distance":
        return f"{t['distance'] + 2 * t['tolerance'] + 0.05:.4f}"
    raise ValueError(f"unknown category {cat!r}")


def reply(value: str) -> str:
    return f"Worked from the file.\nANSWER: {value}"


def oracle_check(questions: list[dict[str, Any]]) -> list[str]:
    """Every key answer scores right and every wrong one scores wrong; returns
    the failures ([] when the key and the scorer agree)."""
    bad = []
    for q in questions:
        got = score(q, reply(key_value(q)))["verdict"]
        if got != "right":
            bad.append(f"{q['id']}: key answer {key_value(q)!r} scored {got}")
        got = score(q, reply(wrong_value(q)))["verdict"]
        if got != "wrong":
            bad.append(f"{q['id']}: wrong answer {wrong_value(q)!r} scored {got}")
    return bad


# --- shortcuts ---------------------------------------------------------------


def _ask(pattern: str, text: str) -> tuple[str, ...]:
    m = re.search(pattern, text)
    if m is None:
        raise ValueError(f"question text does not match {pattern!r}: {text[:120]}")
    return m.groups()


def nearest_by_jc(q: dict[str, Any]) -> str:
    (query,) = _ask(r"cherry\) with (t\d+)\?", q["question"])
    names, d = T.jc_distance_matrix(read_fasta(q["files"]["alignment.fasta"]))
    i = names.index(query)
    row = d[i].copy()
    row[i] = np.inf
    return names[int(np.argmin(row))]


def nj_cherry_partner(q: dict[str, Any]) -> str:
    """The query's partner if the NJ tree on JC distances puts it in a cherry,
    else UNKNOWN."""
    (query,) = _ask(r"cherry\) with (t\d+)\?", q["question"])
    names, d = T.jc_distance_matrix(read_fasta(q["files"]["alignment.fasta"]))
    for a, b in T.cherries(T.neighbor_joining(names, d)):
        if query in (a, b):
            return b if a == query else a
    return "UNKNOWN"


def _homology_ask(q: dict[str, Any]) -> tuple[str, int, str]:
    dst, pos, src = _ask(
        r"which residue of (t\d+) is in the same column as residue (\d+) of (t\d+)\?",
        q["question"],
    )
    return dst, int(pos), src


def same_residue_number(q: dict[str, Any]) -> str:
    return str(_homology_ask(q)[1])


def same_relative_position(q: dict[str, Any]) -> str:
    dst, pos, src = _homology_ask(q)
    seqs = read_fasta(q["files"]["sequences.fasta"])
    return str(round(pos * len(seqs[dst]) / len(seqs[src])))


def uncorrected_p(q: dict[str, Any]) -> str:
    a, b = _ask(r"distance between (t\d+) and (t\d+),", q["question"])
    seqs = read_fasta(q["files"]["alignment.fasta"])
    return f"{T.p_distance(seqs[a], seqs[b]):.4f}"


def variable_sites(q: dict[str, Any]) -> str:
    rows = list(read_fasta(q["files"]["alignment.fasta"]).values())
    return str(sum(len(set(col)) > 1 for col in zip(*rows, strict=True)))


# Chi-square critical value, df = 3, alpha = 0.01.
_CHI2_3_01 = 11.345
# Share of pairwise differences that are transitions: about 1/3 under kappa = 1
# (0.32 with the F81 frequencies used here) and 0.56 under kappa = 2.5 at small
# distances. Halfway-ish, fixed before it was run.
_TS_SHARE = 0.42


def composition_and_ts_rule(q: dict[str, Any]) -> str:
    seqs = list(read_fasta(q["files"]["alignment.fasta"]).values())
    counts = np.array([sum(s.count(b) for s in seqs) for b in T.DNA], dtype=float)
    expected = counts.sum() / 4
    unequal = float(((counts - expected) ** 2 / expected).sum()) > _CHI2_3_01
    ts = tv = 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            for x, y in zip(seqs[i], seqs[j], strict=True):
                if x != y:
                    if frozenset((x, y)) in T._TRANSITIONS:
                        ts += 1
                    else:
                        tv += 1
    kappa = ts / (ts + tv) > _TS_SHARE
    return {(False, False): "JC69", (False, True): "K80", (True, False): "F81"}.get(
        (unequal, kappa), "HKY85"
    )


SHORTCUTS: dict[str, tuple[str, Callable[[dict[str, Any]], str]]] = {
    "nearest taxon by JC distance": ("sister_taxon", nearest_by_jc),
    "neighbour-joining cherry partner": ("sister_taxon", nj_cherry_partner),
    "same residue number": ("homologous_position", same_residue_number),
    "same relative position": ("homologous_position", same_relative_position),
    "uncorrected p-distance": ("pairwise_distance", uncorrected_p),
    "variable sites": ("informative_sites", variable_sites),
    "composition test + transition share": ("model_selection", composition_and_ts_rule),
}


def best_constant(qs: list[dict[str, Any]]) -> tuple[str, list[str]]:
    """The fixed answer right most often across ``qs`` (same kind), and the
    verdicts it gets. Candidates are every key answer in the kind, plus a grid
    for distances; ties go to the first candidate."""
    cands = [key_value(q) for q in qs]
    if qs[0]["category"] == "pairwise_distance":
        cands += [f"{x:.3f}" for x in np.linspace(0.01, 5.0, 500)]
    best, best_v, best_n = cands[0], [], -1
    for c in cands:
        v = [score(q, reply(c))["verdict"] for q in qs]
        n = v.count("right")
        if n > best_n:
            best, best_v, best_n = c, v, n
    return best, best_v


def baselines(questions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One row per (baseline, question): name, category, qid, value, verdict."""
    rows = []
    by_cat: dict[str, list[dict[str, Any]]] = {}
    for q in questions:
        by_cat.setdefault(q["category"], []).append(q)
    for cat, qs in by_cat.items():
        value, verdicts = best_constant(qs)
        for q, v in zip(qs, verdicts, strict=True):
            rows.append(
                {"baseline": "best constant", "category": cat, "qid": q["id"],
                 "value": value, "verdict": v}
            )  # fmt: skip
    for name, (cat, fn) in SHORTCUTS.items():
        for q in by_cat.get(cat, []):
            value = fn(q)
            rows.append(
                {"baseline": name, "category": cat, "qid": q["id"], "value": value,
                 "verdict": score(q, reply(value))["verdict"]}
            )  # fmt: skip
    return rows
