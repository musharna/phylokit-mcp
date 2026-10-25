"""Score one benchmark answer against its question's truth.

Every run ends its reply with ``ANSWER: <value>`` (or ``ANSWER: UNKNOWN``);
the verdict is decided from that line alone, so the same reply always scores
the same way. Verdicts:

- ``right``: the value matches the truth under the category's rule.
- ``wrong``: a value that does not.
- ``declined``: ``ANSWER: UNKNOWN``.
- ``no_answer``: no ``ANSWER:`` line at all (a format failure, kept apart so it
  is never silently counted as either of the others).

Rules by category, and why each is as strict as it is:

- ``sister_taxon``: the taxon name, without case. The data were checked to
  give the generating tree under both neighbour joining and ML, so there is
  one right name.
- ``rf_distance``, ``informative_sites``: the integer, exactly. Both are
  deterministic functions of the text in the question.
- ``homologous_position``: the integer, exactly. Positions were chosen with no
  gap in either sequence within 6 columns, and MAFFT agrees with the true
  alignment on every one, so the homology is not in doubt.
- ``model_selection``: any name of the right model (``aliases`` in the key
  lists every name of all four; IQ-TREE writes HKY85 as ``HKY+F`` and F81 as
  ``F81+F``).
- ``pairwise_distance``: within the question's ``tolerance``, 3 large-sample
  SDs of the JC69 distance estimate at the true distance. An estimate from the
  data (pairwise or tip-to-tip on the ML tree) differs from the generating
  value by sampling error; exact equality would score correct work as wrong.

Field families (``field_2026-10-24/``, PREREG_2026-10-24_field-families.md):

- ``lba_sister``: as ``sister_taxon`` (draws were kept only where the ML tree
  under the generating model, and under ModelFinder's BIC pick, has the
  cherry and the JC tree does not).
- ``supported_clades``: the reply is a Newick tree holding exactly the clades
  reported; its non-trivial splits are scored against the generating tree's.
  Right iff precision >= ``precision_min`` and recall >= ``recall_min`` (the
  reference pipeline's minus 3 binomial SE, fixed at build time). A tree that
  does not hold exactly the alignment's taxa is wrong. Nothing reported has
  precision 1 and recall 0.
- ``align_tree_homology``: ``<residue> | <Newick>``. Right iff the residue
  equals the true alignment's, exactly, AND the tree's normalised RF to the
  generating tree is <= ``nrf_max`` (the reference's mean nRF + max(3 sd,
  one split), 2/(2n-6)). Each part's result is in ``detail``.

A wrong answer that is still a legitimate option (a taxon in the alignment,
one of the four candidate models) is also flagged ``wrong_id``: a well-formed
choice that is not the answer, as opposed to something malformed. The flag
says nothing about why it was chosen.
"""

from __future__ import annotations

import re
from typing import Any

import truth as T

_ANSWER_RE = re.compile(r"^[\s>*_`#-]*answer[*_`]*\s*:\s*(.*?)\s*$", re.IGNORECASE)
_UNKNOWN = {"UNKNOWN", "UNKNOWN.", "N/A", "NONE"}
# A whole number that is not part of a decimal or of a name: "6", "RF = 6",
# "6 splits"; not the 0 or the 5 of "0.5", and not the 6 of "t6" (taxa are
# t1..tn, so "residue 338 of t6" and "t6: 338" must both read as 338).
_INT = re.compile(r"(?<![\w.])(\d[\d,]*)(?![\w.]|,\d)")
# A number that is not part of a name: the 10 of "d(t10, t9) = 0.44" is not it.
_NUMBER = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def extract(reply: str) -> str | None:
    """The value on the reply's last ``ANSWER:`` line, unwrapped; None if none."""
    found = None
    for line in reply.splitlines():
        m = _ANSWER_RE.match(line)
        if m:
            found = m.group(1)
    if found is None:
        return None
    # Some replies wrap the line in <answer>...</answer>; the tags are not the value.
    return re.sub(r"</?answer>", "", found, flags=re.IGNORECASE).strip(" \t*_`.")


def _integer(value: str) -> int | None:
    m = _INT.search(value)
    return int(m.group(1).replace(",", "")) if m else None


def _number(value: str) -> float | None:
    m = _NUMBER.search(value)
    return float(m.group(0)) if m else None


def _model_key(name: str) -> str:
    """Case, spaces and hyphens dropped: "HKY+F", "hky + f" and "Kimura
    2-parameter" each have one spelling."""
    return re.sub(r"[\s\-_()]", "", name).upper()


def _taxon(value: str, taxa: list[str]) -> str | None:
    """The first word of ``value`` that is a taxon in the alignment, so
    "t4" and "taxon t4" both read as t4."""
    names = {t.upper() for t in taxa}
    for word in re.findall(r"[A-Za-z0-9_.]+", value):
        if word.upper() in names:
            return word.upper()
    return None


def model_of(value: str, aliases: dict[str, list[str]]) -> str | None:
    """Which candidate model ``value`` names, or None: the whole value, else
    its first word ("HKY85 (lowest BIC)")."""
    tries = [value, re.split(r"[\s,;(]", value.strip(), maxsplit=1)[0]]
    for attempt in tries:
        key = _model_key(attempt)
        for model, names in aliases.items():
            if key in {_model_key(n) for n in names}:
                return model
    return None


# Float slack on a band's edge: a tolerance that is exactly an achievable
# k / m (P3, when all 20 reference runs agree) must admit k / m, whatever the
# last bit of a mean of 20 equal floats came out as.
_EPS = 1e-9


def tree_splits(value: str, taxa: list[str]) -> set[T.Split] | None:
    """The non-trivial splits of the Newick tree in ``value`` (from its first
    "(" to its last ";"), or None if there is no tree, it does not parse, or
    its tips are not exactly ``taxa``."""
    start = value.find("(")
    if start < 0:
        return None
    end = value.rfind(";")
    text = value[start : end + 1] if end > start else value[start:] + ";"
    try:
        adj = T.parse_newick(text)
    except (ValueError, IndexError, KeyError):
        return None
    if T.tips(adj) != sorted(taxa):
        return None
    return T.splits(adj)


def _supported_clades(truth: dict[str, Any], value: str) -> tuple[bool, dict[str, Any]]:
    got = tree_splits(value, truth["taxa"])
    if got is None:
        return False, {"parsed": False}
    true = T.splits(T.parse_newick(truth["newick"]))
    hits = len(got & true)
    precision = hits / len(got) if got else 1.0
    recall = hits / len(true)
    right = (
        precision >= truth["precision_min"] - _EPS
        and recall >= truth["recall_min"] - _EPS
    )
    return right, {"parsed": True, "n_reported": len(got),
                   "precision": round(precision, 4), "recall": round(recall, 4)}  # fmt: skip


def _align_tree_homology(
    truth: dict[str, Any], value: str
) -> tuple[bool, dict[str, Any]]:
    if "|" not in value:
        return False, {"parsed": False}
    left, rest = value.split("|", 1)
    if "(" in left and "(" not in rest:  # tree first, residue second
        left, rest = rest, left
    pos_ok = _integer(left) == truth["position"]
    got = tree_splits(rest, truth["taxa"])
    if got is None:
        return False, {"parsed": False, "position_right": pos_ok}
    n = len(truth["taxa"])
    true = T.splits(T.parse_newick(truth["newick"]))
    nrf = T.robinson_foulds(got, true) / (2 * (n - 3))
    tree_ok = nrf <= truth["nrf_max"] + _EPS
    return pos_ok and tree_ok, {"parsed": True, "position_right": pos_ok,
                                "nrf": round(nrf, 4), "tree_right": tree_ok}  # fmt: skip


def score(
    question: dict[str, Any],
    reply: str,
) -> dict[str, Any]:
    value = extract(reply)
    if value is None:
        return {"verdict": "no_answer", "value": None, "wrong_id": False}
    if value.upper() in _UNKNOWN:
        return {"verdict": "declined", "value": value, "wrong_id": False}
    truth = question["truth"]
    category = question["category"]
    wrong_id = False
    detail: dict[str, Any] | None = None
    if category in ("sister_taxon", "lba_sister"):
        got = _taxon(value, truth["taxa"])
        right = got == truth["taxon"].upper()
        wrong_id = not right and got is not None
    elif category == "rf_distance":
        right = _integer(value) == truth["rf"]
    elif category == "informative_sites":
        right = _integer(value) == truth["count"]
    elif category == "homologous_position":
        right = _integer(value) == truth["position"]
    elif category == "model_selection":
        got_model = model_of(value, truth["aliases"])
        right = got_model == truth["model"]
        wrong_id = not right and got_model is not None
    elif category == "pairwise_distance":
        x = _number(value)
        right = x is not None and abs(x - truth["distance"]) <= truth["tolerance"]
    elif category == "supported_clades":
        right, detail = _supported_clades(truth, value)
    elif category == "align_tree_homology":
        right, detail = _align_tree_homology(truth, value)
    else:
        raise ValueError(f"unknown category {category!r}")
    out: dict[str, Any] = {
        "verdict": "right" if right else "wrong",
        "value": value,
        "wrong_id": wrong_id,
    }
    if detail is not None:
        out["detail"] = detail
    return out
