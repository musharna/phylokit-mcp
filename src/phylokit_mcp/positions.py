"""Alignment columns and the residue positions they hold.

"Which residue of B is in the same column as residue 332 of A" is a counting
job along two rows of an alignment, and done by eye it slips: in this repo's
agent benchmark (10-16) every wrong answer the server arm gave was a residue
count off by one to ten after `align_sequences` had produced the right
alignment. The count is mechanical, so it is done here.

Positions are 1-based and count residues only; columns are 1-based. A gap is
'-'. '?' and '.' are refused in a mapped row, because whether they occupy a
residue position (missing data at a residue) or not (padding) changes every
count after them, and the file does not say which.
"""

from __future__ import annotations

from typing import Any

from .alignment import AlignmentError

GAP = "-"
_AMBIGUOUS = set("?.")


def _check_rows(aligned: dict[str, str]) -> int:
    lengths = {len(s) for s in aligned.values()}
    if len(lengths) != 1:
        raise AlignmentError(
            "Rows are not all the same length, so this is not an alignment. If "
            "the sequences are unaligned, pass aligned=false and they are "
            "aligned with MAFFT first."
        )
    (n,) = lengths
    if n == 0:
        raise AlignmentError("Rows are empty (zero columns).")
    for nm, row in aligned.items():
        bad = _AMBIGUOUS & set(row)
        if bad:
            raise AlignmentError(
                f"Row {nm!r} contains {sorted(bad)}. Whether those hold a residue "
                "position or are padding changes every position after them, so "
                "this tool will not guess: replace them with '-' (padding) or "
                "'N'/'X' (an unknown residue)."
            )
    return n


def _residue_index(row: str) -> list[int | None]:
    """For each column, the 1-based residue position it holds, or None at a gap."""
    out: list[int | None] = []
    pos = 0
    for ch in row:
        if ch == GAP:
            out.append(None)
        else:
            pos += 1
            out.append(pos)
    return out


def _at_column(row: str, index: list[int | None], col: int) -> dict[str, Any]:
    """What one row holds at a 0-based column: the residue and its position, or
    at a gap the positions of the nearest residues on either side."""
    pos = index[col]
    if pos is not None:
        return {"position": pos, "residue": row[col]}
    before = next(
        (index[i] for i in range(col - 1, -1, -1) if index[i] is not None), None
    )
    after = next(
        (index[i] for i in range(col + 1, len(row)) if index[i] is not None), None
    )
    return {"position": None, "residue": GAP, "between": [before, after]}


def map_column(
    aligned: dict[str, str],
    sequence: str | None,
    position: int | None,
    column: int | None,
) -> dict[str, Any]:
    """The column of `sequence`'s residue `position`, or the given `column`, and
    what every row holds there."""
    n_columns = _check_rows(aligned)
    if (position is None) == (column is None):
        raise AlignmentError(
            "Give exactly one of position (a residue of `sequence`) and column "
            "(an alignment column)."
        )
    if position is not None:
        if sequence is None:
            raise AlignmentError("position needs sequence: the row the residue is in.")
        if sequence not in aligned:
            raise AlignmentError(
                f"No sequence named {sequence!r}. Names: {sorted(aligned)[:20]}"
                + (" ..." if len(aligned) > 20 else "")
            )
        index = _residue_index(aligned[sequence])
        n_res = sum(p is not None for p in index)
        if not isinstance(position, int) or not 1 <= position <= n_res:
            raise AlignmentError(
                f"position must be 1-{n_res}: {sequence!r} has {n_res} residues "
                f"(gaps not counted); got {position!r}."
            )
        col = index.index(position)
    else:
        if sequence is not None and sequence not in aligned:
            raise AlignmentError(
                f"No sequence named {sequence!r}. Names: {sorted(aligned)[:20]}"
            )
        if not isinstance(column, int) or not 1 <= column <= n_columns:
            raise AlignmentError(
                f"column must be 1-{n_columns} (the alignment has {n_columns} "
                f"columns); got {column!r}."
            )
        col = column - 1
    rows = {
        nm: _at_column(row, _residue_index(row), col) for nm, row in aligned.items()
    }
    lo, hi = max(0, col - 5), min(n_columns, col + 6)
    return {
        "column": col + 1,
        "n_columns": n_columns,
        "query": {"sequence": sequence, "position": position, "column": column},
        "rows": rows,
        "window": {
            "columns": [lo + 1, hi],
            "rows": {nm: row[lo:hi] for nm, row in aligned.items()},
        },
    }
