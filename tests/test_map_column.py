"""map_column called directly: the whole result, planted by hand.

Rows are letters, one per column, so every residue, position and window slice
below is read off the row, not recomputed by the code under test. `b` has gaps
at columns 4 and 11.
"""

from __future__ import annotations

import pytest

from phylokit_mcp.alignment import AlignmentError
from phylokit_mcp.positions import map_column

A = "abcdefghijklmnopqrst"
B = "ABC-EFGHIJ-LMNOPQRST"
ALN = {"a": A, "b": B}


def test_a_column_returns_every_row_and_a_window_of_five_either_side():
    assert map_column(ALN, None, None, 10) == {
        "column": 10,
        "n_columns": 20,
        "query": {"sequence": None, "position": None, "column": 10},
        "rows": {
            "a": {"position": 10, "residue": "j"},
            "b": {"position": 9, "residue": "J"},
        },
        "window": {
            "columns": [5, 15],
            "rows": {"a": "efghijklmno", "b": "EFGHIJ-LMNO"},
        },
    }


def test_a_position_is_counted_in_residues_of_its_own_row():
    # b's 10th residue is L, in column 12: two gaps come before it.
    assert map_column(ALN, "b", 10, None) == {
        "column": 12,
        "n_columns": 20,
        "query": {"sequence": "b", "position": 10, "column": None},
        "rows": {
            "a": {"position": 12, "residue": "l"},
            "b": {"position": 10, "residue": "L"},
        },
        "window": {
            "columns": [7, 17],
            "rows": {"a": "ghijklmnopq", "b": "GHIJ-LMNOPQ"},
        },
    }


def test_the_first_and_last_columns_are_valid_and_the_window_stops_at_the_ends():
    first = map_column(ALN, None, None, 1)
    assert first["rows"]["a"] == {"position": 1, "residue": "a"}
    assert first["window"] == {
        "columns": [1, 6],
        "rows": {"a": "abcdef", "b": "ABC-EF"},
    }
    last = map_column(ALN, None, None, 20)
    assert last["rows"]["b"] == {"position": 18, "residue": "T"}
    assert last["window"] == {
        "columns": [15, 20],
        "rows": {"a": "opqrst", "b": "OPQRST"},
    }
    for bad in (0, 21):
        with pytest.raises(AlignmentError) as exc:
            map_column(ALN, None, None, bad)
        assert str(exc.value) == (
            f"column must be 1-20 (the alignment has 20 columns); got {bad}."
        )


def test_the_query_must_name_one_place_and_its_row():
    assert map_column(ALN, "a", 1, None)["column"] == 1
    for position, column in ((None, None), (1, 1)):
        with pytest.raises(AlignmentError) as exc:
            map_column(ALN, "a", position, column)
        assert str(exc.value) == (
            "Give exactly one of position (a residue of `sequence`) and column "
            "(an alignment column)."
        )
    with pytest.raises(AlignmentError) as exc:
        map_column(ALN, None, 1, None)
    assert str(exc.value) == "position needs sequence: the row the residue is in."


def test_a_position_must_be_one_of_the_rows_own_residues():
    # b has 18 residues in 20 columns; its last residue, T, is in column 20.
    assert map_column(ALN, "b", 1, None)["column"] == 1
    assert map_column(ALN, "b", 18, None)["column"] == 20
    for bad in (0, 19, "3"):
        with pytest.raises(AlignmentError) as exc:
            map_column(ALN, "b", bad, None)
        assert str(exc.value) == (
            "position must be 1-18: 'b' has 18 residues (gaps not counted); "
            f"got {bad!r}."
        )


def _named(n: int) -> dict[str, str]:
    return {f"s{i:02d}": "AC" for i in range(n)}


def test_an_unknown_row_for_a_position_lists_twenty_names_then_an_ellipsis():
    twenty = [f"s{i:02d}" for i in range(20)]
    assert map_column(_named(20), "s19", 2, None)["column"] == 2
    with pytest.raises(AlignmentError) as exc:
        map_column(_named(20), "zz", 1, None)
    assert str(exc.value) == f"No sequence named 'zz'. Names: {twenty}"
    with pytest.raises(AlignmentError) as exc:
        map_column(_named(21), "zz", 1, None)
    assert str(exc.value) == f"No sequence named 'zz'. Names: {twenty} ..."


def test_a_named_row_is_checked_for_a_column_too_with_the_same_ellipsis():
    twenty = [f"s{i:02d}" for i in range(20)]
    assert map_column(_named(21), "s20", None, 2)["query"]["sequence"] == "s20"
    with pytest.raises(AlignmentError) as exc:
        map_column(_named(20), "zz", None, 2)
    assert str(exc.value) == f"No sequence named 'zz'. Names: {twenty}"
    # Given in reverse order, the names are still listed sorted.
    backwards = dict(reversed(_named(21).items()))
    with pytest.raises(AlignmentError) as exc:
        map_column(backwards, "zz", None, 2)
    assert str(exc.value) == f"No sequence named 'zz'. Names: {twenty} ..."
