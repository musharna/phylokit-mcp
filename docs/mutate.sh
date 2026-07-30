#!/usr/bin/env bash
# Mutation gate. Asserts GREEN BEFORE mutating — without that, a suite that was
# already red reads as "mutant killed" for every mutant.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
SRC=src/phylokit_mcp

green_check() {
	$PY -m pytest "$1" -q >/dev/null 2>&1
}

run_mutant() {
	local name="$1" file="$2" from="$3" to="$4" tests="$5"
	cp "$file" "$file.bak"
	if ! grep -qF "$from" "$file"; then
		echo "SKIP  $name -- anchor not found, mutation never applied"
		rm -f "$file.bak"
		return
	fi
	# shellcheck disable=SC2016
	$PY - "$file" "$from" "$to" <<'PYEOF'
import sys
p, a, b = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(p).read()
assert a in s, "anchor vanished"
open(p, "w").write(s.replace(a, b, 1))
PYEOF
	if $PY -m pytest $tests -q >/dev/null 2>&1; then
		echo "SURVIVED  $name   <-- NOT COVERED"
	else
		echo "killed    $name"
	fi
	mv "$file.bak" "$file"
}

echo "=== asserting green before mutating ==="
if ! green_check tests/; then
	echo "ABORT: suite is not green before mutation; every result below would be meaningless."
	exit 1
fi
echo "green."
echo
echo "=== mutants ==="

run_mutant "1 support always 1.0" \
	"$SRC/bootstrap.py" \
	"counts[split] = counts.get(split, 0) + 1" \
	"counts[split] = replicates" \
	"tests/test_support.py"

run_mutant "2 splits not canonicalised (no reference anchor)" \
	"$SRC/splits.py" \
	"side = below if ref not in below else tips - below" \
	"side = below" \
	"tests/test_splits.py tests/test_support.py"

run_mutant "3 bootstrap resamples WITHOUT replacement" \
	"$SRC/bootstrap.py" \
	"idx = rng.integers(0, n_sites, size=n_sites)" \
	"idx = rng.permutation(n_sites)" \
	"tests/test_support.py"

run_mutant "4 trivial splits kept (terminal edges counted)" \
	"$SRC/splits.py" \
	"if 2 <= len(side) <= len(tips) - 2:" \
	"if len(side) >= 1:" \
	"tests/test_splits.py"

run_mutant "5 parsimony-informative counts singletons" \
	"$SRC/alignment.py" \
	"if sum(1 for v in tally.values() if v >= 2) >= 2:" \
	"if len(tally) >= 2:" \
	"tests/test_alignment.py"

run_mutant "6 ragged alignment accepted" \
	"$SRC/alignment.py" \
	"if len(lengths) != 1:" \
	"if False:" \
	"tests/test_alignment.py"

run_mutant "7 AIC ignores the parameter penalty" \
	"$SRC/inference.py" \
	"aic = 2 * k - 2 * lnl" \
	"aic = -2 * lnl" \
	"tests/test_models.py"

run_mutant "8 model ranking returns winner only (no runners-up)" \
	"$SRC/inference.py" \
	"ties = [m for m in ranked[1:] if m.delta <= DELTA_INDISTINGUISHABLE]" \
	"ties = []" \
	"tests/test_models.py"

run_mutant "9 unsupported-clade advisory never fires" \
	"$SRC/diagnostics.py" \
	"weak = [c.as_dict() for c in result.clades if c.support < WEAK_SUPPORT]" \
	"weak = []" \
	"tests/test_support.py"

run_mutant "10 server claims bit-exact reproducibility" \
	"$SRC/engine.py" \
	'"bit_exact_on_repeat_within_process": False,' \
	'"bit_exact_on_repeat_within_process": True,' \
	"tests/test_support.py"

echo
echo "=== confirming clean restore ==="
green_check tests/ && echo "green after restore." || echo "NOT GREEN AFTER RESTORE -- files may be dirty"
