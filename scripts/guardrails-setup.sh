#!/usr/bin/env bash
# System prerequisites the Python environment cannot supply. Called by ci.yml
# and, when present and executable, by the template's guardrails workflows.
#
# MAFFT is the engine behind align_sequences. IQ-TREE arrives as a wheel
# (piqtree); MAFFT has no wheel, and the suite runs the real binary rather than
# a stand-in, so a runner without it must fail here, not skip there.
set -euo pipefail

if ! command -v mafft >/dev/null 2>&1; then
	sudo apt-get update -qq
	sudo apt-get install -y -qq mafft
fi
mafft --version 2>&1 | head -1
