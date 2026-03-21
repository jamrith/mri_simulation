#!/usr/bin/env bash
# 260321_demo_pulseq_mprage.sh
# Generate MPRAGE .seq files using the Pulseq demo scripts.
#
# Outputs (in current directory):
#   mprage.seq        — fully sampled 3D MPRAGE  (writeMPRAGE.m)
#   mprage_grappa.seq — GRAPPA R=2, 32 ACS lines (writeMPRAGE_grappa.m)
#
# The demo scripts live in extern/pulseq/matlab/demoSeq/.
# Each is a standalone Matlab script that constructs the sequence, runs a
# timing check, and writes a .seq file.  They end with `return` (which
# exits the batch session), so each must be its own matlab invocation.
#
# Usage (after sourcing env.sh):
#   bash 260321_demo_pulseq_mprage.sh

set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MPATH="addpath('${PROJ_DIR}/extern/pulseq/matlab', '${PROJ_DIR}/extern/pulseq/matlab/demoSeq');"

# Allow Matlab's Qt plotting without a display (headless HPC nodes)
export QT_QPA_PLATFORM=offscreen

echo "=== Fully-sampled MPRAGE (writeMPRAGE.m) ==="
matlab -batch "${MPATH} writeMPRAGE"
echo ""

echo "=== GRAPPA-accelerated MPRAGE (writeMPRAGE_grappa.m) ==="
matlab -batch "${MPATH} writeMPRAGE_grappa"
echo ""

echo "--- Generated files ---"
ls -lh mprage.seq mprage_grappa.seq 2>/dev/null || echo "(no .seq files found — check Matlab output above)"
