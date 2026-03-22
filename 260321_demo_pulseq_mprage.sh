#!/usr/bin/env bash
# 260321_demo_pulseq_mprage.sh
# Generate MPRAGE .seq files using patched copies of the Pulseq demo scripts.
# (Originals in extern/pulseq/matlab/demoSeq/ call seq.plot() which hangs
# on headless HPC nodes.)
#
# Outputs:
#   mprage.seq        - fully sampled 3D MPRAGE
#   mprage_grappa.seq - GRAPPA R=2, 32 ACS lines
#
# Usage (after sourcing env.sh):
#   bash 260321_demo_pulseq_mprage.sh

set -euo pipefail
export QT_QPA_PLATFORM=offscreen

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Fully-sampled MPRAGE ==="
matlab -batch "addpath('${PROJ_DIR}','${PROJ_DIR}/extern/pulseq/matlab'); write_mprage_260321"
echo ""

echo "=== GRAPPA-accelerated MPRAGE ==="
matlab -batch "addpath('${PROJ_DIR}','${PROJ_DIR}/extern/pulseq/matlab'); write_mprage_grappa_260321"
echo ""

echo "--- Generated files ---"
ls -lh mprage.seq mprage_grappa.seq 2>/dev/null || echo "(no .seq files found)"
