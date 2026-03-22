#!/usr/bin/env bash
# 260321_generate_hcp_mprage.sh
# Generate HCP-protocol MPRAGE .seq files using mprage_hcp_260321.m.
# Both variants are produced in a single MATLAB call.
#
# Outputs:
#   mprage_hcp.seq        - fully sampled 3D MPRAGE (0.7 mm iso)
#   mprage_hcp_grappa.seq - GRAPPA R=2, 32 ACS lines
#
# Usage (after sourcing env.sh):
#   bash 260321_generate_hcp_mprage.sh

set -euo pipefail
export QT_QPA_PLATFORM=offscreen

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== HCP MPRAGE (fully-sampled + GRAPPA) ==="
matlab -batch "addpath('${PROJ_DIR}','${PROJ_DIR}/extern/pulseq/matlab'); mprage_hcp_260321"
echo ""

echo "--- Generated files ---"
ls -lh "${PROJ_DIR}/mprage_hcp.seq" "${PROJ_DIR}/mprage_hcp_grappa.seq" 2>/dev/null || echo "(no .seq files found)"
