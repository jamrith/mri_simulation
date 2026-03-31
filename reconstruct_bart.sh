#!/usr/bin/env bash
# reconstruct_bart.sh
#
# BART reconstruction of MPRAGE-GRAPPA k-space exported by
# simulate_mprage_grappa.jl.
#
# Two reconstructions are produced:
#   mag_zerofill  — inverse FFT on zero-filled k-space (shows GRAPPA aliasing)
#   mag_pics      — compressed-sensing (wavelet-regularised) via bart pics
#                   (suppresses aliasing without needing multi-coil data)
#
# Usage:
#   bash reconstruct_bart.sh            # uses sim_output/ in current dir
#   OUT_DIR=/path/to/data bash reconstruct_bart.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/sim_output}"
RECO_DIR="$OUT_DIR/recon"

# Load reconstruction parameters written by the Julia simulation
PARAMS="$OUT_DIR/recon_params.sh"
if [[ ! -f "$PARAMS" ]]; then
    echo "ERROR: $PARAMS not found — run simulate_mprage_grappa.jl first." >&2
    exit 1
fi
source "$PARAMS"

mkdir -p "$RECO_DIR"
cd "$RECO_DIR"

echo "=== BART MPRAGE-GRAPPA Reconstruction ==="
echo "  K-space:    ${NRO} (RO) × ${NPE} (PE) × ${NPAR} (PAR)"
echo "  ACS region: PE lines ${ACS_START} … ${ACS_END}  (${NACS} lines)"
echo "  Output dir: $RECO_DIR"
echo ""

# ── Step 1: FFT modulation ─────────────────────────────────────────────────
# BART's fft convention: DC at index [0,0,0].
# kSpaceCenterLine sits at the array centre, so fftmod applies the
# alternating +/-1 pattern that is equivalent to an fftshift.
echo "Step 1: fftmod (centre DC at corner for BART FFT)..."
bart fftmod 7 "$OUT_DIR/kspace" kspace_mod

# ── Step 2: Zero-filled reconstruction (baseline) ─────────────────────────
echo "Step 2: Zero-filled inverse FFT..."
bart fft -i 7 kspace_mod image_zerofill_cplx
bart rss 8 image_zerofill_cplx mag_zerofill
echo "  -> mag_zerofill.cfl"

# ── Step 3: Compressed-sensing reconstruction via bart pics ───────────────
# With a single receive coil (from KomaMRI), coil-based GRAPPA is degenerate.
# pics with wavelet regularisation is the appropriate single-coil approach:
# it uses the ACS-region signal level to suppress the aliasing that arises
# from the R=2 undersampling.
#
# Flags:
#   -i 100          up to 100 CG iterations
#   -R W:7:0:0.001  wavelet (W) over dims 0+1+2 (7=111b), no joint sparsity,
#                   regularisation weight 0.001
#   -e              enforce data consistency at acquired samples
echo "Step 3: Compressed-sensing reconstruction (bart pics + wavelet)..."

# Sensitivity map: all-ones for single coil (dims: RO × PE × PAR × 1 coil)
bart ones 4 "$NRO" "$NPE" "$NPAR" 1 sens

bart pics \
    -i 100 \
    -R W:7:0:0.001 \
    -e \
    kspace_mod sens image_pics_cplx

bart rss 8 image_pics_cplx mag_pics
echo "  -> mag_pics.cfl"

echo ""
echo "=== Done ==="
echo "  mag_zerofill.cfl  — zero-filled FFT (aliased)"
echo "  mag_pics.cfl      — wavelet CS reconstruction (alias-suppressed)"
echo ""
echo "View with:  bart toimg mag_pics mag_pics.png"
echo "        or: python3 -c \""
echo "    import cfl, matplotlib.pyplot as plt, numpy as np"
echo "    img = np.abs(cfl.readcfl('$RECO_DIR/mag_pics'))"
echo "    plt.imshow(img[:,:,img.shape[2]//2], cmap='gray'); plt.show()\""
