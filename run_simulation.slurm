#!/bin/bash
#SBATCH --job-name=mprage_grappa
#SBATCH --partition=imgpascalq
#SBATCH --qos=img
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=/share/ConnLS/scratch/%u/logs/%x_%j.out
#SBATCH --error=/share/ConnLS/scratch/%u/logs/%x_%j.out
# NOTE: scratch log dir must exist before sbatch is called.
# Submit with:  mkdir -p /share/ConnLS/scratch/$USER/logs && sbatch run_simulation.slurm
set -eo pipefail

# ── Environment ───────────────────────────────────────────────────────────────
# Slurm copies the script to a tmpdir before execution, so BASH_SOURCE[0]
# resolves to that temp copy.  SLURM_SUBMIT_DIR is the directory from which
# sbatch was called — always use that to find project files.
PROJ_DIR="$SLURM_SUBMIT_DIR"
source "$PROJ_DIR/env.sh"

echo "================================================"
echo "Job:        $SLURM_JOB_NAME  ($SLURM_JOB_ID)"
echo "Node:       $SLURMD_NODENAME"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Started:    $(date)"
echo "================================================"

# Verify GPU is visible
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# Note: BART bart-img/0.9.00 is CPU-only (no CUDA linkage).
# KomaMRI Julia simulation uses the P100 GPU via CUDA.jl.

# ── Step 1: Julia simulation ──────────────────────────────────────────────────
# -t $SLURM_CPUS_PER_TASK exposes all allocated CPUs as Julia threads.
# CUDA.functional() inside the script detects the GPU automatically.
echo ""
echo "=== Simulation ==="
julia -t "$SLURM_CPUS_PER_TASK" "$PROJ_DIR/simulate_mprage_grappa.jl"

# ── Step 2: BART reconstruction ───────────────────────────────────────────────
echo ""
echo "=== Reconstruction ==="
OUT_DIR="/share/ConnLS/scratch/$USER/sim_output" bash "$PROJ_DIR/reconstruct_bart.sh"

echo ""
echo "================================================"
echo "Finished:   $(date)"
echo "================================================"
