#!/usr/bin/env bash
# For initialising on the University of Nottingham Ada HPC service
# source this file: . env.sh

# Resolve project root relative to this script regardless of CWD
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

module load julia-uoneasy/1.10.4-linux-x86_64
module load matlab-uon/r2025a
module load bart-img/0.9.00

# Point Julia at the project environment so `using KomaMRI` etc. work
# without an explicit Pkg.activate() in every script.
export JULIA_PROJECT="$PROJ_DIR"

# Julia depot stays in ~/.julia (persists across jobs on Ada's home filesystem).
# If you hit home-directory quota, override with a group/scratch path, e.g.:
#   export JULIA_DEPOT_PATH="/gpfs01/scratch/$USER/julia_depot:$JULIA_DEPOT_PATH"

echo "Loaded Julia, Matlab, and BART"
echo "JULIA_PROJECT=$JULIA_PROJECT"
echo ""
echo "First-time setup: run  julia setup_julia.jl  to install Julia packages."
