# mri_simulation
Using KomaMRI and Pulseq to simulate MRI acquisitions under various conditions.
Reconstruction is handled via BART.

## Setup (University of Nottingham Ada HPC)

### Clone

```bash
git clone --recurse-submodules <repo-url>
cd mri_simulation
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

### Load modules and install Julia packages

```bash
. env.sh
julia setup_julia.jl
```

`env.sh` loads Julia 1.10.4, Matlab, and BART, and sets `JULIA_PROJECT` so all
Julia scripts use the project environment automatically.

`setup_julia.jl` only needs to be run once per user. It dev-links KomaMRI and
its subpackages from the local `extern/KomaMRI.jl` submodule. Julia packages
are cached in `~/.julia` and persist across jobs.

### Subsequent sessions

```bash
. env.sh
julia my_script.jl
```

## Running the MPRAGE-GRAPPA simulation

```bash
# Submit from the project root — logs/ must exist before sbatch is called
mkdir -p logs
sbatch run_simulation.slurm
```

Outputs land in `sim_output/` (k-space CFL files) and `sim_output/recon/` (reconstructed images).

To run interactively instead:

```bash
. env.sh
julia -t 4 simulate_mprage_grappa.jl
bash reconstruct_bart.sh
```

## Dependencies

| Tool | Source |
|------|--------|
| [KomaMRI.jl](https://github.com/JuliaHealth/KomaMRI.jl) | `extern/KomaMRI.jl` (submodule) |
| [pulseq](https://github.com/pulseq/pulseq) | `extern/pulseq` (submodule) |
| CUDA.jl | Julia package (GPU backend for KomaMRI) |
| BART | `module load bart-img/0.9.00` — CPU-only build |
| CUDA | `module load cuda-uoneasy/12.6.0` (Ada module) |
| Matlab | `module load matlab-uon/r2025a` (Ada module) |
