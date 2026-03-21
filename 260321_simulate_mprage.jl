#!/usr/bin/env julia
# 260321_simulate_mprage.jl
# Test the MRI simulation pipeline: Pulseq MPRAGE -> KomaMRI -> reconstruction.
#
# Usage (after sourcing env.sh):
#   julia 260321_simulate_mprage.jl                    # defaults
#   julia 260321_simulate_mprage.jl mprage_hcp.seq 8   # custom seq + subsampling

using KomaMRI
using MAT       # for saving results as .mat

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
seq_file   = length(ARGS) >= 1 ? ARGS[1] : "mprage_hcp.seq"
phantom_ss = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8  # subsampling factor

# Higher ss = coarser phantom = faster simulation.
# ss=8  -> ~4 mm spacing, few thousand spins  (minutes on CPU)
# ss=4  -> ~2 mm spacing, tens of thousands   (hours on CPU, minutes on GPU)
# ss=2  -> ~1 mm spacing, hundreds of thousands (GPU recommended)

# ---------------------------------------------------------------------------
# Load sequence
# ---------------------------------------------------------------------------
println("Loading sequence: $seq_file")
seq = read_seq(seq_file)
println("  Sequence duration: $(round(dur(seq); digits=2)) s")

# ---------------------------------------------------------------------------
# Create phantom
# ---------------------------------------------------------------------------
println("Creating 3D brain phantom (ss=$phantom_ss)...")
obj = brain_phantom3D(; ss=phantom_ss)
println("  Phantom spins: $(length(obj))")

# ---------------------------------------------------------------------------
# Scanner (3 T, matching HCP Connectom)
# ---------------------------------------------------------------------------
sys = Scanner(B0=3.0)

# ---------------------------------------------------------------------------
# Simulate
# ---------------------------------------------------------------------------
println("Starting simulation...")
sim_params = Dict{String,Any}(
    "Nblocks" => 20,            # GPU batch size (increase for more VRAM)
    "gpu"     => false,         # set true if CUDA GPU is available
    "return_type" => "raw",     # RawAcquisitionData (ISMRMRD)
)
raw = simulate(obj, seq, sys; sim_params)
println("Simulation complete.")

# ---------------------------------------------------------------------------
# Reconstruction (simple IFFT via MRIReco)
# ---------------------------------------------------------------------------
println("Reconstructing...")
# MRIReco's reconstruction handles 3D Cartesian data from ISMRMRD raw.
# For GRAPPA-accelerated data, use BART instead (see notes below).
acq = AcquisitionData(raw)
Nx, Ny, Nz = acq.encodingSize
rec_params = Dict{Symbol,Any}(
    :reconSize => (Nx, Ny, Nz),
    :densityWeighting => true,
)
image = reconstruction(acq, rec_params)
println("  Image size: $(size(image))")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
outfile = replace(seq_file, ".seq" => "_sim.mat")
matwrite(outfile, Dict(
    "image"      => abs.(image[:,:,:,1]),
    "seq_file"   => seq_file,
    "phantom_ss" => phantom_ss,
))
println("Saved: $outfile")

# ---------------------------------------------------------------------------
# BART reconstruction (for GRAPPA data)
# ---------------------------------------------------------------------------
# To reconstruct GRAPPA-accelerated data with BART:
#
# 1. Save k-space in CFL format:
#      ksp = raw_to_kspace(raw)  # reshape to [Nro, Npe1, Npe2, Ncoil]
#      write_cfl("kspace", ksp)
#
# 2. Run BART (after `module load bart-img/0.9.00`):
#      bart ecalib -m1 kspace sens         # estimate coil sensitivities
#      bart pics -l1 -r0.01 kspace sens output  # L1-regularised recon
#
# CFL helper (complex float32, column-major):
#   function write_cfl(name, data)
#       open(name * ".hdr", "w") do f
#           println(f, "# Dimensions")
#           println(f, join(size(data), " "))
#       end
#       open(name * ".cfl", "w") do f
#           for v in vec(data)
#               write(f, Float32(real(v)))
#               write(f, Float32(imag(v)))
#           end
#       end
#   end
