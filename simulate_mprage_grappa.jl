#!/usr/bin/env julia
# simulate_mprage_grappa.jl
#
# Simulate mprage_grappa.seq on a 3D brain phantom using KomaMRI,
# sort the acquired k-space lines, and write BART CFL files ready
# for GRAPPA reconstruction.
#
# Usage:
#   julia -t <nthreads> simulate_mprage_grappa.jl [--gpu]
#
# Outputs (in sim_output/):
#   kspace.hdr / .cfl   — zero-filled 3-D k-space  (Nro × Npe × Npar × 1)
#   kspace_acs.hdr/.cfl — ACS calibration region   (Nro × Nacs × Npar × 1)
#   mask.hdr / .cfl     — binary sampling mask      (1 × Npe × 1 × 1)
#   recon_params.sh     — shell vars for reconstruct_bart.sh

using KomaMRI
using CUDA   # triggers KomaCUDAExt in KomaMRICore — GPU used automatically when available

# ─── Configuration ────────────────────────────────────────────────────────────

const SEQ_FILE = joinpath(@__DIR__, "mprage_grappa.seq")
const OUT_DIR  = get(ENV, "SIM_OUT_DIR", "/share/ConnLS/scratch/lpxjr4/sim_output")

# Phantom subsampling along each axis.
# ss=4 → 2 mm voxels, ~130 k spins (recommended for HPC CPU runs)
# ss=6 → 3 mm voxels, ~40 k spins  (quick test)
# ss=2 → 1 mm voxels, ~1 M spins   (high quality, needs GPU)
const PHANTOM_SS = parse(Int, get(ENV, "PHANTOM_SS", "4"))

# GPU: use if CUDA is functional (loaded via cuda-uoneasy/12.6.0 + CUDA.jl)
const USE_GPU = CUDA.functional()

# ─── BART CFL helpers ─────────────────────────────────────────────────────────

"""
Write a complex array to BART's CFL format (<basename>.hdr + <basename>.cfl).

BART stores data in column-major (Fortran / Julia) order.  Dimensions are
padded to at least 4 (READ × PHASE1 × PHASE2 × COILS) per BART convention.
"""
function write_bart_cfl(basename::String, data::AbstractArray{<:Complex})
    dims = collect(size(data))
    while length(dims) < 4
        push!(dims, 1)
    end
    open(basename * ".hdr", "w") do f
        println(f, "# Dimensions")
        println(f, join(dims, " "))
    end
    open(basename * ".cfl", "w") do f
        for x in data
            write(f, Float32(real(x)))
            write(f, Float32(imag(x)))
        end
    end
    println("  Wrote $(basename)  [$(join(dims, "×"))]")
end

# ─── Main ─────────────────────────────────────────────────────────────────────

mkpath(OUT_DIR)

# ── 1. Load sequence ──────────────────────────────────────────────────────────
println("\n=== Step 1: Load sequence ===")
seq = read_seq(SEQ_FILE)
println("  Sequence:       $(get(seq.DEF, "Name", "unknown"))")
println("  Total duration: $(round(dur(seq); digits=1)) s")
fov_m = get(seq.DEF, "FOV", [0.0, 0.0, 0.0])
println("  FOV:            $(round.(fov_m .* 1e3; digits=1)) mm")

k_center_lin = Int(get(seq.DEF, "kSpaceCenterLine", 0.0))
println("  kSpaceCenterLine: $k_center_lin")

# ── 2. Build 3D brain phantom ─────────────────────────────────────────────────
println("\n=== Step 2: Build 3D brain phantom (ss=$PHANTOM_SS) ===")
# brain_phantom3D uses the BrainWeb digital phantom.
# T1/T2/T2s/ρ are set to realistic 3 T values by default.
# start_end selects the axial slab in z (default [160,200] ≈ 1 cm slab).
# Extend the slab to cover the full sequence FOV in z:
# Original volume is 181 slices at 0.5 mm = 90.5 mm; FOV_z = 256 mm here,
# so we take the central 80 slices which at ss=4 (step 4) gives ~10 mm slab.
# For a fuller simulation increase start_end range.
@time obj = brain_phantom3D(; ss=PHANTOM_SS, start_end=[140, 220])
println("  $(length(obj)) spins  ($(round(length(obj)/1e3; digits=1)) k)")

# ── 3. Pre-compute sequence labels ────────────────────────────────────────────
# get_label is fast (O(n_blocks)) and gives us LIN/PAR for every block.
# We do this before simulate() so we can use return_type="mat", which avoids
# the internal signal_to_raw_data() call that recomputes the full k-space
# trajectory at all 13 M ADC sample times (~50 min hang).
println("\n=== Step 3: Pre-compute sequence labels ===")
@time labels = get_label(seq)
println("  $(length(labels)) blocks labelled")

# ── 4. Scanner ────────────────────────────────────────────────────────────────
sys = Scanner()   # Default 1.5 T; MPRAGE sequence was designed for 3 T but
                  # the Bloch simulation is field-independent beyond γ scaling.

# ── 5. Simulate ───────────────────────────────────────────────────────────────
println("\n=== Step 4: Simulate ===")
println("  GPU: $USE_GPU  |  Threads: $(Threads.nthreads())")
println("  NOTE: Full 3D MPRAGE is ~340 s of sequence time.")
println("        Expect wall-clock minutes–hours depending on hardware.")

sim_params = Dict{String,Any}(
    "return_type" => "mat",     # Return raw signal matrix; avoids slow signal_to_raw_data
    "gpu"         => USE_GPU,
    "Nthreads"    => Threads.nthreads(),
    "Nblocks"     => 40,        # Smaller blocks → less peak RAM
    "precision"   => "f32",     # float32 sufficient for MRI SNR levels
)

@time sig = simulate(obj, seq, sys; sim_params)
# sig: (total_ADC_samples, 1) ComplexF32
println("  Total ADC samples: $(size(sig, 1))")

# ── 6. Sort signal into 3D k-space array ──────────────────────────────────────
println("\n=== Step 5: Sort k-space ===")

# Nro: samples per readout — use the first ADC block
first_adc_b = findfirst(b -> is_ADC_on(seq[b]), 1:length(seq))
Nro = seq.ADC[first_adc_b].N

# k-space dimensions from label ranges at ADC blocks
lin_at_adc = [labels[b].LIN for b in 1:length(seq) if is_ADC_on(seq[b])]
par_at_adc = [labels[b].PAR for b in 1:length(seq) if is_ADC_on(seq[b])]
lin_max, par_max = maximum(lin_at_adc), maximum(par_at_adc)
Npe  = lin_max + 1
Npar = par_max + 1

println("  Readout points:  $Nro  (includes 2× oversampling if sequence uses it)")
println("  PE lines full:   $Npe  (LIN 0 … $lin_max)")
println("  Partitions:      $Npar (PAR 0 … $par_max)")

# Zero-filled k-space: shape (Nro, Npe, Npar, 1)  — 1 coil
kspace = zeros(ComplexF32, Nro, Npe, Npar, 1)
# Sampling mask: shape (1, Npe, 1, 1) — just PE dimension (partitions fully sampled)
mask = zeros(ComplexF32, 1, Npe, 1, 1)

let sample_offset = 1
    for b in 1:length(seq)
        if is_ADC_on(seq[b])
            N   = seq.ADC[b].N
            lin = labels[b].LIN + 1   # 1-indexed
            par = labels[b].PAR + 1
            if 1 <= lin <= Npe && 1 <= par <= Npar
                kspace[:, lin, par, 1] .= sig[sample_offset:sample_offset+N-1, 1]
                mask[1, lin, 1, 1] = 1f0 + 0im
            end
            sample_offset += N
        end
    end
end

n_acq  = Int(real(sum(mask))) * Npar   # acquired PE × PAR points
n_full = Npe * Npar
pct    = round(100.0 * n_acq / n_full; digits=1)
println("  Acquired $n_acq / $n_full PE×PAR combinations ($pct%)")

# ── 7. Extract ACS calibration region ─────────────────────────────────────────
# The GRAPPA sequence fully samples 32 ACS lines centred on k-space DC.
# kSpaceCenterLine is the LIN index of the DC component (0-indexed).
n_acs = 32
acs_half = n_acs ÷ 2
# 0-indexed range, clamped to array bounds
acs_start_0 = max(0, k_center_lin - acs_half)
acs_end_0   = min(Npe - 1, acs_start_0 + n_acs - 1)
n_acs_actual = acs_end_0 - acs_start_0 + 1

# 1-indexed slice for Julia array
acs_idx = (acs_start_0 + 1):(acs_end_0 + 1)
kspace_acs = kspace[:, acs_idx, :, :]   # (Nro, Nacs, Npar, 1)

println("  ACS region: LIN $acs_start_0 … $acs_end_0 ($n_acs_actual lines)")

# ── 8. Write BART CFL files ───────────────────────────────────────────────────
println("\n=== Step 6: Write BART files ===")

write_bart_cfl(joinpath(OUT_DIR, "kspace"),     kspace)
write_bart_cfl(joinpath(OUT_DIR, "kspace_acs"), kspace_acs)
write_bart_cfl(joinpath(OUT_DIR, "mask"),       mask)

# Reconstruction parameters for the shell script
open(joinpath(OUT_DIR, "recon_params.sh"), "w") do f
    println(f, "# Auto-generated by simulate_mprage_grappa.jl — do not edit")
    println(f, "NRO=$Nro")
    println(f, "NPE=$Npe")
    println(f, "NPAR=$Npar")
    println(f, "K_CENTER_LIN=$k_center_lin")
    println(f, "ACS_START=$acs_start_0")
    println(f, "ACS_END=$acs_end_0")
    println(f, "NACS=$n_acs_actual")
end
println("  Wrote recon_params.sh")

println("\n=== Done ===")
println("Output directory: $OUT_DIR")
println("Run the reconstruction:")
println("  bash reconstruct_bart.sh")
