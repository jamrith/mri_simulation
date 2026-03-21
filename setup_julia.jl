#!/usr/bin/env julia
# Run once to install Julia dependencies into the project environment.
# Requires Julia 1.10.4 (julia-uoneasy/1.10.4-linux-x86_64 on Ada).
# julia-img/latest (1.12.3) segfaults during Pkg operations on this system.
#
# Usage (after sourcing env.sh):
#   julia setup_julia.jl
#
# This devs KomaMRI from the local submodule so the source is pinned to
# whatever commit the submodule tracks rather than the registry version.

using Pkg

println("Activating project environment...")
Pkg.activate(@__DIR__)

# Dev all workspace subpackages from local source first so the resolver
# doesn't try to pull them from the registry (which may have older versions).
koma_root = joinpath(@__DIR__, "extern", "KomaMRI.jl")
for sub in ["KomaMRIBase", "KomaMRICore", "KomaMRIFiles", "KomaMRIPlots"]
    println("Dev-ing $sub...")
    Pkg.develop(path=joinpath(koma_root, sub))
end
println("Dev-ing KomaMRI (top-level)...")
Pkg.develop(path=koma_root)

println("\nDone. Verify with: using KomaMRI")
