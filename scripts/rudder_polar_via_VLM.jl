#!/usr/bin/env julia
#
# Rudder polar via LiftingSurfaces.jl VLM. Sweeps rudder angle δ from
# -25° to +25° and tabulates CL, CD, side-force, and the slope dCL/dα
# in the linear range. No WaterLily flow — uniform inflow only — but
# uses the exact LiftingSurfaces.Rudder + rudder_forces path that the
# coupled scripts use.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)
println("packages added"); flush(stdout)

using LiftingSurfaces
using Printf

const Vinf  = 1.0
const chord = 1.0
const span  = 2.0
const ns    = 16
const nc    = 8

rudder = Rudder(; chord=chord, span=span, ns=ns, nc=nc)
@printf "=== Rudder polar via LiftingSurfaces.jl ===\n"
@printf "  chord=%.1f, span=%.1f (AR=%.1f), ns=%d, nc=%d\n\n" chord span span/chord ns nc

δ_list = -25.0:5.0:25.0
@printf "  δ (°)   CL         CD         CY         L/D\n"
@printf "  ────   ────────   ────────   ────────   ──────\n"
results = []
for δ_deg in δ_list
    r = rudder_forces(rudder, deg2rad(δ_deg), Vinf)
    LoD = abs(r.CD) > 1e-9 ? r.CL / r.CD : Inf
    @printf "  %+5.1f  %+8.4f   %+8.5f  %+8.5f  %+6.2f\n" δ_deg r.CL r.CD r.CY LoD
    push!(results, (δ_deg, r.CL, r.CD, r.CY))
end

# Linear-range slope: CL/α at ±5°
CL_p5  = results[findfirst(r -> r[1] ==  5.0, results)][2]
CL_m5  = results[findfirst(r -> r[1] == -5.0, results)][2]
slope  = (CL_p5 - CL_m5) / deg2rad(10.0)
elliptic = 2π * (span/chord) / (span/chord + 2)
@printf "\n  dCL/dα (linear range, ±5°)        = %.3f per rad\n" slope
@printf "  Prandtl-LLT (elliptic loading)     = %.3f per rad\n" elliptic
@printf "  Ratio VLM/LLT                      = %.3f  (rectangular ≈ 0.78)\n" slope/elliptic

# Sanity: drag minima at δ=0
@printf "\n  CD at δ=0                          = %.6f  (should be ~0)\n" results[findfirst(r -> r[1]==0.0, results)][3]
@printf "  CD at δ=±25                         = %.5f, %.5f  (induced drag ∝ CL²)\n" results[1][3] results[end][3]
