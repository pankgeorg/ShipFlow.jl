#!/usr/bin/env julia
#
# Propeller open-water polar via LiftingSurfaces.BladedRotor. Sweeps
# the advance ratio J = V∞ / (n D) and tabulates KT, KQ, η. No
# WaterLily flow — uniform inflow only.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl"))]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)
println("packages added"); flush(stdout)

using LiftingSurfaces
using Printf

const Vinf  = 1.0
const R     = 1.0
const D     = 2 * R
const Rhub  = 0.2
const N_blades = 3

rotor = BladedRotor(; N_blades=N_blades, R=R, R_hub=Rhub,
    chord = (0.25, 0.18),
    twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 16, nc = 6)

@printf "=== Propeller polar via LiftingSurfaces.jl ===\n"
@printf "  %d blades, R=%.1f, D=%.1f, R_hub=%.1f, chord(%.2f→%.2f), twist(%d°→%d°)\n\n" N_blades R D Rhub rotor.chord[1] rotor.chord[2] round(Int, rad2deg(rotor.twist[1])) round(Int, rad2deg(rotor.twist[2]))

# Sweep J = V∞/(nD). n = Ω/2π → Ω = 2π·n = 2π·V∞/(J·D).
J_list = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)

@printf "  J     |CT|     |CQ|     KT      10·KQ    η_VLM\n"
@printf "  ───   ──────   ──────   ─────   ──────   ──────\n"
for J in J_list
    n  = Vinf / (J * D)
    Ω  = 2π * n
    r  = rotor_forces(rotor, Vinf, Ω)
    # ship-prop conventions: KT = T/(ρ n² D⁴), KQ = Q/(ρ n² D⁵)
    Sref = π * R^2
    thrust = abs(r.CT) * 0.5 * Vinf^2 * Sref          # ρ=1
    torque = abs(r.CQ) * 0.5 * Vinf^2 * Sref * R
    KT = thrust / (n^2 * D^4)
    KQ = torque / (n^2 * D^5)
    η  = (J * KT) / (2π * KQ)
    @printf "  %4.2f  %.4f   %.4f   %.4f  %.4f   %.4f\n" J abs(r.CT) abs(r.CQ) KT 10*KQ η
end

@printf "\nNote: VLM is inviscid + linear. KT/KQ/η values are\n"
@printf "ideal-fluid bounds; real-prop measurements include viscous\n"
@printf "drag and cavitation losses that shift KT down ~5-15%% and KQ\n"
@printf "up ~10-25%% in a Wageningen B4-70 sense.\n"
