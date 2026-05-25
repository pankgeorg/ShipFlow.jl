#!/usr/bin/env julia
#
# Combined rudder × propeller polar via VLM. The pure-VLM tier (no
# WaterLily). Sweeps δ over rudder angles, J over advance ratios,
# tabulates the (CL, CD) of the rudder + (CT, CQ, η) of the propeller
# in one table — the inputs a maneuvering-model coupler needs.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl"))]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)

using LiftingSurfaces
using Printf

const Vinf = 1.0
const R    = 1.0
const D    = 2 * R

rudder = Rudder(; chord=1.0, span=2.0, ns=12, nc=6)
rotor  = BladedRotor(; N_blades=3, R=R, R_hub=0.2,
                    chord=(0.25, 0.18),
                    twist=(deg2rad(35.0), deg2rad(15.0)),
                    ns=12, nc=4)

δ_list = (-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0)
J_list = (0.5, 0.7, 0.9, 1.1)

@printf "=== Combined rudder × propeller polar ===\n"
@printf "  Rudder: AR=2 rectangular flat-plate, %d×%d panels\n" 12 6
@printf "  Propeller: 3×VLM blades, %d×%d panels/blade\n\n" rotor.ns rotor.nc

@printf "\nRudder polar (J-independent; depends only on δ):\n"
@printf "    δ    CL       CD       L/D\n"
@printf "    ──   ──────   ──────   ──────\n"
for δ in δ_list
    r = rudder_forces(rudder, deg2rad(δ), Vinf)
    LoD = abs(r.CD) > 1e-9 ? r.CL / r.CD : Inf
    @printf "  %+5.1f  %+8.4f  %+8.5f  %+7.2f\n" δ r.CL r.CD LoD
end

@printf "\nPropeller polar (δ-independent; depends only on J):\n"
@printf "   J     |CT|      |CQ|      KT      10·KQ    η_VLM\n"
@printf "  ────   ──────    ──────    ──────  ──────   ──────\n"
for J in J_list
    n = Vinf / (J * D)
    Ω = 2π * n
    r = rotor_forces(rotor, Vinf, Ω)
    Sref = π * R^2
    thrust = abs(r.CT) * 0.5 * Vinf^2 * Sref
    torque = abs(r.CQ) * 0.5 * Vinf^2 * Sref * R
    KT = thrust / (n^2 * D^4)
    KQ = torque / (n^2 * D^5)
    η  = (J * KT) / (2π * KQ)
    @printf "  %4.2f  %8.4f   %.4f   %.4f  %.4f   %.4f\n" J abs(r.CT) abs(r.CQ) KT 10*KQ η
end

@printf "\nIn the pure-VLM tier the two surfaces are decoupled (each sees\n"
@printf "only the freestream Vinf=1). In a coupled WaterLily run, the\n"
@printf "rudder VLM samples the actual local flow (rotor race + hull wake)\n"
@printf "via trilinear_inflow, so the table above becomes a baseline that\n"
@printf "the coupled CL/CD/CT/CQ get *perturbed* from. See\n"
@printf "wigley_full_stack_VLM.jl for the coupled path.\n"
