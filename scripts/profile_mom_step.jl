#!/usr/bin/env julia
#
# J5: profile per-step allocations for ShipFlow's actual usage pattern.
# WaterLily's `test/alloctest.jl` covers a basic 2D scalar-ν cylinder
# case; this script covers the 3D + array-ν + VoF + udf path that the
# integrated stack actually uses.
#
# Reports the median allocation count (allocs) and bytes allocated
# per step under five configurations:
#
#   A. baseline 3D scalar-ν, no body force
#   B. 3D + VoF (array-ν via vof.ν)
#   C. 3D + VoF + smear_force udf (rotor proxy)
#   D. 3D + VoF + smear_force + smear_torque udf (rotor with swirl)
#   E. 3D + VoF + Wigley body + rotor + rudder udf

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "BenchmarkTools"]; io=devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, BenchmarkTools

const NX = 64; const NY = 32; const NZ = 32   # small grid — we care about allocations, not wall time
const L_c = 24f0; const B_c = 5f0; const T_c = 3f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.25f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 1000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/3); const hull_yc = Float32(NY/2); const hull_zc = H_w_c
const R_prop = 0.8f0 * T_c / 2
const prop_xc = Float32(hull_xc + L_c/2 + T_c/2)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)
const δ = 5f0; const J_op = 0.32
const Ω_rot = Float64(π) * U∞ / (J_op * R_prop)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

# ----------------------------------------------------------------------------
# Case A: baseline 3D scalar-ν, no body
# ----------------------------------------------------------------------------
function setupA()
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=ν_w_c, Δt=0.25f0, ϵ=1, perdir=(2,), exitBC=true, U=U∞)
    return sim
end

# ----------------------------------------------------------------------------
# Case B: 3D + VoF (array-ν)
# ----------------------------------------------------------------------------
function setupB()
    vof = VoFFlow((NX, NY, NZ); α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=vof.ν, Δt=0.25f0, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor, U=U∞,
        g=(i, x, t) -> i == 3 ? -G_c : 0f0)
    return sim, vof
end

# udfs of various complexity
function udf_smear(flow, t; kwargs...)
    smear_force!(flow.f, SVector(50f0, 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    return nothing
end
function udf_smear_torque(flow, t; kwargs...)
    smear_force!(flow.f, SVector(50f0, 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    smear_torque!(flow.f, 5f0,
                  SVector(prop_xc, prop_yc, prop_zc),
                  SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                  N=8, ε=2.0f0)
    return nothing
end
function udf_full(rotor, rudder, thrust, torque)
    function (flow, t; kwargs...)
        smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                     SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
        smear_torque!(flow.f, Float32(torque),
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                      N=8, ε=2.0f0)
        r_rud = rudder_forces(rudder, deg2rad(δ), U∞)
        side = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        drag = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                     SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
        return nothing
    end
end

println("=== Mom_step! allocation profile ===")
@printf "Grid: %d × %d × %d = %.2fM cells\n\n" NX NY NZ NX*NY*NZ/1e6
flush(stdout)

# Case A
simA = setupA()
WaterLily.mom_step!(simA.flow, simA.pois)   # warm-up
bA = @benchmarkable WaterLily.mom_step!($simA.flow, $simA.pois) samples=20 evals=1
rA = run(bA)
@printf "A  scalar-ν, no udf, no body:                %4.1f KiB,  %5d allocs\n" rA.memory/1e3 rA.allocs
flush(stdout)

# Case B
simB, vofB = setupB()
WaterLily.mom_step!(simB.flow, simB.pois)
bB = @benchmarkable WaterLily.mom_step!($simB.flow, $simB.pois) samples=20 evals=1
rB = run(bB)
@printf "B  + VoF (array-ν), no udf:                  %4.1f KiB,  %5d allocs\n" rB.memory/1e3 rB.allocs
flush(stdout)

# Case C
WaterLily.mom_step!(simB.flow, simB.pois; udf=udf_smear)
bC = @benchmarkable WaterLily.mom_step!($simB.flow, $simB.pois; udf=$udf_smear) samples=20 evals=1
rC = run(bC)
@printf "C  + smear_force udf:                        %4.1f KiB,  %5d allocs\n" rC.memory/1e3 rC.allocs
flush(stdout)

# Case D
WaterLily.mom_step!(simB.flow, simB.pois; udf=udf_smear_torque)
bD = @benchmarkable WaterLily.mom_step!($simB.flow, $simB.pois; udf=$udf_smear_torque) samples=20 evals=1
rD = run(bD)
@printf "D  + smear_force + smear_torque:             %4.1f KiB,  %5d allocs\n" rD.memory/1e3 rD.allocs
flush(stdout)

# Case E (full path with body, rotor, rudder)
function setupE()
    vof = VoFFlow((NX, NY, NZ); α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=vof.ν, Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor, U=U∞,
        g=(i, x, t) -> i == 3 ? -G_c : 0f0)
    rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.2*R_prop,
        chord=(0.25*R_prop, 0.18*R_prop),
        twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
    rudder = Rudder(; chord=2.0, span=3.0, ns=12, nc=6)
    r_rot = rotor_forces(rotor, U∞, Ω_rot)
    Sref = π * R_prop^2
    thrust = abs(r_rot.CT) * 0.5 * U∞^2 * Sref
    torque = r_rot.CQ * 0.5 * U∞^2 * Sref * R_prop
    return sim, vof, udf_full(rotor, rudder, thrust, torque)
end

simE, vofE, udfE = setupE()
WaterLily.mom_step!(simE.flow, simE.pois; udf=udfE)
bE = @benchmarkable WaterLily.mom_step!($simE.flow, $simE.pois; udf=$udfE) samples=20 evals=1
rE = run(bE)
@printf "E  + Wigley body + rotor + rudder (VLM):     %4.1f KiB,  %5d allocs\n" rE.memory/1e3 rE.allocs

# Step_vof_mules separately
bV = @benchmarkable step_vof_mules!($vofE, $simE; dt=$simE.flow.Δt[end-1], perdir=(2,)) samples=20 evals=1
rV = run(bV)
@printf "    step_vof_mules! (separate from mom_step): %4.1f KiB,  %5d allocs\n" rV.memory/1e3 rV.allocs

println("\nNotes:")
println("- WaterLily upstream baseline (2D scalar-ν cylinder) is ~2 KiB per mom_step!.")
println("- VLM-via-VortexLattice in udf is the major allocator: each rudder_forces call")
println("  builds a fresh VortexLattice System. Future optimisation: cache the System.")
