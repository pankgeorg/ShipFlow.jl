#!/usr/bin/env julia
#
# Numerical comparison: Propellers.SwirlingDisk vs LiftingSurfaces.BladedRotor
# in the same Wigley + free-surface setup. Same prescribed thrust;
# compare hull drag and free-surface wave structure.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Propellers.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes, Propellers, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics

const NX, NY, NZ = 128, 64, 32
const NSTEPS = 80
const L_c = 32f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/4)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

const R_prop = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J = 0.7
const Ω = Float64(π) * U∞ / (J * R_prop)

# Calibrate SwirlingDisk thrust to match what BladedRotor's VLM produces,
# so we compare LIKE WITH LIKE.
rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)),
    ns=12, nc=4)
r_blade  = rotor_forces(rotor, U∞, Ω)
thrust_T = abs(r_blade.CT * 0.5 * U∞^2 * π * R_prop^2)
torque_T = abs(r_blade.CQ * 0.5 * U∞^2 * π * R_prop^2 * R_prop)
@printf "BladedRotor at J=%.2f:  |CT|=%.4f  → thrust=%.3f  |CQ|=%.4f  → torque=%.3f\n" J abs(r_blade.CT) thrust_T abs(r_blade.CQ) torque_T

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

function run_case(label, udf)
    @info "Running case: $label"
    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T = Float32, ν = vof.ν,
        g = (i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = U∞)
    drag_hist = Float64[]
    t0 = time()
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; udf=udf, pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
        push!(drag_hist, -(Float64(Fp[1]) + Float64(Fv[1])))
        if step % 20 == 0
            @printf "  %s  step=%3d  |u|=%.3f  drag=%.3f  elapsed=%.1fs\n" label step maximum(abs, sim.flow.u) drag_hist[end] time()-t0
            flush(stdout)
        end
    end
    last_third = drag_hist[2*NSTEPS÷3+1:end]
    drag_m = mean(last_third); drag_s = std(last_third)
    # Wave RMS past stern
    nx, ny, nz = size(vof.α)
    eta = Float64[]
    for j in 2:ny-1, i in round(Int, hull_xc + L_c/2 + 2):nx-1
        for k in nz-1:-1:2
            if vof.α[i,j,k] ≥ 0.5 && vof.α[i,j,k+1] < 0.5
                Δa = vof.α[i,j,k] - vof.α[i,j,k+1]
                t = abs(Δa) > 1e-9 ? (0.5 - vof.α[i,j,k+1]) / Δa : 0.5
                z = (k+1 - 1.5) - t
                push!(eta, z - H_w_c)
                break
            end
        end
    end
    η_rms = isempty(eta) ? NaN : sqrt(mean(abs2, eta .- mean(eta)))
    η_pp  = isempty(eta) ? NaN : maximum(eta) - minimum(eta)
    return (; label, drag_m, drag_s, η_rms, η_pp)
end

# SwirlingDisk udf with matched thrust + torque
disk = SwirlingDisk(center=SVector(prop_xc, prop_yc, prop_zc),
    axis=SVector(1f0, 0f0, 0f0), R=R_prop, w=1.5f0,
    thrust=Float32(thrust_T), torque=Float32(torque_T))
udf_swirl = (flow, t; kw...) -> disk(flow, t; kw...)

# BladedRotor udf (smear thrust at disk centre)
function udf_bladed(flow, t; kw...)
    r = rotor_forces(rotor, U∞, Ω)
    thrust = abs(r.CT * 0.5 * U∞^2 * π * R_prop^2)
    smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    return nothing
end

r_swirl = run_case("SWIRL  ", udf_swirl)
r_blade = run_case("BLADED ", udf_bladed)

println()
@printf "═══════════════════════════════════════════════════════════════\n"
@printf "                          SWIRL          BLADED         Δ%%\n"
@printf "  hull drag (mean ± σ)   %.3f±%.3f   %.3f±%.3f   %+.1f%%\n" r_swirl.drag_m r_swirl.drag_s r_blade.drag_m r_blade.drag_s (r_blade.drag_m/r_swirl.drag_m - 1)*100
@printf "  wave RMS (post-stern)  %.4f         %.4f         %+.1f%%\n" r_swirl.η_rms r_blade.η_rms (r_blade.η_rms/r_swirl.η_rms - 1)*100
@printf "  wave peak-peak         %.4f         %.4f         %+.1f%%\n" r_swirl.η_pp  r_blade.η_pp  (r_blade.η_pp/r_swirl.η_pp - 1)*100
@printf "═══════════════════════════════════════════════════════════════\n"
