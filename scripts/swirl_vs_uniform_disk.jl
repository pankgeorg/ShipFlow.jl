#!/usr/bin/env julia
#
# Compare ActuatorDisk (uniform thrust) vs SwirlingDisk (thrust+torque)
# behind a Wigley hull, same C_T. Question: does swirl change the
# converged hull drag and free-surface wave amplitude?

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Propellers.jl")),
]; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes, Turbulence, Propellers
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "180"))

const L_c = 48.0; const B_c = 10.0; const T_c = 6.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞ = 1.0
const Fr = 0.30
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c

const H_w_c   = NZ/2
const hull_xc = NX/3
const hull_yc = NY/2
const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = NY/2
const prop_zc = H_w_c - T_c/2
const R_prop = T_c/2 * 0.6
const W_prop = 1.5
const C_T = 0.6
const thrust = 0.5 * π * R_prop^2 * U∞^2 * C_T
const torque = 0.5 * thrust * R_prop   # Q/(T·R) = 0.5

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== ActuatorDisk vs SwirlingDisk (same C_T) ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Hull       = L=%.1f at xc=%.1f\n" L_c hull_xc
@printf "  Disk       = R=%.2f at (%.1f, %.1f, %.1f)\n" R_prop prop_xc prop_yc prop_zc
@printf "  C_T=%.2f  thrust=%.3f  torque=%.3f  Q/(T·R)=0.5\n" C_T thrust torque
@printf "  Steps      = %d\n\n" NSTEPS
flush(stdout)

function run_case(label, disk_obj)
    @info "Running case: $label"
    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float64)
    turb = WALE((NX, NY, NZ); Cw = 0.5, ν₀ = 0.0)

    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end

    sim = WaterLily.Simulation((NX, NY, NZ),
        (U∞, 0.0, 0.0), L_c;
        T = Float64, ν = turb.ν,
        g = (i, x, t) -> i == 3 ? -G_c : 0.0,
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = U∞,
    )
    disk_udf = (flow, t; kwargs...) -> disk_obj(flow, t; kwargs...)

    drag_hist = Float64[]
    t0 = time()
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf,
            pois_tol=1e-7, pois_itmx=80)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        update_νt!(turb, sim.flow.u, vof.ν)
        Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
        push!(drag_hist, -(Fp[1] + Fv[1]))
        if step % 30 == 0
            @printf "  %s  step=%3d  |u|=%.3f  drag=%.3f  elapsed=%.1fs\n" label step maximum(abs, sim.flow.u) drag_hist[end] time()-t0
            flush(stdout)
        end
    end
    last_third = drag_hist[2*NSTEPS÷3+1:end]
    drag_mean = mean(last_third); drag_std = std(last_third)
    # Wave RMS in domain past stern
    α = vof.α
    nx, ny, nz = size(α)
    eta_vals = Float64[]
    for j in 2:ny-1, i in round(Int, hull_xc + L_c/2 + 2):nx-1
        for k in nz-1:-1:2
            if α[i,j,k] ≥ 0.5 && α[i,j,k+1] < 0.5
                Δα = α[i,j,k] - α[i,j,k+1]
                t = abs(Δα) > 1e-9 ? (0.5 - α[i,j,k+1]) / Δα : 0.5
                z = (k+1 - 1.5) - t
                push!(eta_vals, z - H_w_c)
                break
            end
        end
    end
    η_rms = isempty(eta_vals) ? NaN : sqrt(mean(abs2, eta_vals .- mean(eta_vals)))
    η_pp  = isempty(eta_vals) ? NaN : maximum(eta_vals) - minimum(eta_vals)
    return (; label, drag_mean, drag_std, η_rms, η_pp, n_eta=length(eta_vals))
end

# Build both disks (axis=+x)
disk_uniform = ActuatorDisk(
    center = SVector(prop_xc, prop_yc, prop_zc),
    axis   = SVector(1.0, 0.0, 0.0),
    R = R_prop, w = W_prop, thrust = thrust,
)
disk_swirl = SwirlingDisk(
    center = SVector(prop_xc, prop_yc, prop_zc),
    axis   = SVector(1.0, 0.0, 0.0),
    R = R_prop, w = W_prop,
    thrust = thrust, torque = torque,
)

r_uni = run_case("UNIFORM", disk_uniform)
r_swl = run_case("SWIRL  ", disk_swirl)

println()
@printf "═══════════════════════════════════════════════════════════════\n"
@printf "                          UNIFORM        SWIRL          Δ%%\n"
@printf "  hull drag (mean ±σ)   %+.3f±%.3f   %+.3f±%.3f    %+.1f%%\n" r_uni.drag_mean r_uni.drag_std r_swl.drag_mean r_swl.drag_std (r_swl.drag_mean/r_uni.drag_mean - 1)*100
@printf "  wave RMS (post-stern) %.4f         %.4f          %+.1f%%\n" r_uni.η_rms r_swl.η_rms (r_swl.η_rms/r_uni.η_rms - 1)*100
@printf "  wave peak-peak        %.4f         %.4f          %+.1f%%\n" r_uni.η_pp  r_swl.η_pp  (r_swl.η_pp/r_uni.η_pp - 1)*100
@printf "  thrust (prescribed)   %.3f          %.3f\n" thrust thrust
@printf "  torque (prescribed)   0.000          %.3f  (Q/(T·R)=0.5)\n" torque
@printf "═══════════════════════════════════════════════════════════════\n"
