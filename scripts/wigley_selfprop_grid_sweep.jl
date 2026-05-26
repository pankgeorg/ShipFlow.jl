#!/usr/bin/env julia
#
# J4: grid-resolution sensitivity for J_self. Re-runs the Wigley
# self-propulsion sweep (G1 / NEW_PLAN F1 procedure) at three grid
# resolutions: 64×32×32, 96×48×48 (baseline), 128×64×64. Reports
# how J_self shifts with grid; tight clustering = grid-converged.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf

const NSTEPS_PER = 80
const AVG_FRAC = 0.25
const J_LIST = [0.20, 0.25, 0.30, 0.35, 0.45]

# Grid set — proportional, hull occupies same fraction of each.
const GRIDS = [
    (NX=64,  NY=32, NZ=32, L_c=24.0, B_c=5.0, T_c=3.0),
    (NX=96,  NY=48, NZ=48, L_c=48.0, B_c=10.0, T_c=6.0),
    (NX=128, NY=64, NZ=64, L_c=64.0, B_c=14.0, T_c=8.0),
]

const ρ_w = 10.0; const ρ_a = 1.0
const U∞ = 1.0
const Fr = 0.25
const Re = 1000.0

function run_one(grid, J)
    NX, NY, NZ, L_c, B_c, T_c = grid.NX, grid.NY, grid.NZ, grid.L_c, grid.B_c, grid.T_c
    H_w_c = NZ/2
    hull_xc = NX/3; hull_yc = NY/2; hull_zc = H_w_c
    R_prop = 0.8 * T_c / 2
    prop_xc = hull_xc + L_c/2 + T_c/2
    prop_yc = hull_yc; prop_zc = H_w_c - T_c/2
    G_c = U∞^2 / (Fr^2 * L_c)
    ν_w_c = U∞ * L_c / Re
    μ_w_c = ρ_w * ν_w_c
    μ_a_c = ρ_a * 18 * ν_w_c
    Ω = π * U∞ / (J * R_prop)

    α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
    hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
    hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
    rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.2*R_prop,
        chord=(0.25*R_prop, 0.18*R_prop),
        twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)

    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float64)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0.0, 0.0), L_c;
        T=Float64, ν=vof.ν,
        g=(i, x, t) -> i == 3 ? -G_c : 0.0,
        Δt=0.25, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor, U=U∞,
    )
    r_rot = rotor_forces(rotor, U∞, Ω)
    Sref = π * R_prop^2
    thrust = abs(r_rot.CT) * 0.5 * U∞^2 * Sref
    torque = r_rot.CQ * 0.5 * U∞^2 * Sref * R_prop
    function rotor_udf(flow, t; kwargs...)
        smear_force!(flow.f, SVector(thrust, 0.0, 0.0),
                     SVector(prop_xc, prop_yc, prop_zc); ε=2.5)
        smear_torque!(flow.f, torque,
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1.0, 0.0, 0.0), 0.7 * R_prop;
                      N=8, ε=2.0)
        return nothing
    end
    drags = Float64[]
    for step in 1:NSTEPS_PER
        WaterLily.mom_step!(sim.flow, sim.pois; udf=rotor_udf, pois_tol=1e-6, pois_itmx=50)
        step_vof!(vof, sim; dt=sim.flow.Δt[end-1])
        if step > NSTEPS_PER * (1 - AVG_FRAC)
            Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
            push!(drags, -(Fp[1] + Fv[1]))
        end
    end
    drag_avg = sum(drags) / length(drags)
    return abs(r_rot.CT), thrust, drag_avg
end

function find_sign_change(rows)
    for i in 1:length(rows)-1
        d1, d2 = rows[i].T - rows[i].D, rows[i+1].T - rows[i+1].D
        if d1 * d2 < 0
            t = -d1 / (d2 - d1)
            return rows[i].J + t * (rows[i+1].J - rows[i].J)
        end
    end
    return nothing
end

results = []
for grid in GRIDS
    @printf "\n=== Grid %d×%d×%d  (L=%.1f, B=%.1f, T=%.1f) ===\n" grid.NX grid.NY grid.NZ grid.L_c grid.B_c grid.T_c
    rows = []
    t0 = time()
    for J in J_LIST
        CT, T_val, D_val = run_one(grid, J)
        push!(rows, (J=J, CT=CT, T=T_val, D=D_val))
        @printf "  J=%.2f  CT=%.3f  T=%.2f  D=%.2f  T-D=%+.2f  elapsed=%.1fs\n" J CT T_val D_val (T_val-D_val) (time()-t0)
        flush(stdout)
    end
    J_sp = find_sign_change(rows)
    push!(results, (grid=grid, J_sp=J_sp, rows=rows))
    @printf "  J_self ≈ %s\n" (J_sp === nothing ? "no bracket" : @sprintf("%.4f", J_sp))
end

println("\n" * "=" ^ 70)
@printf "  Grid             cells   J_self\n"
println("=" ^ 70)
for r in results
    g = r.grid
    @printf "  %d×%d×%d     %.2fM   %s\n" g.NX g.NY g.NZ (g.NX*g.NY*g.NZ/1e6) (
        r.J_sp === nothing ? "no bracket" : @sprintf("%.4f", r.J_sp))
end
println("=" ^ 70)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "selfprop_grid_sweep"))
mkpath(OUTDIR)
open(joinpath(OUTDIR, "summary.csv"), "w") do io
    println(io, "NX,NY,NZ,cells,J_self")
    for r in results
        g = r.grid
        Js = r.J_sp === nothing ? "NaN" : string(r.J_sp)
        println(io, "$(g.NX),$(g.NY),$(g.NZ),$(g.NX*g.NY*g.NZ),$Js")
    end
end
println("Summary CSV written to $(joinpath(OUTDIR, "summary.csv"))")
