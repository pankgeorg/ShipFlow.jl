#!/usr/bin/env julia
#
# K4: Re sweep with WALE LES enabled. Re-runs the integrated VLM
# stack (Wigley + rotor + rudder) at Re ∈ {5000, 20000, 50000},
# reporting stability (max |u|, settled drag) and a final-frame
# snapshot per Re for visual comparison.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "CairoMakie"]; io=devnull)

using WaterLily, VoF, Turbulence, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics, CairoMakie

const NX = 128; const NY = 64; const NZ = 32
const NSTEPS = 60
const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5); const hull_yc = Float32(NY/2); const hull_zc = H_w_c

const R_prop = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J_op = 0.32; const Ω = Float64(π) * U∞ / (J_op * R_prop)
const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)
const δ_DEG = 5f0

const RE_LIST = [5000.0, 20000.0, 50000.0]

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)
r_rot  = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop

function run_re(Re)
    @info "Run Re=$Re"
    ν_w_c = U∞ * L_c / Re
    μ_w_c = ρ_w * ν_w_c
    μ_a_c = ρ_a * 18 * ν_w_c
    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=Float32(μ_w_c), μ_a=Float32(μ_a_c), T=Float32)
    turb = WALE((NX, NY, NZ); Cw=0.5f0, ν₀=0f0)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν = Turbulence.viscosity(turb),
        g=(i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor, U=U∞,
    )
    function vlm_udf(flow, t; kwargs...)
        smear_blades!(flow.f, thrust, torque,
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1.0, 0.0, 0.0),
                      Float64(R_prop), 0.2*Float64(R_prop);
                      N_blades=3, N_sections=4, ε=1.5)
        r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞)
        side = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        drag = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                     SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
        return nothing
    end
    drag_hist = Float64[]
    umax_hist = Float64[]
    νt_max_hist = Float64[]
    t0 = time()
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        update_νt!(turb, sim.flow.u, vof.ν)
        Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
        push!(drag_hist, -(Float64(Fp[1]) + Float64(Fv[1])))
        push!(umax_hist, maximum(abs, sim.flow.u))
        # ν_t - ν₀ from VoF gives the eddy contribution
        push!(νt_max_hist, maximum(turb.ν))
    end
    @printf "  Re=%-7g  elapsed=%.1fs  D_mean=%.2f  |u|_max=%.3f  ν_t_max=%.4e\n" Re (time()-t0) mean(drag_hist[end-15:end]) maximum(umax_hist) maximum(νt_max_hist)
    return (; Re, drag_hist, umax_hist, νt_max_hist, vof, sim)
end

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "headline_re_sweep"))
mkpath(OUTDIR)

results = [run_re(Re) for Re in RE_LIST]

# Time-series plot
fig = Figure(size=(900, 700))
ax1 = Axis(fig[1, 1]; xlabel="step", ylabel="drag",
    title="K4: Re sweep with WALE LES")
ax2 = Axis(fig[2, 1]; xlabel="step", ylabel="|u|_max")
ax3 = Axis(fig[3, 1]; xlabel="step", ylabel="ν_t_max")
cols = [:steelblue, :tomato, :purple]
for (r, c) in zip(results, cols)
    lines!(ax1, 1:length(r.drag_hist), r.drag_hist;
        color=c, linewidth=2, label="Re=$(Int(r.Re))")
    lines!(ax2, 1:length(r.umax_hist), r.umax_hist;
        color=c, linewidth=2)
    lines!(ax3, 1:length(r.νt_max_hist), r.νt_max_hist;
        color=c, linewidth=2)
end
axislegend(ax1)
save(joinpath(OUTDIR, "Re_sweep.png"), fig)
println("Wrote $(joinpath(OUTDIR, "Re_sweep.png"))")

# Summary
println("\n" * "=" ^ 60)
@printf "  Re        D_mean   |u|_max   ν_t_max\n"
println("=" ^ 60)
for r in results
    @printf "  %-7g   %.2f    %.3f     %.4e\n" r.Re mean(r.drag_hist[end-15:end]) maximum(r.umax_hist) maximum(r.νt_max_hist)
end
println("=" ^ 60)
