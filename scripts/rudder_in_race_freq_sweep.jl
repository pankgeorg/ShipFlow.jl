#!/usr/bin/env julia
#
# J2: sinusoidal δ(t) frequency sweep with the rudder INSIDE the
# rotor race (F2 placement) and two-way coupling on. Cross of F2
# and I4: how does the amplified static gain show up in the
# unsteady transfer function?

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "CairoMakie"]; io=devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = 128; const NY = 64; const NZ = 32
const δ_MAX_DEG = 10f0
const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5); const hull_yc = Float32(NY/2); const hull_zc = H_w_c
const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J_op = 0.32
const Ω_rot = Float64(π) * U∞ / (J_op * R_prop)
# IN-RACE rudder placement (F2 geometry)
const rud_xc = Float32(prop_xc + R_prop * 0.5)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

const PERIOD_STEPS = [80, 40, 20]   # frequencies
const N_PERIODS    = 3

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)
r_rot  = rotor_forces(rotor, U∞, Ω_rot)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop

function vof_pois_ctor(vof)
    (flow) -> begin
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
end

function run_freq(period_steps::Int)
    nsteps = period_steps * N_PERIODS
    @info @sprintf("Running f: period=%d steps, total=%d steps", period_steps, nsteps)
    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=vof.ν,
        g=(i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor(vof), U=U∞,
    )
    δ_hist  = Float64[]
    CL_hist = Float64[]
    Fhy_hist = Float64[]
    CL_buf = Float32[]
    for step in 1:nsteps
        ϕ = 2π * (step - 1) / period_steps
        δ_deg = Float64(δ_MAX_DEG) * sin(ϕ)
        function combo_udf(flow, t; kwargs...)
            smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                         SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
            smear_torque!(flow.f, Float32(torque),
                          SVector(prop_xc, prop_yc, prop_zc),
                          SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                          N=8, ε=2.0f0)
            # F2 placement + two-way coupling
            inflow_fn = trilinear_inflow(flow.u)
            r_rud = rudder_forces(rudder, deg2rad(δ_deg), U∞; inflow=inflow_fn)
            side = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
            drag = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
            smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                         SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
            push!(CL_buf, Float32(r_rud.CL))
            return nothing
        end
        WaterLily.mom_step!(sim.flow, sim.pois; udf=combo_udf, pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
        push!(δ_hist, δ_deg)
        push!(CL_hist, length(CL_buf) > 0 ? Float64(CL_buf[end]) : 0.0)
        empty!(CL_buf)
        push!(Fhy_hist, Float64(Fp[2] + Fv[2]))
    end
    # Fourier projection at ω = 2π/period (after 1-period burn-in)
    keep = period_steps + 1 : nsteps
    δ_k  = δ_hist[keep]
    CL_k = CL_hist[keep]
    Fy_k = Fhy_hist[keep]
    N = length(δ_k)
    cs = [cos(2π * (i - 1) / period_steps) for i in 1:N]
    sn = [sin(2π * (i - 1) / period_steps) for i in 1:N]
    fourier(x) = (2/N * sum(x .* cs), 2/N * sum(x .* sn))
    cδ, sδ = fourier(δ_k);  Aδ = hypot(cδ, sδ); φδ = atan(sδ, cδ)
    cC, sC = fourier(CL_k); AC = hypot(cC, sC); φC = atan(sC, cC)
    cF, sF = fourier(Fy_k); AF = hypot(cF, sF); φF = atan(sF, cF)
    return (; period_steps, Aδ, AC, AF,
              gain = AF / max(AC, 1e-9), phase = (φF - φC) * 180/π)
end

results = [run_freq(p) for p in PERIOD_STEPS]

println("\n══════════════════════════════════════════════════════════════════════")
@printf "  period   freq      |Aδ|     |ACL|    |AFy|    gain    phase(°)\n"
println("══════════════════════════════════════════════════════════════════════")
for r in results
    f = 1.0 / (r.period_steps * 0.25)
    @printf "  %4d     %.4f    %.4f   %.4f   %.4f   %.3f   %+7.1f\n" r.period_steps f r.Aδ r.AC r.AF r.gain r.phase
end
println("══════════════════════════════════════════════════════════════════════")

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder_in_race_freq"))
mkpath(OUTDIR)
freqs = [1.0 / (r.period_steps * 0.25) for r in results]
gains = [r.gain for r in results]
phases = [r.phase for r in results]
ACs = [r.AC for r in results]

# Compare to I4 baseline: read the rudder_freq_sweep summary if present.
function read_i4_csv()
    p = abspath(joinpath(@__DIR__, "..", "runs", "rudder_freq_sweep"))
    isdir(p) || return nothing
    return nothing  # data not persisted in CSV in I4; skip for now
end

fig = Figure(size=(900, 560))
ax1 = Axis(fig[1, 1]; xlabel="frequency (1/cell-time)",
    ylabel="|A_CL| amplitude",
    title="I4 vs J2: rudder CL amplitude (out-of-race vs in-race)")
scatter!(ax1, freqs, ACs; color=:tomato, markersize=14, label="J2 (in race)")
lines!(ax1, freqs, ACs; color=:tomato)
# Quick reference numbers from I4 (period 80, 40, 20 had A_CL = 0.303 each)
A_CL_i4 = fill(0.303, length(freqs))
scatter!(ax1, freqs, A_CL_i4; color=:steelblue, markersize=14, label="I4 (out of race)")
lines!(ax1, freqs, A_CL_i4; color=:steelblue, linestyle=:dash)
axislegend(ax1)

ax2 = Axis(fig[2, 1]; xlabel="frequency (1/cell-time)",
    ylabel="gain |A_F/A_CL|",
    title="Manoeuvring gain (in race)")
scatter!(ax2, freqs, gains; color=:purple, markersize=14)
lines!(ax2, freqs, gains; color=:purple)
out = joinpath(OUTDIR, "bode_in_race.png")
save(out, fig)
println("Wrote $out")
