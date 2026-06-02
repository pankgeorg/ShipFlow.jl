#!/usr/bin/env julia
#
# I4: sinusoidal rudder frequency sweep. Drive δ(t) = δ_max · sin(ω t)
# for several ω, record F_hull,y(t) and CL_rud(t), and compute the
# manoeuvring transfer function (gain + phase) by single-frequency
# Fourier projection.
#
# Frequency convention: period given in STEPS (the cell-time integrator
# step). dt = 0.25 cell-time-unit, so period_in_cells = period_steps · dt.

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
const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

const PERIOD_STEPS = [80, 40, 20, 12]   # → frequencies decreasing
const N_PERIODS    = 3                    # cycles per frequency
const BURNIN_FRAC  = 1                    # discard first period

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
        T=Float32, ν = VoF.viscosity(vof),
        g=(i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor(vof), U=U∞,
    )
    δ_hist  = Float64[]
    CL_hist = Float64[]
    Fhy_hist = Float64[]
    CL_buf = Float32[]
    for step in 1:nsteps
        # Phase ∈ [0, 2π) over `period_steps`.
        ϕ = 2π * (step - 1) / period_steps
        δ_deg = Float64(δ_MAX_DEG) * sin(ϕ)
        function combo_udf(flow, t; kwargs...)
            smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                         SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
            smear_torque!(flow.f, Float32(torque),
                          SVector(prop_xc, prop_yc, prop_zc),
                          SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                          N=8, ε=2.0f0)
            r_rud = rudder_forces(rudder, deg2rad(δ_deg), U∞)
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
        push!(δ_hist,  δ_deg)
        push!(CL_hist, length(CL_buf) > 0 ? Float64(CL_buf[end]) : 0.0)
        empty!(CL_buf)
        push!(Fhy_hist, Float64(Fp[2] + Fv[2]))
    end
    # Discard first period for transient.
    keep = (period_steps * BURNIN_FRAC + 1):nsteps
    δ_k  = δ_hist[keep]
    CL_k = CL_hist[keep]
    Fy_k = Fhy_hist[keep]
    # Single-frequency Fourier projection at ω = 2π/period.
    N = length(δ_k)
    cs = [cos(2π * (i - 1) / period_steps) for i in 1:N]
    sn = [sin(2π * (i - 1) / period_steps) for i in 1:N]
    fourier(x) = (2/N * sum(x .* cs), 2/N * sum(x .* sn))
    cδ, sδ = fourier(δ_k);  Aδ = hypot(cδ, sδ); φδ = atan(sδ, cδ)
    cC, sC = fourier(CL_k); AC = hypot(cC, sC); φC = atan(sC, cC)
    cF, sF = fourier(Fy_k); AF = hypot(cF, sF); φF = atan(sF, cF)
    gain    = AF / AC
    phase   = (φF - φC) * 180 / π   # degrees, F relative to CL
    return (; period_steps, Aδ, AC, AF, gain, phase,
            δ_hist=δ_k, CL_hist=CL_k, Fy_hist=Fy_k)
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

# Bode-style plot
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder_freq_sweep"))
mkpath(OUTDIR)
freqs = [1.0 / (r.period_steps * 0.25) for r in results]
gains = [r.gain for r in results]
phases = [r.phase for r in results]

fig = Figure(size=(900, 600))
ax1 = Axis(fig[1, 1]; xlabel="frequency (1/cell-time)",
    ylabel="|F_hull,y| / |CL_rud|",
    title="Manoeuvring response — gain")
scatter!(ax1, freqs, gains; color=:tomato, markersize=14)
lines!(ax1, freqs, gains; color=:tomato)
ax2 = Axis(fig[2, 1]; xlabel="frequency (1/cell-time)",
    ylabel="phase (°)", title="Phase of F_hull,y relative to CL_rud")
scatter!(ax2, freqs, phases; color=:steelblue, markersize=14)
lines!(ax2, freqs, phases; color=:steelblue)
hlines!(ax2, [0]; color=:grey, linestyle=:dash)
out = joinpath(OUTDIR, "bode.png")
save(out, fig)
println("Wrote $out")

# Time series figure for each frequency, stacked
fig2 = Figure(size=(1100, 1000))
for (i, r) in enumerate(results)
    ax = Axis(fig2[i, 1]; xlabel="step", ylabel="signal",
        title=@sprintf("period=%d steps  f=%.4f  gain=%.3f  phase=%.1f°",
            r.period_steps, 1.0/(r.period_steps*0.25), r.gain, r.phase))
    xs = collect(1:length(r.δ_hist))
    lines!(ax, xs, r.δ_hist;
        color=:black, linewidth=1.5, label="δ (°)")
    lines!(ax, xs, r.CL_hist .* 5;
        color=:tomato, linewidth=1.5, label="5 · CL_rud")
    lines!(ax, xs, r.Fy_hist;
        color=:steelblue, linewidth=1.5, label="F_hull,y")
    hlines!(ax, [0]; color=:grey, linestyle=:dash)
    axislegend(ax; position=:rt)
end
out2 = joinpath(OUTDIR, "timeseries.png")
save(out2, fig2)
println("Wrote $out2")
