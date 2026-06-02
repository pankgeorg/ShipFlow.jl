#!/usr/bin/env julia
#
# L3: hydrostatic pressure sanity check. With the hull fixed and the
# inflow essentially zero, the only force on the body should be the
# Archimedes upthrust ρ_w·g·V_sub. This script verifies that
# `WaterLily.pressure_force(sim)[3]` converges to that value (with
# the body→fluid sign convention: +Archimedes on body ⇒ negative
# Fp[3]).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)

using WaterLily, VoF, ShipShapes
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie

const NX = 128; const NY = 64; const NZ = 32
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "200"))
const L_c = 36.0; const B_c = 8.0; const T_c = 5.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞  = 1e-4   # essentially zero inflow
const Fr  = 0.25
const G_c = 1.0 / (Fr^2 * L_c)   # use U=1 for gravity scaling (consistent)
const Re  = 5000.0
const ν_w_c = 1.0 * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c

const V0 = 4 * L_c * B_c * T_c / 9
const F_Archimedes_expected = ρ_w * G_c * V0

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L=L_c, B=B_c, T=T_c, map=hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== L3: hydrostatic sanity check ===\n"
@printf "  Grid     = %d × %d × %d  NSTEPS=%d\n" NX NY NZ NSTEPS
@printf "  L=%.1f B=%.1f T=%.1f  V₀=%.2f  ρ_w·g·V₀=%.2f (expected)\n" L_c B_c T_c V0 F_Archimedes_expected
@printf "  Inflow U∞=%.0e (≈ zero)\n\n" U∞
flush(stdout)

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
    T=Float64, ν = VoF.viscosity(vof),
    g=(i, x, t) -> i == 3 ? -G_c : 0.0,
    Δt=0.25, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=1.0,    # reference scale only
)

Fp_hist = Float64[]
err_pct_hist = Float64[]
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1e-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    Fp = WaterLily.pressure_force(sim)
    F_on_body = -Fp[3]    # body→fluid → fluid→body
    err = (F_on_body - F_Archimedes_expected) / F_Archimedes_expected * 100
    push!(Fp_hist, F_on_body)
    push!(err_pct_hist, err)
    if step % 25 == 0
        @printf "  step=%3d  F_on_body=%.2f  expected=%.2f  err=%+7.2f%%  |u|=%.4e\n" step F_on_body F_Archimedes_expected err maximum(abs, sim.flow.u)
        flush(stdout)
    end
end

# Final values from last 25 %
tail = (NSTEPS - NSTEPS÷4 + 1):NSTEPS
F_final = sum(Fp_hist[tail]) / length(tail)
err_final = (F_final - F_Archimedes_expected) / F_Archimedes_expected * 100
@printf "\nFinal F_on_body (last 25%%) = %.2f\n" F_final
@printf "Expected ρ_w·g·V₀          = %.2f\n" F_Archimedes_expected
@printf "Relative error             = %+.2f %%\n" err_final
verdict = abs(err_final) < 5 ? "PASS" : "FAIL"
@printf "Verdict (< 5%% target):     = %s\n" verdict

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "hydrostatic_sanity"))
mkpath(OUTDIR)
fig = Figure(size=(800, 480))
ax = Axis(fig[1, 1]; xlabel="step", ylabel="F_z on body",
    title="L3: hydrostatic Archimedes verification")
lines!(ax, 1:NSTEPS, Fp_hist; color=:steelblue, linewidth=2, label="measured")
hlines!(ax, [F_Archimedes_expected]; color=:tomato, linestyle=:dash, linewidth=2, label="ρ_w·g·V₀")
axislegend(ax)
save(joinpath(OUTDIR, "hydrostatic.png"), fig)
println("Wrote $(joinpath(OUTDIR, "hydrostatic.png"))")
