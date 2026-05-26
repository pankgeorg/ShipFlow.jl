#!/usr/bin/env julia
#
# I5: Wigley wave-resistance curve at the current substrate.
# Sweeps Fr ∈ {0.20, 0.25, 0.30, 0.35, 0.40, 0.45}, computes total
# resistance coefficient Cw_total = D / (0.5·ρ·U²·S_wet), and writes
# the curve. The classical Wigley result (Bai 1979, Inui 1980) shows a
# hump near Fr ≈ 0.30–0.32 with Cw ≈ 4×10⁻³ at the peak.
#
# We're at Re = 5000, so the viscous component is large compared to a
# real Re = 10⁹ ship. The reported number is *total* drag including
# friction — a clean wave-resistance comparison would require either
# subtracting the ITTC friction line or running a no-wave (deep-water)
# reference.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
]; io=devnull)
Pkg.add(["CairoMakie", "DelimitedFiles"]; io=devnull)

using WaterLily, VoF, ShipShapes
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics, CairoMakie

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "120"))
const FR_LIST = let s = get(ENV, "WL_FR_LIST", "0.20,0.25,0.30,0.35,0.40,0.45")
    [parse(Float64, x) for x in split(s, ",")]
end

const L_c = 36.0; const B_c = 8.0; const T_c = 5.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞  = 1.0
const Re  = 5000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/5; const hull_yc = NY/2; const hull_zc = H_w_c
const S_wet = containership_volume(L_c, B_c, T_c) / T_c  # crude — V/T ≈ planform area
# A better estimate: full Wigley wetted area ≈ 1.18·L·(B + 2T) for half-form.
# Using V/T here gives the right order of magnitude (~150-200 cells²).

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c, map = hull_map)

@printf "=== I5: Wigley resistance vs Fr ===\n"
@printf "  Grid       = %d × %d × %d   NSTEPS=%d\n" NX NY NZ NSTEPS
@printf "  Hull       = L=%.1f, B=%.1f, T=%.1f   S_ref≈%.2f\n" L_c B_c T_c S_wet
@printf "  Re=%.0f   FR_LIST=%s\n\n" Re join(FR_LIST, ", ")
flush(stdout)

function run_fr(Fr)
    G_c = U∞^2 / (Fr^2 * L_c)
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
    Dp_hist = Float64[]
    Dv_hist = Float64[]
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1e-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        if step > NSTEPS - NSTEPS÷4
            Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
            push!(Dp_hist, -Fp[1])
            push!(Dv_hist, -Fv[1])
        end
    end
    Dp = mean(Dp_hist); Dv = mean(Dv_hist)
    Cwp = Dp / (0.5 * ρ_w * U∞^2 * S_wet)
    Cwv = Dv / (0.5 * ρ_w * U∞^2 * S_wet)
    return Dp, Dv, Cwp, Cwv
end

results = []
csv = open(joinpath(abspath(@__DIR__), "..", "runs", "containership_fr_resistance.csv"), "w")
println(csv, "Fr,Dp,Dv,D_total,Cwp,Cwv,Cw_total")
for Fr in FR_LIST
    @printf "  Fr=%.2f …  " Fr
    flush(stdout)
    Dp, Dv, Cwp, Cwv = run_fr(Fr)
    push!(results, (Fr=Fr, Dp=Dp, Dv=Dv, Cwp=Cwp, Cwv=Cwv,
                    Cw_tot=Cwp + Cwv))
    @printf "Dp=%.2f  Dv=%.2f  D_tot=%.2f  Cw_p=%.4e  Cw_v=%.4e\n" Dp Dv (Dp+Dv) Cwp Cwv
    println(csv, "$Fr,$Dp,$Dv,$(Dp+Dv),$Cwp,$Cwv,$(Cwp+Cwv)")
    flush(csv)
end
close(csv)

println()
@printf "═══════════════════════════════════════════════════════════════════════════\n"
@printf "  Fr     Dp        Dv      D_total   Cwp           Cwv           Cw_total\n"
@printf "═══════════════════════════════════════════════════════════════════════════\n"
for r in results
    @printf "  %.2f   %7.3f   %7.3f   %7.3f    %.4e    %.4e    %.4e\n" r.Fr r.Dp r.Dv (r.Dp+r.Dv) r.Cwp r.Cwv r.Cw_tot
end
@printf "═══════════════════════════════════════════════════════════════════════════\n"

# Plot
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "containership_fr_resistance"))
mkpath(OUTDIR)
fig = Figure(size=(800, 480))
ax = Axis(fig[1, 1]; xlabel="Fr", ylabel="Cw",
    title="Wigley resistance curve (Re=5000, ν₀ MULES, no propeller)")
Frs = [r.Fr for r in results]
Cwp = [r.Cwp for r in results]
Cwv = [r.Cwv for r in results]
Cwt = [r.Cw_tot for r in results]
lines!(ax, Frs, Cwp; color=:tomato, linewidth=2, label="Cw_pressure")
lines!(ax, Frs, Cwv; color=:steelblue, linewidth=2, label="Cw_viscous")
lines!(ax, Frs, Cwt; color=:black, linewidth=2.5, linestyle=:dash, label="Cw_total")
axislegend(ax)
save(joinpath(OUTDIR, "Cw_vs_Fr.png"), fig)
println("Wrote $(joinpath(OUTDIR, "Cw_vs_Fr.png"))")
