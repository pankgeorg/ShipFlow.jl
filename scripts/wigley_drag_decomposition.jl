#!/usr/bin/env julia
#
# L2: drag decomposition. Same Wigley hull, two configurations:
#   (a) With VoF free-surface — total drag includes wave-making
#   (b) Deep water (α₀ = 1 everywhere, no free surface) — viscous + form only
# Difference (a − b) is the wave-resistance contribution.
#
# Caveat: the deep-water case still has a "phantom waterline" because
# the Wigley SDF is defined relative to body-z=0; with α=1 everywhere
# the hull behaves as a fully-submerged streamlined body, no surface
# piercing. So we're comparing surface-piercing vs fully-submerged of
# the same shape, not strict free-surface-on vs off.

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
using Printf, Statistics, CairoMakie

const NX = 128; const NY = 64; const NZ = 32
const NSTEPS = 120
const L_c = 36.0; const B_c = 8.0; const T_c = 5.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞ = 1.0
const Re = 5000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c

const FR_LIST = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45]

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L=L_c, B=B_c, T=T_c, map=hull_map)

function run_case(Fr, with_free_surface::Bool)
    G_c = U∞^2 / (Fr^2 * L_c)
    α₀_fn = with_free_surface ? (i, x) -> (x[3] ≤ H_w_c ? 1.0 : 0.0) :
                                 (i, x) -> 1.0   # all water = deep
    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀_fn, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float64)
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
    Dp = Float64[]; Dv = Float64[]
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1e-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        if step > NSTEPS - NSTEPS÷4
            Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
            push!(Dp, -Fp[1]); push!(Dv, -Fv[1])
        end
    end
    return mean(Dp), mean(Dv)
end

println("=== L2: drag decomposition (free surface vs deep) ===")
println("Fr     D_total  D_deep   D_wave (= D_total - D_deep)  Cw_ratio")
println("─" ^ 70)

results = []
for Fr in FR_LIST
    @printf "Fr=%.2f " Fr
    flush(stdout)
    Dp_fs, Dv_fs = run_case(Fr, true)
    @printf "(fs done) "
    flush(stdout)
    Dp_d,  Dv_d  = run_case(Fr, false)
    D_total = Dp_fs + Dv_fs
    D_deep  = Dp_d  + Dv_d
    D_wave  = D_total - D_deep
    ratio = D_wave / max(D_total, 1e-9)
    push!(results, (Fr=Fr, D_total=D_total, D_deep=D_deep, D_wave=D_wave, ratio=ratio))
    @printf "  D_t=%6.2f  D_d=%6.2f  D_w=%6.2f  wave/total=%.2f\n" D_total D_deep D_wave ratio
end

println("─" ^ 70)
println("Note: D_deep is the fully-submerged-hull drag (no free surface).")
println("      D_wave is the difference, attributed to surface wave-making.")

# Plot
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "drag_decomp"))
mkpath(OUTDIR)
Frs = [r.Fr for r in results]
Dts = [r.D_total for r in results]
Dds = [r.D_deep for r in results]
Dws = [r.D_wave for r in results]
fig = Figure(size=(900, 480))
ax = Axis(fig[1, 1]; xlabel="Fr", ylabel="Drag (cell units)",
    title="L2: Wigley drag decomposition (Re=5000)")
lines!(ax, Frs, Dts; color=:black, linewidth=2.5, label="D_total (free-surface)")
lines!(ax, Frs, Dds; color=:steelblue, linewidth=2, label="D_deep (no free surface)")
lines!(ax, Frs, Dws; color=:tomato, linewidth=2, linestyle=:dash, label="D_wave = D_total − D_deep")
hlines!(ax, [0]; color=:grey, linestyle=:dot)
axislegend(ax)
save(joinpath(OUTDIR, "decomp.png"), fig)
println("Wrote $(joinpath(OUTDIR, "decomp.png"))")
