#!/usr/bin/env julia
#
# Run the headline simulation (Wigley + AD + VoF + WALE) and dump
# midplane slices of α (free surface), u_x (velocity), as PNGs.
# Visualises the wave pattern and the propeller wake.

import Pkg
# Build a temp env with the ShipFlow stack + CairoMakie.
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Propellers.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)

using WaterLily
using VoF
using ShipShapes
using Turbulence
using Propellers
const SVector = Propellers.StaticArrays.SVector
using Printf
using CairoMakie

const NX = 128; const NY = 48; const NZ = 64
const L_c = 56; const B_c = 10; const T_c = 8
const ρ_w = 10.0; const ρ_a = 1.0
const U∞ = 1.0
const Fr = 0.25
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c
const T_NUM = Float64

const H_w_c   = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = NY/2
const prop_zc = H_w_c - T_c/2
const R_prop  = 1.5 * T_c/2
const W_prop  = 1.5
# Run at the self-propulsion C_T found earlier
const C_T = 2.27
const thrust = 0.5 * 1.0 * π * R_prop^2 * U∞^2 * C_T

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
turb = WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0))
disk = ActuatorDisk(
    center = SVector(prop_xc, prop_yc, prop_zc),
    axis   = SVector(1.0, 0.0, 0.0),
    R = R_prop, w = W_prop, thrust = thrust,
)

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end
sim = WaterLily.Simulation((NX, NY, NZ),
    (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
    T = T_NUM, ν = turb.ν,
    g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
    Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = T_NUM(U∞),
)
disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "100"))
@info "Running $NSTEPS steps for snapshot…"
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois;
        udf = disk_udf,
        pois_tol = T_NUM(1e-8), pois_itmx = 100)
    step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
    update_νt!(turb, sim.flow.u, vof.ν)
    if mod(step, 20) == 0
        @info "step=$step"
    end
end

# Side-view slice at midplane y = NY/2+1 (1-indexed)
j_mid = Int(round(NY/2)) + 1
α_slice = vof.α[2:NX+1, j_mid, 2:NZ+1]
u_slice = sim.flow.u[2:NX+1, j_mid, 2:NZ+1, 1]
println("α range: ", extrema(α_slice), "  u_x range: ", extrema(u_slice))

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_snapshot"))
mkpath(OUTDIR)

# --- α (free surface) heatmap ---
fig1 = Figure(size=(1000, 500))
ax1 = Axis(fig1[1, 1]; aspect = DataAspect(),
    xlabel="x [cells]", ylabel="z [cells]",
    title="α (water=1, air=0) — midplane y = $j_mid, step=$NSTEPS")
hm1 = heatmap!(ax1, 1:NX, 1:NZ, transpose(α_slice'), colormap=:RdBu, colorrange=(0,1))
Colorbar(fig1[1, 2], hm1)
# overlay waterline
hlines!(ax1, [H_w_c]; color=:black, linewidth=1, linestyle=:dash)
# overlay propeller location
scatter!(ax1, [prop_xc], [prop_zc]; marker=:diamond, color=:lime, markersize=12)
save(joinpath(OUTDIR, "alpha_midplane.png"), fig1)
println("Wrote ", joinpath(OUTDIR, "alpha_midplane.png"))

# --- u_x velocity heatmap ---
fig2 = Figure(size=(1000, 500))
ax2 = Axis(fig2[1, 1]; aspect = DataAspect(),
    xlabel="x [cells]", ylabel="z [cells]",
    title="u_x velocity — midplane y = $j_mid, step=$NSTEPS")
hm2 = heatmap!(ax2, 1:NX, 1:NZ, transpose(u_slice'), colormap=:viridis)
Colorbar(fig2[1, 2], hm2)
hlines!(ax2, [H_w_c]; color=:white, linewidth=1, linestyle=:dash)
scatter!(ax2, [prop_xc], [prop_zc]; marker=:diamond, color=:lime, markersize=12)
save(joinpath(OUTDIR, "u_midplane.png"), fig2)
println("Wrote ", joinpath(OUTDIR, "u_midplane.png"))

# --- Top view: alpha + u in horizontal plane just below waterline ---
k_below = Int(round(H_w_c - 1))   # 1 cell below waterline
u_top = sim.flow.u[2:NX+1, 2:NY+1, k_below, 1]
fig3 = Figure(size=(1000, 500))
ax3 = Axis(fig3[1, 1]; aspect = DataAspect(),
    xlabel="x [cells]", ylabel="y [cells]",
    title="u_x — horizontal slice at z=$k_below (just below waterline)")
hm3 = heatmap!(ax3, 1:NX, 1:NY, transpose(u_top'), colormap=:viridis)
Colorbar(fig3[1, 2], hm3)
save(joinpath(OUTDIR, "u_top.png"), fig3)
println("Wrote ", joinpath(OUTDIR, "u_top.png"))

println("\nSnapshot done.")
