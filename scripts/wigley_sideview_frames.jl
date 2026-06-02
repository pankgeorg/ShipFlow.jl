#!/usr/bin/env julia
#
# Side-view frame dump: vertical x-z slice through y=hull_yc.
# Shows axial velocity u_x as heatmap, with the free-surface
# (α=0.5 isocontour) and hull cross-section overlaid.
# WL_PROP=1 WL_SWIRL=1 enable the propeller.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Propellers.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes, Turbulence, Propellers
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

println("packages loaded"); flush(stdout)

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "96"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "120"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.40"))

const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const G_c = U∞^2 / (FR^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

@printf "=== Side-view frame dump ===\n"
@printf "  Grid       = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Frames     = %d  (≈ %.1fs @ 30 fps)\n" NFRAMES NFRAMES/30
flush(stdout)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)

const HAS_PROP  = get(ENV, "WL_PROP",  "0") == "1"
const HAS_SWIRL = get(ENV, "WL_SWIRL", "0") == "1"
const R_prop  = Float32(T_c / 2 * 0.6)
const prop_xc = Float32(hull_xc + L_c/2 + T_c/2)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const C_T     = parse(Float32, get(ENV, "WL_CT", "0.6"))
const thrust  = 0.5f0 * π * R_prop^2 * U∞^2 * C_T
const torque  = parse(Float32, get(ENV, "WL_TORQUE_RATIO", "0.5")) * thrust * R_prop
const disk = HAS_SWIRL ?
    SwirlingDisk(center=SVector(prop_xc, prop_yc, prop_zc),
                 axis=SVector(1f0, 0f0, 0f0), R=R_prop, w=1.5f0,
                 thrust=thrust, torque=torque) :
    ActuatorDisk(center=SVector(prop_xc, prop_yc, prop_zc),
                 axis=SVector(1f0, 0f0, 0f0), R=R_prop, w=1.5f0,
                 thrust=thrust)
disk_udf = HAS_PROP ? ((flow, t; kwargs...) -> disk(flow, t; kwargs...)) : nothing
if HAS_PROP
    label = HAS_SWIRL ? "SwirlingDisk" : "ActuatorDisk"
    @printf "  Propeller  = %s, R=%.2f at (%.1f, %.1f, %.1f)\n" label R_prop prop_xc prop_yc prop_zc
end

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end

sim = WaterLily.Simulation((NX, NY, NZ),
    (U∞, 0f0, 0f0), L_c;
    T = Float32, ν = VoF.viscosity(vof),
    g = (i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)

# Slice index for y=hull_yc
const j_slice = round(Int, hull_yc + 1.5)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_sideview", "frames"))
mkpath(OUTDIR)

@info "Burn-in 15 steps…"; flush(stdout)
for _ in 1:15
    if HAS_PROP
        WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    else
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    end
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
end

# Wigley parabolic profile at y=0 (centerplane): z_keel(x) = T (within hull length)
# Actually Wigley B*S(x,z) with S=(1-(2x/L)²)(1-(z/T)²); the surface is at S=1/2 for example.
# For visualization just sketch hull as a rectangle at the centerplane.
@info "Rendering $NFRAMES frames…"; flush(stdout)
t0 = time()
for frame in 1:NFRAMES
    if HAS_PROP
        WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    else
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    end
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    # x-z slice at y=hull_yc: take u_x (component 1).
    # Mask air cells (α<0.5) to NaN so the heatmap shows water-only velocity.
    NX_α, NY_α, NZ_α = size(vof.α)
    ux_slice = Array{Float32}(undef, NX_α, NZ_α)
    α_slice  = Array{Float32}(undef, NX_α, NZ_α)
    @inbounds for i in 1:NX_α, k in 1:NZ_α
        # u_x lives on the x-face — interpolate to cell centre
        u_here = 0.5f0 * (sim.flow.u[i, j_slice, k, 1] + sim.flow.u[min(i+1, NX_α), j_slice, k, 1])
        α_slice[i, k]  = vof.α[i, j_slice, k]
        ux_slice[i, k] = α_slice[i, k] >= 0.5f0 ? u_here : NaN32
    end

    fig = Figure(size=(1100, 500))
    ax = Axis(fig[1, 1]; aspect = DataAspect(),
        xlabel="x (cells, +x = downstream)", ylabel="z (cells, +z = up)",
        title=@sprintf("Side view at y=%.1f (centerline)  —  frame %3d  Fr=%.2f", hull_yc, frame, FR))
    # u_x heatmap
    hm = heatmap!(ax, 1:size(ux_slice,1), 1:size(ux_slice,2), ux_slice;
        colormap=:roma, colorrange=(0f0, 2f0),
        nan_color=RGBAf(0.85, 0.92, 1.0, 1.0))   # light sky-blue for air
    Colorbar(fig[1, 2], hm, label="u_x / U∞")
    # Free surface contour at α=0.5
    contour!(ax, 1:size(α_slice,1), 1:size(α_slice,2), α_slice;
        levels=[0.5f0], color=:cyan, linewidth=2.5)
    # Wetted hull on the centerplane (y=0): the Wigley SDF reduces to the
    # solid rectangle |x-hull_xc| ≤ L/2, hull_zc-T ≤ z ≤ hull_zc.
    poly!(ax, [Point2f(hull_xc - L_c/2, hull_zc - T_c),
               Point2f(hull_xc + L_c/2, hull_zc - T_c),
               Point2f(hull_xc + L_c/2, hull_zc),
               Point2f(hull_xc - L_c/2, hull_zc)];
        color=(:grey20, 0.7), strokecolor=:black, strokewidth=1.5)
    # Cosmetic above-water freeboard / deck (NOT in the simulation; the
    # Wigley model is wetted-only). Tall enough to clear wake amplitude.
    fb = 1.5f0 * T_c
    poly!(ax, [Point2f(hull_xc - L_c*0.48, hull_zc),
               Point2f(hull_xc + L_c*0.48, hull_zc),
               Point2f(hull_xc + L_c*0.42, hull_zc + fb),
               Point2f(hull_xc - L_c*0.42, hull_zc + fb)];
        color=(:grey55, 0.95), strokecolor=:black, strokewidth=1.2)
    # Small deckhouse for ship-recognition
    poly!(ax, [Point2f(hull_xc - L_c*0.15, hull_zc + fb),
               Point2f(hull_xc + L_c*0.10, hull_zc + fb),
               Point2f(hull_xc + L_c*0.08, hull_zc + fb + T_c*0.6),
               Point2f(hull_xc - L_c*0.13, hull_zc + fb + T_c*0.6)];
        color=(:grey70, 0.95), strokecolor=:black, strokewidth=1.0)
    # Still-water reference line
    hlines!(ax, [hull_zc]; color=:grey50, linestyle=:dot, linewidth=1)
    # Disk marker
    if HAS_PROP
        scatter!(ax, [Point2f(prop_xc, prop_zc)];
            color=:orange, marker=:cross, markersize=14, strokewidth=2)
    end

    fname = joinpath(OUTDIR, @sprintf("frame_%05d.png", frame))
    save(fname, fig)
    if frame % 10 == 0
        @printf "frame %3d / %d (elapsed %.1fs)\n" frame NFRAMES time()-t0
        flush(stdout)
    end
end
@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
