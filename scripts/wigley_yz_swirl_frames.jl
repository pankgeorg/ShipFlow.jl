#!/usr/bin/env julia
#
# Y-Z cross-section frame dump: show the swirling wake of the
# SwirlingDisk by rendering u_θ on a y-z plane downstream of the
# propeller, every frame.

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
const R_prop  = Float32(T_c / 2 * 0.6)
const prop_xc = Float32(hull_xc + L_c/2 + T_c/2)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const C_T     = parse(Float32, get(ENV, "WL_CT", "0.6"))
const thrust  = 0.5f0 * π * R_prop^2 * U∞^2 * C_T
const torque  = parse(Float32, get(ENV, "WL_TORQUE_RATIO", "0.5")) * thrust * R_prop

# Slice at 2R behind disk
const X_SLICE = round(Int, prop_xc + 2*R_prop)

@printf "=== Y-Z swirl cross-section ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Disk       = R=%.1f at (%.1f, %.1f, %.1f)\n" R_prop prop_xc prop_yc prop_zc
@printf "  Slice plane= x = %d  (= prop_xc + 2R)\n" X_SLICE
@printf "  Frames     = %d\n\n" NFRAMES
flush(stdout)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)

disk = SwirlingDisk(
    center = SVector(prop_xc, prop_yc, prop_zc),
    axis   = SVector(1f0, 0f0, 0f0),
    R      = R_prop, w = 1.5f0,
    thrust = thrust, torque = torque,
)
disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)

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

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_yz_swirl", "frames"))
mkpath(OUTDIR)

@info "Burn-in 15 steps…"; flush(stdout)
for _ in 1:15
    WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
end

@info "Rendering $NFRAMES y-z slices at x=$X_SLICE …"; flush(stdout)
t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    # Slice at i=X_SLICE: sample u_θ over (j, k) — j is y, k is z.
    NX_α, NY_α, NZ_α = size(vof.α)
    uθ_slice = Array{Float32}(undef, NY_α, NZ_α)
    α_slice  = Array{Float32}(undef, NY_α, NZ_α)
    @inbounds for j in 1:NY_α, k in 1:NZ_α
        y = Float32(j - 1.5) - prop_yc
        z = Float32(k - 1.5) - prop_zc
        r = sqrt(y*y + z*z)
        α_slice[j, k] = vof.α[X_SLICE, j, k]
        if r < 1f-3 || α_slice[j, k] < 0.5f0
            uθ_slice[j, k] = NaN32
            continue
        end
        # Interpolate u_y, u_z to cell centres
        uy = 0.5f0 * (sim.flow.u[X_SLICE, j, k, 2] + sim.flow.u[X_SLICE, min(j+1,NY_α), k, 2])
        uz = 0.5f0 * (sim.flow.u[X_SLICE, j, k, 3] + sim.flow.u[X_SLICE, j, min(k+1,NZ_α), 3])
        # Tangential: ê_θ = axis × r̂, with axis=+x, r̂=(0, y/r, z/r) ⇒ ê_θ=(0, -z/r, y/r)
        uθ_slice[j, k] = (-z * uy + y * uz) / r
    end

    # Symmetric colour range
    finite_uθ = filter(isfinite, vec(uθ_slice))
    uθmax = isempty(finite_uθ) ? 0.3f0 : max(0.05f0, Float32(quantile(abs.(finite_uθ), 0.98)))

    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1]; aspect = DataAspect(),
        xlabel="y (cells)", ylabel="z (cells, +z = up)",
        title=@sprintf("Y-Z slice at x=%d  (2R behind disk)  —  frame %3d  Fr=%.2f", X_SLICE, frame, FR))
    hm = heatmap!(ax, 1:size(uθ_slice,1), 1:size(uθ_slice,2), uθ_slice;
        colormap=:RdBu, colorrange=(-uθmax, uθmax),
        nan_color=RGBAf(0.85, 0.92, 1.0, 1.0))
    Colorbar(fig[1, 2], hm, label="u_θ  (right-hand swirl about +x is positive)")
    # Disk outline (projected): a circle of radius R_prop centered at (prop_yc, prop_zc)
    θ = range(0f0, 2f0π, length=64)
    disk_y = prop_yc .+ R_prop .* cos.(θ)
    disk_z = prop_zc .+ R_prop .* sin.(θ)
    lines!(ax, disk_y, disk_z; color=:black, linestyle=:dash, linewidth=1.5, label="disk perimeter")
    # Waterline reference
    hlines!(ax, [H_w_c]; color=:cyan, linestyle=:dot, linewidth=1.5)

    fname = joinpath(OUTDIR, @sprintf("frame_%05d.png", frame))
    save(fname, fig)
    if frame % 10 == 0
        @printf "frame %3d / %d (elapsed %.1fs)\n" frame NFRAMES time()-t0
        flush(stdout)
    end
end
@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
