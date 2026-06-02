#!/usr/bin/env julia
#
# Dump per-frame PNGs of the Kelvin wave-pattern evolution.
# Use ffmpeg externally to assemble: e.g.
#   ffmpeg -framerate 30 -i frame_%05d.png -c:v libx264 -pix_fmt yuv420p out.mp4

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
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "150"))   # 5 s @ 30 fps
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

@printf "=== Kelvin animation frames ===\n"
@printf "  Grid          = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Frames        = %d  (≈ %.1fs @ 30 fps)\n" NFRAMES NFRAMES/30
flush(stdout)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)

# Optional actuator-disk propeller behind the stern (WL_PROP=1)
const HAS_PROP = get(ENV, "WL_PROP", "0") == "1"
const R_prop  = Float32(T_c / 2 * 0.6)
const prop_xc = Float32(hull_xc + L_c/2 + T_c/2)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const C_T     = parse(Float32, get(ENV, "WL_CT", "0.6"))
const thrust  = 0.5f0 * π * R_prop^2 * U∞^2 * C_T
const HAS_SWIRL = get(ENV, "WL_SWIRL", "0") == "1"
const torque  = parse(Float32, get(ENV, "WL_TORQUE_RATIO", "0.5")) * thrust * R_prop
const disk = HAS_SWIRL ?
    SwirlingDisk(
        center = SVector(prop_xc, prop_yc, prop_zc),
        axis   = SVector(1f0, 0f0, 0f0),
        R      = R_prop, w = 1.5f0,
        thrust = thrust, torque = torque,
    ) :
    ActuatorDisk(
        center = SVector(prop_xc, prop_yc, prop_zc),
        axis   = SVector(1f0, 0f0, 0f0),
        R      = R_prop, w = 1.5f0,
        thrust = thrust,
    )
disk_udf = HAS_PROP ? ((flow, t; kwargs...) -> disk(flow, t; kwargs...)) : nothing
if HAS_PROP
    if HAS_SWIRL
        @printf "  Propeller     = R=%.2f at (%.1f, %.1f, %.1f), C_T=%.2f, thrust=%.3f, torque=%.3f (SWIRL)\n" R_prop prop_xc prop_yc prop_zc C_T thrust torque
    else
        @printf "  Propeller     = R=%.2f at (%.1f, %.1f, %.1f), C_T=%.2f, thrust=%.3f\n" R_prop prop_xc prop_yc prop_zc C_T thrust
    end
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

function eta_surface!(η, α)
    nx, ny, nz = size(α)
    fill!(η, NaN32)
    @inbounds for j in 1:ny, i in 1:nx
        x_hull = Float32(i - 1.5) - hull_xc
        y_hull = Float32(j - 1.5) - hull_yc
        if abs(x_hull) ≤ L_c/2 && abs(y_hull) ≤ B_c/2
            continue
        end
        if i < hull_xc - L_c
            continue
        end
        for k in (nz-1):-1:2
            if α[i, j, k] >= 0.5f0 && α[i, j, k+1] < 0.5f0
                Δα = α[i, j, k] - α[i, j, k+1]
                t = abs(Δα) > 1f-9 ? (0.5f0 - α[i, j, k+1]) / Δα : 0.5f0
                z = Float32(k+1 - 1.5) - t
                η[i, j] = z - H_w_c
                break
            end
        end
    end
    return η
end

η = fill(NaN32, size(vof.α, 1), size(vof.α, 2))

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_kelvin", "frames"))
mkpath(OUTDIR)

# Burn-in
@info "Burn-in 15 steps…"; flush(stdout)
for _ in 1:15
    if HAS_PROP
        WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    else
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    end
    if get(ENV, "WL_MULES", "0") == "1"
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
    else
        step_vof!(vof, sim; dt=sim.flow.Δt[end-1])
    end
end

# Set color scale from the post-burn-in state
eta_surface!(η, vof.α)
ηfin = filter(isfinite, vec(η))
η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
ηmax = isempty(ηfin) ? 1f0 : Float32(quantile(abs.(ηfin .- η_med), 0.95)) * 2
ηmax = max(ηmax, 0.3f0)
@printf "Colour range: η_median=%.3f, ηmax=%.3f\n" η_med ηmax
flush(stdout)

ang = deg2rad(19.47)
x_stern = hull_xc + L_c/2

t0 = time()
for frame in 1:NFRAMES
    if HAS_PROP
        WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf, pois_tol=1f-6, pois_itmx=50)
    else
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    end
    if get(ENV, "WL_MULES", "0") == "1"
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
    else
        step_vof!(vof, sim; dt=sim.flow.Δt[end-1])
    end
    eta_surface!(η, vof.α)
    η_disp = η .- η_med

    fig = Figure(size=(900, 500))
    ax = Axis(fig[1, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="y (cells)",
        title=@sprintf("frame %3d  step=%d  Fr=%.2f", frame, frame+15, FR))
    hm = heatmap!(ax, 1:size(η, 1), 1:size(η, 2), η_disp;
        colormap=:RdBu, colorrange=(-ηmax, ηmax),
        nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
    Colorbar(fig[1, 2], hm)
    poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
               Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
        color=:transparent, strokecolor=:black, strokewidth=1)
    if HAS_PROP
        scatter!(ax, [Point2f(prop_xc, prop_yc)];
            color=:orange, marker=:cross, markersize=14, strokewidth=1.5)
    end
    xs = collect(Float32(0):Float32(NX) - x_stern)
    lines!(ax, xs .+ x_stern, hull_yc .+ xs .* Float32(tan(ang));
        color=:lime, linestyle=:dash, linewidth=1)
    lines!(ax, xs .+ x_stern, hull_yc .- xs .* Float32(tan(ang));
        color=:lime, linestyle=:dash, linewidth=1)

    fname = joinpath(OUTDIR, @sprintf("frame_%05d.png", frame))
    save(fname, fig)
    if frame % 10 == 0
        @printf "frame %3d / %d (elapsed %.1fs)\n" frame NFRAMES time()-t0
        flush(stdout)
    end
end

@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
@printf "To assemble: ffmpeg -framerate 30 -i %s/frame_%%05d.png -c:v libx264 -pix_fmt yuv420p out.mp4\n" OUTDIR
