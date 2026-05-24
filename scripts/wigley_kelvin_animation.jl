#!/usr/bin/env julia
#
# Animation: Kelvin wave pattern developing behind a Wigley hull.
# Same physics as wigley_kelvin.jl, but dumps a frame each step and
# assembles them into an MP4 via Makie's record().

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "ShipFlow.jl"))
# CairoMakie is added to an ephemeral overlay so we don't permanently
# pollute the ShipFlow Manifest.
let
    # Only add CairoMakie if it isn't already loadable.
    try
        @eval using CairoMakie
    catch
        Pkg.add("CairoMakie"; io=devnull)
    end
end

using WaterLily, VoF, ShipShapes, Turbulence
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

println("=== Animation script started, packages loaded ==="); flush(stdout)

const NX = parse(Int, get(ENV, "WL_NX", "256"))
const NY = parse(Int, get(ENV, "WL_NY", "128"))
const NZ = parse(Int, get(ENV, "WL_NZ", "64"))
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "300"))   # ≈10 s @ 30 fps
const FPS = parse(Int, get(ENV, "WL_FPS", "30"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.40"))
const CW = parse(Float32, get(ENV, "WL_CW", "0.0"))

const L_c = 48f0; const B_c = 8f0; const T_c = 6f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const G_c = U∞^2 / (FR^2 * L_c)
const Re = 10000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

@printf "=== Kelvin animation ===\n"
@printf "  Grid          = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Hull          = L=%.0f, B=%.0f, T=%.0f at (%.0f, %.0f)\n" L_c B_c T_c hull_xc hull_yc
@printf "  Frames        = %d @ %d fps  (%.1f s)\n" NFRAMES FPS NFRAMES/FPS
@printf "  Fr=%.2f  Re=%.0f  CW=%.2f\n\n" FR Re CW

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
turb = WALE((NX, NY, NZ); Cw = CW, ν₀ = 0f0, T = Float32)

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end

sim = WaterLily.Simulation((NX, NY, NZ),
    (U∞, 0f0, 0f0), L_c;
    T = Float32, ν = turb.ν,
    g = (i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)

# η_surface with sub-cell linear interpolation.
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

# Pre-allocate η buffer
η = fill(NaN32, size(vof.α, 1), size(vof.α, 2))

# Burn-in: a few steps to let the initial transient settle so the
# animation doesn't start from a flat field.
@info "Burn-in 20 steps…"
for _ in 1:20
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    if CW > 0
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
        update_νt!(turb, sim.flow.u, vof.ν)
    else
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
        copyto!(turb.ν, vof.ν)
    end
end

# Set up the figure once and update inside record().
fig = Figure(size=(1200, 600))
ang = deg2rad(19.47)
x_stern = hull_xc + L_c/2
xs = collect(Float32(0):Float32(NX) - x_stern)
ys_plus  = hull_yc .+ xs .* Float32(tan(ang))
ys_minus = hull_yc .- xs .* Float32(tan(ang))

# Initial frame for setting scales
eta_surface!(η, vof.α)
ηfin = filter(isfinite, vec(η))
η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
ηmax  = isempty(ηfin) ? 1f0 : Float32(quantile(abs.(ηfin .- η_med), 0.95)) * 2
ηmax  = max(ηmax, 0.3f0)
η_obs = Observable(η .- η_med)
title_obs = Observable("Wigley hull — Fr=$FR, Re=$(Int(Re)), step=20")

ax = Axis(fig[1, 1]; aspect = DataAspect(),
    xlabel = "x (cells, +x = downstream)",
    ylabel = "y (cells)",
    title = title_obs)
hm = heatmap!(ax, 1:size(η, 1), 1:size(η, 2), η_obs;
    colormap = :RdBu, colorrange = (-ηmax, ηmax),
    nan_color = :grey80, highclip = :darkblue, lowclip = :darkred)
Colorbar(fig[1, 2], hm, label="η (cells from waterline)")
poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
           Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
    color=:transparent, strokecolor=:black, strokewidth=1)
lines!(ax, xs .+ x_stern, ys_plus;  color=:lime, linestyle=:dash, linewidth=1.5)
lines!(ax, xs .+ x_stern, ys_minus; color=:lime, linestyle=:dash, linewidth=1.5)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_kelvin"))
mkpath(OUTDIR)
outfile = joinpath(OUTDIR, "kelvin_evolution.mp4")

@info "Recording $NFRAMES frames at $FPS fps → $outfile"
t0 = time()
record(fig, outfile, 1:NFRAMES; framerate=FPS) do step
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
    if CW > 0
        update_νt!(turb, sim.flow.u, vof.ν)
    else
        copyto!(turb.ν, vof.ν)
    end
    eta_surface!(η, vof.α)
    η_obs[] = η .- η_med
    title_obs[] = @sprintf("Wigley hull — Fr=%.2f, Re=%d, step=%d  (sim t=%.2f)",
        FR, Int(Re), step + 20, (step + 20) * sim.flow.Δt[end-1])
    if mod(step, 25) == 0
        @info @sprintf("frame %3d / %d  elapsed=%.1fs", step, NFRAMES, time() - t0)
    end
end
@info @sprintf("Done — wrote %s in %.1fs", outfile, time() - t0)
