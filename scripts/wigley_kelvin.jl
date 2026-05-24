#!/usr/bin/env julia
#
# Dedicated wide+long domain to actually capture the Kelvin wave wedge
# behind a Wigley hull.  Earlier 96×48×48 snapshots showed a wake but
# no clear V-pattern — the spanwise width is too small (NY=48 vs hull
# length L=48) to fit the 19.47° Kelvin half-angle wedge.
#
# Setup: NX=256 (≥4 hull lengths downstream), NY=128 (~2 hull lengths
# wide), NZ=48 (vertical).  Float32 throughout for speed.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)

using WaterLily, VoF, ShipShapes
using ShipShapes: StaticArrays
using Turbulence
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "256"))
const NY = parse(Int, get(ENV, "WL_NY", "128"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))

const L_c = 48f0     # hull length (cells)
const B_c = 8f0      # narrower beam so wave pattern dominates over wake
const T_c = 6f0
const ρ_w = 10f0
const ρ_a = 1f0
const U∞  = 1f0
const Fr  = parse(Float32, get(ENV, "WL_FR", "0.30"))
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = parse(Float32, get(ENV, "WL_RE", "10000"))
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c

# Hull at 1/5 along x — plenty of downstream room
const H_w_c   = Float32(NZ/2)
const hull_xc = Float32(NX/5)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

@printf "=== Kelvin pattern domain ===\n"
@printf "  Grid          = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Hull          = L=%.0f, B=%.0f, T=%.0f at (%.0f, %.0f)\n" L_c B_c T_c hull_xc hull_yc
@printf "  Downstream    = %.1f hull lengths\n" (NX - hull_xc)/L_c
@printf "  Spanwise      = %.1f hull lengths\n" Float32(NY)/L_c
@printf "  Fr=%.2f  Re=%.0f  ρ_w/ρ_a=%.0f\n\n" Fr Re ρ_w/ρ_a

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
const CW = parse(Float32, get(ENV, "WL_CW", "0.5"))
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

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "200"))
@info "Running $NSTEPS steps on $(NX*NY*NZ) cells…"
t0 = time()
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois;
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof!(vof, sim; dt = sim.flow.Δt[end-1],
        mass_repair = parse(Bool, get(ENV, "WL_MASS_REPAIR", "true")))
    if CW > 0
        update_νt!(turb, sim.flow.u, vof.ν)
    else
        # Cw=0 → tracking only the per-cell molecular ν from VoFFlow.
        copyto!(turb.ν, vof.ν)
    end
    if mod(step, 25) == 0
        u_max = maximum(abs, sim.flow.u)
        @info @sprintf("step=%3d  Δt=%.3f  |u|=%.3f  elapsed=%.1fs",
            step, sim.flow.Δt[end-1], u_max, time()-t0)
    end
end

# Top view: free-surface elevation η(x,y). For each (i,j) find the
# highest k with α>0.5 → η = z_surf - H_w_c.
function eta_surface(α)
    nx, ny, nz = size(α)
    η = fill(NaN32, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        # Skip cells whose centre is inside the hull box (so the
        # hull's footprint doesn't contaminate the surface map).
        x_hull = Float32(i - 1.5) - hull_xc
        y_hull = Float32(j - 1.5) - hull_yc
        if abs(x_hull) ≤ L_c/2 && abs(y_hull) ≤ B_c/2
            continue
        end
        for k in nz:-1:1
            if α[i, j, k] > 0.5f0
                η[i, j] = Float32(k - 1.5) - H_w_c
                break
            end
        end
    end
    return η
end

η = eta_surface(vof.α)
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_kelvin"))
mkpath(OUTDIR)

# Use a divergent colormap centred on 0 to highlight peaks/troughs
fig = Figure(size=(1200, 600))
ax = Axis(fig[1, 1]; aspect = DataAspect(),
    xlabel="x (cells, +x = downstream)",
    ylabel="y (cells)",
    title=@sprintf("Free-surface elevation η — Wigley hull at Fr=%.2f, Re=%.0f, step=%d", Fr, Re, NSTEPS))
# Subtract the median (= "the baseline VoF mass-loss drift") so what
# we render is the WAVE pattern, not the global surface drop.
ηfinite = filter(isfinite, vec(η))
η_med = isempty(ηfinite) ? 0f0 : Float32(quantile(ηfinite, 0.5))
η_disp = η .- η_med
@printf("subtracted baseline η_median = %.3f\n", η_med)
# Clip color range to twice the 90th-percentile of |η_displayed|.
ηd_finite = filter(isfinite, vec(η_disp))
ηmax = isempty(ηd_finite) ? 1f0 : Float32(quantile(abs.(ηd_finite), 0.90)) * 2
ηmax = max(ηmax, 0.5f0)
hm = heatmap!(ax, 1:NX, 1:NY, η_disp; colormap=:RdBu, colorrange=(-ηmax, ηmax),
    nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
Colorbar(fig[1, 2], hm, label="η (cells from waterline)")
# Annotate hull location
poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
           Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
    color=:transparent, strokecolor=:black, strokewidth=1)
# Kelvin half-angle reference lines (19.47° from hull stern)
ang = deg2rad(19.47)
x_stern = hull_xc + L_c/2
xs = collect(Float32(0):Float32(NX) - x_stern)
ys_plus  = hull_yc .+ xs .* Float32(tan(ang))
ys_minus = hull_yc .- xs .* Float32(tan(ang))
lines!(ax, xs .+ x_stern, ys_plus;  color=:lime, linestyle=:dash, linewidth=1.5)
lines!(ax, xs .+ x_stern, ys_minus; color=:lime, linestyle=:dash, linewidth=1.5)

save(joinpath(OUTDIR, "eta_top.png"), fig)
println("\nWrote $(joinpath(OUTDIR, "eta_top.png"))")
ηfin = filter(isfinite, vec(η)); @printf("η range (excl hull): min=%.3f, max=%.3f, p1=%.3f, p99=%.3f, std=%.3f\n", minimum(ηfin), maximum(ηfin), quantile(ηfin, 0.01), quantile(ηfin, 0.99), Statistics.std(ηfin))
