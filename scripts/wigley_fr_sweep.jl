#!/usr/bin/env julia
#
# Fr sweep: run 4 Wigley + free-surface (no propeller) cases at
# different Froude numbers and composite a 2×2 panel showing the
# Kelvin wake at each. Classical-CFD validation: the wedge half-angle
# should be Fr-invariant at the Kelvin value 19.47°.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "96"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "120"))

const FR_LIST = (0.20f0, 0.30f0, 0.40f0, 0.50f0)

const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

function eta_surface!(η, α)
    nx, ny, nz = size(α)
    fill!(η, NaN32)
    @inbounds for j in 1:ny, i in 1:nx
        x_hull = Float32(i - 1.5) - hull_xc
        y_hull = Float32(j - 1.5) - hull_yc
        if abs(x_hull) ≤ L_c/2 && abs(y_hull) ≤ B_c/2; continue; end
        if i < hull_xc - L_c; continue; end
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

@printf "=== Wigley Kelvin Fr sweep ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Hull       = L=%.0f, B=%.0f, T=%.0f at xc=%.0f\n" L_c B_c T_c hull_xc
@printf "  Fr list    = %s\n"  FR_LIST
@printf "  Steps each = %d\n\n" NSTEPS
flush(stdout)

snapshots = Vector{Matrix{Float32}}()
medians   = Float32[]

for FR in FR_LIST
    G_c = U∞^2 / (FR^2 * L_c)
    @info "Running Fr=$FR  (G_c=$(round(G_c, digits=4))) …"

    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)

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

    t0 = time()
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
        if step % 30 == 0
            @info @sprintf("  Fr=%.2f step=%3d/%d  |u|=%.3f  elapsed=%.1fs",
                FR, step, NSTEPS, maximum(abs, sim.flow.u), time()-t0)
        end
    end

    η = fill(NaN32, size(vof.α, 1), size(vof.α, 2))
    eta_surface!(η, vof.α)
    ηfin = filter(isfinite, vec(η))
    η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
    push!(snapshots, copy(η))
    push!(medians, η_med)
    @info @sprintf("  Fr=%.2f  η_median=%.3f  range=[%.3f, %.3f]",
        FR, η_med, isempty(ηfin) ? 0 : minimum(ηfin), isempty(ηfin) ? 0 : maximum(ηfin))
end

# Composite 2x2 figure
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_kelvin"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "fr_sweep.png")

# Common color scale from the 0.40 case (most active)
η40 = snapshots[3] .- medians[3]
ηfin40 = filter(isfinite, vec(η40))
ηmax = isempty(ηfin40) ? 1f0 : Float32(quantile(abs.(ηfin40), 0.95)) * 1.5f0
ηmax = max(ηmax, 0.3f0)

fig = Figure(size=(1500, 800))
ang = deg2rad(19.47)
x_stern = hull_xc + L_c/2

for (k, FR) in enumerate(FR_LIST)
    row = (k-1) ÷ 2 + 1
    col = (k-1) % 2 + 1
    ax = Axis(fig[row, col]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="y (cells)",
        title=@sprintf("Fr = %.2f  (Kelvin wedge half-angle = 19.47°)", FR))
    η_disp = snapshots[k] .- medians[k]
    hm = heatmap!(ax, 1:size(η_disp,1), 1:size(η_disp,2), η_disp;
        colormap=:RdBu, colorrange=(-ηmax, ηmax),
        nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
    poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
               Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
        color=:transparent, strokecolor=:black, strokewidth=1)
    xs = collect(Float32(0):Float32(NX) - x_stern)
    lines!(ax, xs .+ x_stern, hull_yc .+ xs .* Float32(tan(ang));
        color=:lime, linestyle=:dash, linewidth=1.5)
    lines!(ax, xs .+ x_stern, hull_yc .- xs .* Float32(tan(ang));
        color=:lime, linestyle=:dash, linewidth=1.5)
end
Label(fig[0, :], "Wigley hull — Kelvin wave pattern across Fr (wedge angle is Fr-invariant)";
    fontsize=18, halign=:center)
Colorbar(fig[1:2, 3], colormap=:RdBu, limits=(-ηmax, ηmax),
    label="η - ⟨η⟩  (cells)")

save(out, fig)
@printf "\nWrote %s\n" out
