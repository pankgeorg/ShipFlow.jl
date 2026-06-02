#!/usr/bin/env julia
#
# Frame-dump variant of wigley_full_stack_VLM.jl. Same simulation;
# saves a top-down η heatmap per frame so you can assemble a GIF
# of the integrated VLM full stack evolving in time.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["CairoMakie", "VortexLattice"]; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "60"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.40"))
const δ_DEG  = parse(Float32, get(ENV, "WL_DELTA", "10"))

const L_c = 32f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const G_c = U∞^2 / (FR^2 * L_c)
const Re = parse(Float32, get(ENV, "WL_RE", "5000"))
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5)
const hull_yc = Float32(NY/2)
const hull_zc = H_w_c

const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J = 0.7
const Ω = Float64(π) * U∞ / (J * R_prop)

const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

@printf "=== Integrated VLM frames ===\n"
@printf "  Grid       = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Frames     = %d  (≈ %.1fs @ 30 fps)\n" NFRAMES NFRAMES/30
@printf "  Fr=%.2f, Re=%.0f, δ=%.1f°\n\n" FR Re δ_DEG
flush(stdout)

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)

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

function vlm_udf(flow, t; kwargs...)
    r_rot = rotor_forces(rotor, U∞, Ω)
    thrust = abs(r_rot.CT * 0.5 * U∞^2 * π * R_prop^2)
    smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞)
    side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                 SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
    return nothing
end

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

η = fill(NaN32, size(vof.α, 1), size(vof.α, 2))

# Burn-in
@info "Burn-in 15 steps…"
for _ in 1:15
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
end

# Colour-scale once
eta_surface!(η, vof.α)
ηfin = filter(isfinite, vec(η))
η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
ηmax = isempty(ηfin) ? 1f0 : Float32(quantile(abs.(ηfin .- η_med), 0.95)) * 2
ηmax = max(ηmax, 0.3f0)
@printf "η_median=%.3f, ηmax=%.3f\n" η_med ηmax

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "vlm_full", "frames"))
mkpath(OUTDIR)

@info "Rendering $NFRAMES frames…"
t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    eta_surface!(η, vof.α)

    fig = Figure(size=(1200, 540))
    ax = Axis(fig[1, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="y (cells)",
        title=@sprintf("Wigley + VLM rotor + VLM rudder (δ=%.1f°), frame %3d", δ_DEG, frame))
    hm = heatmap!(ax, 1:size(η,1), 1:size(η,2), η .- η_med;
        colormap=:RdBu, colorrange=(-ηmax, ηmax),
        nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
    Colorbar(fig[1, 2], hm, label="η - ⟨η⟩")
    poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
               Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
               Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
        color=:transparent, strokecolor=:black, strokewidth=1.5)
    scatter!(ax, [Point2f(prop_xc, prop_yc)];
        color=:orange, marker=:cross, markersize=16, strokewidth=2)
    scatter!(ax, [Point2f(rud_xc, rud_yc)];
        color=:lime, marker=:dtriangle, markersize=14, strokewidth=2)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 10 == 0
        @printf "frame %3d / %d (elapsed %.1fs)\n" frame NFRAMES time()-t0
        flush(stdout)
    end
end
@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
