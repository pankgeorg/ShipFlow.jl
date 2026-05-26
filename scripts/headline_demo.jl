#!/usr/bin/env julia
#
# F3: headline animation. Integrated VLM stack with every recent
# upgrade switched on:
#   - Containership hull (Cb ≈ 0.75, parallel midbody)
#   - WALE LES (Turbulence.jl)
#   - Sectional smear (smear_blades! — 3 blades × 4 sections)
#   - Two-way coupling (rudder samples flow.u via trilinear_inflow)
#
# Per frame, dumps a 2-panel figure:
#   top:   η elevation heatmap (free-surface)
#   side:  u_x velocity slice at y = NY/2 (the rotor race + hull wake)

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["CairoMakie", "VortexLattice"]; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, Turbulence, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "80"))
const BURNIN  = parse(Int, get(ENV, "WL_BURNIN", "15"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.30"))
const δ_DEG  = parse(Float32, get(ENV, "WL_DELTA", "10"))
const PAR_FRAC = parse(Float64, get(ENV, "WL_PAR_FRAC", "0.5"))

const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const G_c = U∞^2 / (FR^2 * L_c)
const Re = parse(Float32, get(ENV, "WL_RE", "20000"))
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
const J = 0.32
const Ω = Float64(π) * U∞ / (J * R_prop)

const rud_xc = Float32(prop_xc + R_prop * 0.5)   # inside the race
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
    par_frac = PAR_FRAC, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

@printf "=== F3 headline animation ===\n"
@printf "  Hull       = Containership (par_frac=%.2f, Cb≈%.2f)\n" PAR_FRAC (1+PAR_FRAC)/2
@printf "  Grid       = %d × %d × %d (%.2fM cells)\n" NX NY NZ NX*NY*NZ/1e6
@printf "  Frames     = %d  burn-in=%d\n" NFRAMES BURNIN
@printf "  Fr=%.2f, Re=%.0f, J=%.2f, δ=%.1f°\n" FR Re J δ_DEG
@printf "  features   = WALE LES + sectional smear + two-way rudder\n"
flush(stdout)

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
turb = WALE((NX, NY, NZ); Cw = 0.5f0, ν₀ = 0f0)
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
    T = Float32, ν = turb.ν,
    g = (i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)

# Cache rotor thrust+torque once (no two-way on rotor in this demo).
r_rot = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop
@printf "  rotor    |CT|=%.3f thrust=%.2f |CQ|=%.3f torque=%.2f\n" abs(r_rot.CT) thrust abs(r_rot.CQ) torque

function vlm_udf(flow, t; kwargs...)
    # G4: sectional smear for the rotor body force.
    smear_blades!(flow.f, thrust, torque,
                  SVector(prop_xc, prop_yc, prop_zc),
                  SVector(1.0, 0.0, 0.0),
                  Float64(R_prop), 0.2 * Float64(R_prop);
                  N_blades=3, N_sections=4, ε=1.5)
    # F2: two-way coupling — rudder samples flow.u.
    inflow_fn = trilinear_inflow(flow.u)
    r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞; inflow=inflow_fn)
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
        # mask out hull silhouette in plan view
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

η  = fill(NaN32, size(vof.α, 1), size(vof.α, 2))

# Burn-in (with WALE update so ν_t spins up).
@info @sprintf("Burn-in %d steps…", BURNIN)
t0 = time()
for s in 1:BURNIN
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    update_νt!(turb, sim.flow.u, vof.ν)
end
@printf "  burnin %d done in %.1fs\n" BURNIN time()-t0

# Colour-scale once on the post-burnin state.
eta_surface!(η, vof.α)
ηfin = filter(isfinite, vec(η))
η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
ηmax = isempty(ηfin) ? 0.3f0 : Float32(quantile(abs.(ηfin .- η_med), 0.95)) * 2
ηmax = max(ηmax, 0.3f0)
u_range = (-1.5f0, 3.0f0)
@printf "  η_median=%.3f ηmax=%.3f u_range=%.1f..%.1f\n" η_med ηmax u_range[1] u_range[2]

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "headline_demo", "frames"))
mkpath(OUTDIR)

@info @sprintf("Rendering %d frames…", NFRAMES)
t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    update_νt!(turb, sim.flow.u, vof.ν)
    eta_surface!(η, vof.α)

    # Side view: u_x at j = NY÷2 (centreline)
    j_mid = NY ÷ 2
    ux_side = @view sim.flow.u[:, j_mid, :, 1]

    fig = Figure(size=(1200, 720))
    ax_top = Axis(fig[1, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="y (cells)",
        title=@sprintf("Headline demo  frame %3d/%d  (Containership + WALE + sectional + two-way)",
                       frame, NFRAMES))
    hm_top = heatmap!(ax_top, 1:size(η,1), 1:size(η,2), η .- η_med;
        colormap=:RdBu, colorrange=(-ηmax, ηmax),
        nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
    Colorbar(fig[1, 2], hm_top, label="η - ⟨η⟩")
    poly!(ax_top, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
                   Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
                   Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
                   Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
        color=:transparent, strokecolor=:black, strokewidth=1.5)
    scatter!(ax_top, [Point2f(prop_xc, prop_yc)];
        color=:orange, marker=:cross, markersize=16, strokewidth=2)
    scatter!(ax_top, [Point2f(rud_xc, rud_yc)];
        color=:lime, marker=:dtriangle, markersize=14, strokewidth=2)

    ax_side = Axis(fig[2, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="z (cells)",
        title="centreline u_x  (side view)")
    hm_side = heatmap!(ax_side, 1:size(ux_side,1), 1:size(ux_side,2), ux_side;
        colormap=:viridis, colorrange=u_range)
    Colorbar(fig[2, 2], hm_side, label="u_x")
    hlines!(ax_side, [H_w_c]; color=:cyan, linestyle=:dash)  # nominal waterline
    lines!(ax_side, [hull_xc - L_c/2, hull_xc + L_c/2,
                     hull_xc + L_c/2, hull_xc - L_c/2,
                     hull_xc - L_c/2],
                    [hull_zc - T_c, hull_zc - T_c,
                     hull_zc, hull_zc, hull_zc - T_c];
        color=:black, linewidth=1.5)
    scatter!(ax_side, [Point2f(prop_xc, prop_zc)];
        color=:orange, marker=:cross, markersize=16, strokewidth=2)
    scatter!(ax_side, [Point2f(rud_xc, rud_zc)];
        color=:lime, marker=:dtriangle, markersize=14, strokewidth=2)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 10 == 0
        @printf "frame %3d / %d  |u|=%.3f  elapsed=%.1fs\n" frame NFRAMES maximum(abs, sim.flow.u) time()-t0
        flush(stdout)
    end
end
@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
