#!/usr/bin/env julia
#
# Integrated headline demo: Wigley hull (BDIM) + BladedRotor (VLM) +
# Rudder (VLM) + VoF/MULES free surface + (no LES; keep moving parts
# minimal for the first integrated run).
#
# Each WaterLily step:
#   1. Solve VLM for the bladed rotor at J=0.7  → CT → axial thrust
#      smeared at the disk centre
#   2. Solve VLM for the rudder at the prescribed δ → CL → side-force
#      smeared at the rudder centre
#   3. mom_step! with the combined udf
#   4. step_vof_mules!
#
# Output: top-down η heatmap with both lifting surfaces marked.

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

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "96"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "100"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.40"))
const δ_DEG  = parse(Float32, get(ENV, "WL_DELTA", "10"))
const TWOWAY = get(ENV, "WL_TWOWAY", "0") == "1"

const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
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

# Rotor: just behind stern, depth half-draft
const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J = 0.7
const Ω = Float64(π) * U∞ / (J * R_prop)

# Rudder: behind the rotor, span below waterline
const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

@printf "=== Integrated VLM full stack ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Hull       = L=%.0f, B=%.0f, T=%.0f at xc=%.1f\n" L_c B_c T_c hull_xc
@printf "  Rotor      = 3×VLM at (%.1f,%.1f,%.1f), J=%.2f\n" prop_xc prop_yc prop_zc J
@printf "  Rudder     = AR=2 at (%.1f,%.1f,%.1f), δ=%.1f°\n" rud_xc rud_yc rud_zc δ_DEG
@printf "  Fr=%.2f, Re=%.0f, ρ_w/ρ_a=%.0f\n\n" FR Re ρ_w/ρ_a
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
    T = Float32, ν = vof.ν,
    g = (i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)

# Combined udf — both rotor and rudder per step.
function vlm_udf(flow, t; kwargs...)
    # Rotor: forward thrust at the rotor centre.
    r_rot = rotor_forces(rotor, U∞, Ω)
    thrust = abs(r_rot.CT * 0.5 * U∞^2 * π * R_prop^2)
    smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    # Rudder: side force at the rudder centre (CL acts in body-z; we
    # interpret as +y in ship-world via mapping). When TWOWAY=true,
    # the rudder's `additional_velocity` callback samples flow.u
    # at panel control points so VLM sees the actual local inflow
    # (rotor race + hull wake). Subtract V∞ since VortexLattice
    # treats `additional_velocity` as a perturbation on Vinf.
    inflow_cb = if TWOWAY
        u_sample = trilinear_inflow(flow.u)
        (xv) -> let v = u_sample(SVector(Float32(xv[1] + rud_xc),
                                         Float32(xv[2] + rud_yc),
                                         Float32(xv[3] + rud_zc)))
            SVector(v[1] - U∞, v[2], v[3])
        end
    else
        nothing
    end
    r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞; inflow=inflow_cb)
    side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    smear_force!(flow.f,
                 SVector(-Float32(drag), Float32(side), 0f0),
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

@info "Running $NSTEPS steps with hull + rotor (VLM) + rudder (VLM) + VoF…"
t0 = time()
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    if step % 20 == 0
        u_max = maximum(abs, sim.flow.u)
        @info @sprintf("step=%3d  |u|=%.3f  elapsed=%.1fs", step, u_max, time()-t0)
        flush(stdout)
    end
end

# Render top-down η with hull + rotor + rudder marked.
η = fill(NaN32, size(vof.α, 1), size(vof.α, 2))
eta_surface!(η, vof.α)
ηfin = filter(isfinite, vec(η))
η_med = isempty(ηfin) ? 0f0 : Float32(quantile(ηfin, 0.5))
ηmax = isempty(ηfin) ? 1f0 : Float32(quantile(abs.(ηfin .- η_med), 0.95)) * 2
ηmax = max(ηmax, 0.3f0)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "vlm_full"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "headline.png")

fig = Figure(size=(1500, 700))
ax = Axis(fig[1, 1]; aspect=DataAspect(),
    xlabel="x (cells)", ylabel="y (cells)",
    title=@sprintf("Wigley + VLM rotor + VLM rudder (δ=%.1f°), Fr=%.2f, %d steps", δ_DEG, FR, NSTEPS))
hm = heatmap!(ax, 1:size(η,1), 1:size(η,2), η .- η_med;
    colormap=:RdBu, colorrange=(-ηmax, ηmax),
    nan_color=:grey80, highclip=:darkblue, lowclip=:darkred)
Colorbar(fig[1, 2], hm, label="η - ⟨η⟩  (cells)")
poly!(ax, [Point2f(hull_xc - L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc - B_c/2),
           Point2f(hull_xc + L_c/2, hull_yc + B_c/2),
           Point2f(hull_xc - L_c/2, hull_yc + B_c/2)],
    color=:transparent, strokecolor=:black, strokewidth=1.5)
scatter!(ax, [Point2f(prop_xc, prop_yc)];
    color=:orange, marker=:cross, markersize=18, strokewidth=2)
scatter!(ax, [Point2f(rud_xc, rud_yc)];
    color=:lime, marker=:dtriangle, markersize=18, strokewidth=2)

save(out, fig)
@printf "\nWrote %s\n" out

# Also a side-view: u_x on the y = hull_yc x-z slice, water-only.
j_slice = round(Int, hull_yc + 1.5)
NX_α, NY_α, NZ_α = size(vof.α)
ux_slice = Array{Float32}(undef, NX_α, NZ_α)
@inbounds for i in 1:NX_α, k in 1:NZ_α
    αv = vof.α[i, j_slice, k]
    if αv < 0.5
        ux_slice[i, k] = NaN32
    else
        ux_slice[i, k] = 0.5f0 * (sim.flow.u[i, j_slice, k, 1] +
                                  sim.flow.u[min(i+1, NX_α), j_slice, k, 1])
    end
end

out2 = joinpath(OUTDIR, "headline_side.png")
fig2 = Figure(size=(1500, 600))
ax2 = Axis(fig2[1,1]; aspect=DataAspect(),
    xlabel="x (cells)", ylabel="z (cells)",
    title=@sprintf("Side view (y=%.1f): u_x/U∞,  δ=%.1f°, %d steps", hull_yc, δ_DEG, NSTEPS))
hm2 = heatmap!(ax2, 1:size(ux_slice,1), 1:size(ux_slice,2), ux_slice;
    colormap=:roma, colorrange=(0f0, 2f0),
    nan_color=RGBAf(0.85, 0.92, 1.0, 1.0))
Colorbar(fig2[1,2], hm2, label="u_x / U∞")
# Hull silhouette on centerplane (rectangle below z=hull_zc)
poly!(ax2, [Point2f(hull_xc - L_c/2, hull_zc - T_c),
            Point2f(hull_xc + L_c/2, hull_zc - T_c),
            Point2f(hull_xc + L_c/2, hull_zc),
            Point2f(hull_xc - L_c/2, hull_zc)];
    color=(:grey20, 0.7), strokecolor=:black, strokewidth=1.2)
scatter!(ax2, [Point2f(prop_xc, prop_zc)];
    color=:orange, marker=:cross, markersize=16)
lines!(ax2, [rud_xc, rud_xc],
        [rud_zc - rudder.span/2, rud_zc + rudder.span/2];
    color=:lime, linewidth=3)
hlines!(ax2, [hull_zc]; color=:grey50, linestyle=:dot, linewidth=1)
save(out2, fig2)
@printf "Wrote %s\n" out2
