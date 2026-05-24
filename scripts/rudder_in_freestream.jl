#!/usr/bin/env julia
#
# End-to-end Rudder ↔ WaterLily coupling demo using LiftingSurfaces.jl.
# Setup: uniform inflow, one rudder at δ angle, smeared body force.
# Output: a still image showing the wake deflection behind the rudder.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["CairoMakie", "VortexLattice"]; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, LiftingSurfaces
using LiftingSurfaces: SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "120"))
const δ_DEG  = parse(Float32, get(ENV, "WL_DELTA", "10"))

const Vinf = 1f0
const Re   = 5000f0
const ν_c  = Vinf * Float32(NX/3) / Re

const rudder_x = Float32(NX / 3)
const rudder_y = Float32(NY / 2)
const rudder_z = Float32(NZ / 2)
const rudder = Rudder(; chord=4.0, span=8.0, ns=12, nc=6)

@printf "=== Rudder + WaterLily coupling demo ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Rudder     = chord=%.1f, span=%.1f at (%.1f, %.1f, %.1f), δ=%.1f°\n" rudder.chord rudder.span rudder_x rudder_y rudder_z δ_DEG
@printf "  Steps      = %d\n" NSTEPS
flush(stdout)

sim = WaterLily.Simulation((NX, NY, NZ),
    (Vinf, 0f0, 0f0), Float32(NX/3);
    T = Float32, ν = ν_c, Δt = 0.5f0, ϵ = 1, U = Vinf,
)

# Convert rudder polar to force (Newton-scale, cell-units) then smear.
# Side force ≈ CY * 0.5·ρ·V² · Sref (ρ=1 cell). Sref = chord × span.
function rudder_udf(flow, t; kwargs...)
    Sref = rudder.chord * rudder.span
    # No inflow-sampling on this first cut — pass nothing; the rudder
    # sees pure freestream Vinf. (Two-way coupling is the next step.)
    r = rudder_forces(rudder, deg2rad(δ_DEG), Vinf)
    side  = r.CL * 0.5f0 * Vinf^2 * Sref     # CL in wind frame → z-side
    drag  = r.CD * 0.5f0 * Vinf^2 * Sref
    # CL points in wind-frame z, which is the rudder's chord-perpendicular
    # direction; with rudder spanning +y and chord along +x, CL acts in z.
    force = SVector(-Float32(drag), 0f0, Float32(side))
    smear_force!(flow.f, force, SVector(rudder_x, rudder_y, rudder_z); ε=2.0f0)
    return nothing
end

@info "Running $NSTEPS steps with rudder forcing…"
t0 = time()
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf=rudder_udf,
        pois_tol=1f-6, pois_itmx=50)
    if step % 20 == 0
        u_max = maximum(abs, sim.flow.u)
        @info @sprintf("step=%3d  |u|=%.3f  elapsed=%.1fs", step, u_max, time()-t0)
        flush(stdout)
    end
end

# Side-view slice through rudder y-z midplane → vertical x-z plane.
# Actually the wake bends in z (because force is in +z) — easier to see
# on the y=rudder_y plane: x-z slice showing u_z.
j_slice = round(Int, rudder_y + 1.5)
NX_α, NY_α, NZ_α = size(sim.flow.p)
uz_slice = Array{Float32}(undef, NX_α, NZ_α)
for i in 1:NX_α, k in 1:NZ_α
    uz_slice[i, k] = 0.5f0 * (sim.flow.u[i, j_slice, k, 3] +
                              sim.flow.u[i, j_slice, min(k+1, NZ_α), 3])
end

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, @sprintf("wake_delta%02d.png", round(Int, δ_DEG)))

fig = Figure(size=(1000, 500))
ax = Axis(fig[1, 1]; aspect=DataAspect(),
    xlabel="x (cells, +x downstream)", ylabel="z (cells)",
    title=@sprintf("Rudder wake: u_z on y=%.1f, δ=%.1f° after %d steps", rudder_y, δ_DEG, NSTEPS))
uz_max = quantile(abs.(vec(uz_slice)), 0.98)
hm = heatmap!(ax, 1:NX_α, 1:NZ_α, uz_slice;
    colormap=:RdBu, colorrange=(-uz_max, uz_max))
Colorbar(fig[1, 2], hm, label="u_z")
scatter!(ax, [Point2f(rudder_x, rudder_z)];
    color=:lime, marker=:dtriangle, markersize=14)
lines!(ax, [rudder_x, rudder_x], [rudder_z - 4, rudder_z + 4];
    color=:lime, linewidth=2)

save(out, fig)
@printf "\nWrote %s\n" out
