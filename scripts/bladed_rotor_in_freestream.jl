#!/usr/bin/env julia
#
# BladedRotor ↔ WaterLily coupling demo. Uniform inflow, one VLM
# blade-resolved propeller, smeared body force. Reports thrust and
# observed wake; compare visually to the SwirlingDisk wake.

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
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "80"))

const Vinf = 1f0
const Re   = 5000f0
const L_c  = 32f0
const ν_c  = Vinf * L_c / Re

const prop_x = Float32(NX / 4)
const prop_y = Float32(NY / 2)
const prop_z = Float32(NZ / 2)
const R_prop = 6.0
const rotor = BladedRotor(;
    N_blades = 3, R = R_prop, R_hub = 1.0,
    chord = (1.5, 0.9),
    twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 12, nc = 4,
)
const J = 0.7
const Ω = Float64(π) * Vinf / (J * R_prop)

@printf "=== BladedRotor + WaterLily coupling ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Rotor      = %d blades, R=%.1f, J=%.2f, Ω=%.3f\n" rotor.N_blades R_prop J Ω
@printf "  Position   = (%.1f, %.1f, %.1f)\n" prop_x prop_y prop_z
@printf "  Steps      = %d\n\n" NSTEPS
flush(stdout)

sim = WaterLily.Simulation((NX, NY, NZ),
    (Vinf, 0f0, 0f0), L_c;
    T = Float32, ν = ν_c, Δt = 0.5f0, ϵ = 1, U = Vinf,
)

# VLM-resolved thrust per step → smear into +x at the disk centre.
# (For now, ignore the swirl component; smearing the radial force
# distribution as a single body-force is the simplest one-way coupling.)
function rotor_udf(flow, t; kwargs...)
    r = rotor_forces(rotor, Vinf, Ω)
    # CT and CQ are normalised; reconstruct the actual force/torque.
    # Flip sign because the (twist > 0, Ω > 0) combination produces
    # CT < 0 in VortexLattice's convention — physically the propeller
    # is propelling forward, so we apply the magnitude in +x.
    Sref = π * R_prop^2
    thrust = abs(r.CT * 0.5 * Vinf^2 * Sref)
    force = SVector(Float32(thrust), 0f0, 0f0)
    smear_force!(flow.f, force, SVector(prop_x, prop_y, prop_z); ε=3.0f0)
    return nothing
end

@info "Running $NSTEPS steps with rotor forcing…"
t0 = time()
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf=rotor_udf,
        pois_tol=1f-6, pois_itmx=50)
    if step % 10 == 0
        u_max = maximum(abs, sim.flow.u)
        u_axial_at_prop = sim.flow.u[round(Int, prop_x + 1.5),
                                     round(Int, prop_y + 1.5),
                                     round(Int, prop_z + 1.5), 1]
        @info @sprintf("step=%3d  |u|=%.3f  u_x(prop)=%.3f  elapsed=%.1fs",
            step, u_max, u_axial_at_prop, time()-t0)
        flush(stdout)
    end
end

# Centerplane slice (y = prop_y): u_x heatmap.
j_slice = round(Int, prop_y + 1.5)
NX_α, NY_α, NZ_α = size(sim.flow.p)
ux_slice = Array{Float32}(undef, NX_α, NZ_α)
for i in 1:NX_α, k in 1:NZ_α
    ux_slice[i, k] = 0.5f0 * (sim.flow.u[i, j_slice, k, 1] +
                              sim.flow.u[min(i+1, NX_α), j_slice, k, 1])
end

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "bladed_rotor"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "wake_uniform_J07.png")

fig = Figure(size=(1100, 480))
ax = Axis(fig[1, 1]; aspect=DataAspect(),
    xlabel="x (cells, +x downstream)", ylabel="z (cells)",
    title=@sprintf("BladedRotor wake: u_x on y=%.1f, J=%.2f, %d steps", prop_y, J, NSTEPS))
ux_range = (0.8f0, 1.5f0)
hm = heatmap!(ax, 1:NX_α, 1:NZ_α, ux_slice;
    colormap=:roma, colorrange=ux_range)
Colorbar(fig[1, 2], hm, label="u_x / U∞")
θ = range(0f0, 2f0π, length=64)
disk_z = prop_z .+ R_prop .* sin.(θ)
# Project disk onto x-z plane (it's a vertical circle in y-z, so on
# this x-z slice it shows as a line segment).
lines!(ax, [prop_x, prop_x], [prop_z - R_prop, prop_z + R_prop];
    color=:orange, linewidth=3, label="rotor disk")

save(out, fig)
@printf "\nWrote %s\n" out
