#!/usr/bin/env julia
#
# Headline F3 demo, refactored to use ShipMakie recipes instead of
# bespoke heatmap / scatter / hand-rolled `eta_surface!` code.
#
# Side-by-side comparison with `headline_demo.jl`:
#   - This script:               ~120 lines, mostly simulation setup
#   - Original `headline_demo.jl`: 225 lines, with custom viz helpers
#
# Both produce equivalent output. This one is meant as a real-world
# showcase of ShipMakie's value.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipMakie.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "CairoMakie"]; io=devnull)

using WaterLily, VoF, Turbulence, ShipShapes, LiftingSurfaces, ShipMakie
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie
CairoMakie.activate!()

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NFRAMES = parse(Int, get(ENV, "WL_NFRAMES", "60"))
const BURNIN  = parse(Int, get(ENV, "WL_BURNIN", "15"))
const FR = parse(Float32, get(ENV, "WL_FR", "0.30"))
const δ_DEG = parse(Float32, get(ENV, "WL_DELTA", "10"))
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
const hull_xc = Float32(NX/5); const hull_yc = Float32(NY/2); const hull_zc = H_w_c

const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J = 0.32
const Ω = Float64(π) * U∞ / (J * R_prop)
const rud_xc = Float32(prop_xc + R_prop * 0.5)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

# --- Sim setup -------------------------------------------------------
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
    par_frac = PAR_FRAC, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
vof  = VoFFlow((NX, NY, NZ); α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
turb = WALE((NX, NY, NZ); Cw=0.5f0, ν₀=0f0)
rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)
r_rot = rotor_forces(rotor, U∞, Ω)
thrust = abs(r_rot.CT * 0.5 * U∞^2 * π * R_prop^2)
torque = r_rot.CQ * 0.5 * U∞^2 * π * R_prop^2 * R_prop

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
end
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T=Float32, ν = Turbulence.viscosity(turb),
    g=(i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=U∞,
)

function vlm_udf(flow, t; kwargs...)
    smear_blades!(flow.f, thrust, torque,
        SVector(prop_xc, prop_yc, prop_zc), SVector(1.0, 0.0, 0.0),
        Float64(R_prop), 0.2 * Float64(R_prop);
        N_blades=3, N_sections=4, ε=1.5)
    inflow_fn = trilinear_inflow(flow.u)
    r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞; inflow=inflow_fn)
    side = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    drag = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
        SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
    return nothing
end

@info "Burn-in $BURNIN steps…"
for _ in 1:BURNIN
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    update_νt!(turb, sim.flow.u, vof.ν)
end

# --- Render with ShipMakie recipes -----------------------------------
mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c/2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c/2

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "headline_demo_shipmakie", "frames"))
mkpath(OUTDIR)
@info "Rendering $NFRAMES frames using ShipMakie recipes…"
t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf=vlm_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
    update_νt!(turb, sim.flow.u, vof.ν)

    fig = Figure(size=(1200, 720))
    ax_top = Axis(fig[1, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="y (cells)",
        title=@sprintf("Headline (ShipMakie)  frame %d/%d", frame, NFRAMES))
    etaheatmap!(ax_top, vof.α;
        waterline_z = H_w_c, mask = mask_inside,
        hull_box = (hull_xc, hull_yc, L_c, B_c),
        colorrange = (-0.4, 0.4))
    Makie.scatter!(ax_top, [Point2f(prop_xc, prop_yc)];
        color=:orange, marker=:cross, markersize=16, strokewidth=2)
    Makie.scatter!(ax_top, [Point2f(rud_xc, rud_yc)];
        color=:lime, marker=:dtriangle, markersize=14, strokewidth=2)

    ax_side = Axis(fig[2, 1]; aspect=DataAspect(),
        xlabel="x (cells)", ylabel="z (cells)",
        title="centreline u_x  (side view)")
    velocityslice!(ax_side, sim.flow.u;
        slice_axis = :y, index = NY ÷ 2, component = 1,
        colorrange = (-1.5, 3.0))
    Makie.hlines!(ax_side, [H_w_c]; color=:cyan, linestyle=:dash)
    Makie.scatter!(ax_side, [Point2f(prop_xc, prop_zc)];
        color=:orange, marker=:cross, markersize=16, strokewidth=2)
    Makie.scatter!(ax_side, [Point2f(rud_xc, rud_zc)];
        color=:lime, marker=:dtriangle, markersize=14, strokewidth=2)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 10 == 0
        @printf "frame %3d / %d  elapsed=%.1fs\n" frame NFRAMES (time() - t0)
        flush(stdout)
    end
end
@printf "Done: %d frames in %s\n" NFRAMES OUTDIR
