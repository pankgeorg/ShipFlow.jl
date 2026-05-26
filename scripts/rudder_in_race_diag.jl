#!/usr/bin/env julia
#
# I3: diagnostic for the F2 result. The F2 trial showed CL_rudder
# +300 % when the rudder sits inside the rotor race, but the hull
# integrated F_hull,y went slightly DOWN (-2.49 → -2.19). This
# script visualises the per-x-station side-force distribution on the
# hull for the two cases, to localise where the deficit appears.
#
# Idea: integrate pressure × normal · ŷ over hull cells at each x
# slice. The hull is masked via the AutoBody μ₀ kernel (zero inside,
# 1 outside) — local pressure × ∂μ₀/∂y gives the y-direction body force
# density. Summing over y, z at fixed x gives F_y(x).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "CairoMakie"]; io=devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics, CairoMakie

const NX = 128; const NY = 64; const NZ = 32
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "80"))
const δ_DEG  = 10f0
const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5); const hull_yc = Float32(NY/2); const hull_zc = H_w_c

const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J_op = 0.32
const Ω = Float64(π) * U∞ / (J_op * R_prop)
const rud_xc = Float32(prop_xc + R_prop * 0.5)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)
r_rot = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop

function vof_pois_ctor(vof)
    (flow) -> begin
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
end

# Per-x-station hull side force F_y(x) from the BDIM kernel + pressure.
# `μ₀[I, d]` is the BDIM kernel in direction d at face d of cell I.
# The body force on the fluid is -p · ∂μ₀/∂y so the force ON the body is
# +p · ∂μ₀/∂y · cell_area. We approximate ∂μ₀/∂y via central difference
# on the cell-centred μ₀ (which sees the SDF).
function per_x_force_y(sim)
    p = sim.flow.p
    μ₀ = sim.flow.μ₀
    nx, ny, nz, _ = size(μ₀)
    Fy_of_x = zeros(Float64, nx)
    for i in 2:nx-1, j in 2:ny-1, k in 2:nz-1
        # ∂μ₀_y / ∂y at cell centre via central diff
        dμy = (μ₀[i, j+1, k, 2] - μ₀[i, j-1, k, 2]) * 0.5
        # Force ON the body = pressure × ∂μ₀/∂y (sign convention check)
        Fy_of_x[i] += p[i, j, k] * dμy
    end
    return Fy_of_x
end

function run_case(label; twoway::Bool)
    @info "Running $label"
    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=vof.ν,
        g=(i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor(vof), U=U∞,
    )
    function combo_udf(flow, t; kwargs...)
        smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                     SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
        smear_torque!(flow.f, Float32(torque),
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                      N=8, ε=2.0f0)
        inflow_fn = twoway ? trilinear_inflow(flow.u) : nothing
        r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞; inflow=inflow_fn)
        side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                     SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
        return nothing
    end
    Fy_history = []   # collect last-25% per-x forces
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; udf=combo_udf, pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        if step > NSTEPS - NSTEPS÷4
            push!(Fy_history, per_x_force_y(sim))
        end
        if step % 20 == 0
            @info @sprintf("  [%s] step=%3d  |u|=%.3f", label, step, maximum(abs, sim.flow.u))
        end
    end
    # Average over the last 25% of steps
    Fy_mean = mean(Fy_history)
    return Fy_mean
end

Fy_off = run_case("TWOWAY=0"; twoway=false)
Fy_on  = run_case("TWOWAY=1"; twoway=true)

# Plot
xs = collect(1:length(Fy_off))
fig = Figure(size=(950, 480))
ax = Axis(fig[1, 1]; xlabel="x (cells)", ylabel="dF_y / dx",
    title="Hull side force distribution along x: TWOWAY=0 vs TWOWAY=1")
lines!(ax, xs, Fy_off; color=:steelblue, linewidth=2, label="TWOWAY=0")
lines!(ax, xs, Fy_on;  color=:tomato, linewidth=2, label="TWOWAY=1")
hlines!(ax, [0]; color=:grey, linestyle=:dash)
# Mark hull extent
vlines!(ax, [hull_xc - L_c/2, hull_xc + L_c/2]; color=:black, linestyle=:dot)
# Mark rotor + rudder
vlines!(ax, [prop_xc]; color=:orange, linestyle=:dash, label="rotor")
vlines!(ax, [rud_xc]; color=:lime, linestyle=:dash, label="rudder")
axislegend(ax)
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder_in_race"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "Fy_per_x.png")
save(out, fig)

# Numerics
F_off_total = sum(Fy_off)
F_on_total = sum(Fy_on)
F_off_hull_only = sum(Fy_off[round(Int, hull_xc - L_c/2):round(Int, hull_xc + L_c/2)])
F_on_hull_only  = sum(Fy_on[round(Int, hull_xc - L_c/2):round(Int, hull_xc + L_c/2)])
@printf "\nTotal F_y on body:  TWOWAY=0 %+.3f   TWOWAY=1 %+.3f\n" F_off_total F_on_total
@printf "F_y over hull span: TWOWAY=0 %+.3f   TWOWAY=1 %+.3f\n" F_off_hull_only F_on_hull_only
@printf "Wrote %s\n" out
