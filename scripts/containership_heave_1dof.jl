#!/usr/bin/env julia
#
# J1: heave 1-DOF. The hull is no longer rigidly fixed in z — it
# floats according to the balance of:
#   * gravity                 (constant, M·g)
#   * buoyancy                (Archimedes on the submerged Wigley volume)
#   * vertical hydrodynamic   (BDIM-kernel measured each step)
#
# Equations of motion (semi-implicit Euler):
#   z̈_h = (F_b(z_h) − M·g + F_hydro_z) / M
#   ż_h += z̈_h · dt
#   z_h += ż_h · dt
#
# The hull's world-frame z-offset is `z_h`; the map function reads it
# from a Ref so the WaterLily Simulation sees the moving body.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)

using WaterLily, VoF, ShipShapes
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics, CairoMakie

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "240"))

const L_c = 36.0; const B_c = 8.0; const T_c = 5.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞  = 1.0
const Fr  = 0.25
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = 5000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/5; const hull_yc = NY/2; const hull_zc0 = H_w_c

# Containership submerged volume at z_h = 0 (par_frac=0.5 ⇒ Cb=0.75).
# We do NOT use the analytical formula in the EOM — instead we use
# the BDIM-measured F_hyz directly. But we need V₀ to set M_ship.
const V0 = 0.75 * L_c * B_c * T_c
const M_analytical = ρ_w * V0
# Per R1 (RESULTS-heave-containership.md), the BDIM Archimedes
# under-measures the analytical buoyancy for the Containership by
# ~63 % (sharp corners). To make the hull float at z_h=0 we set
# `M_ship = bias · M_analytical` where bias is measured by a
# hydrostatic pre-pass below.
const BIAS_FACTOR = parse(Float64, get(ENV, "WL_BIAS", "1.0"))
const M_ship = BIAS_FACTOR * M_analytical

# Mutable heave state.
z_h    = Ref(0.0)
zdot_h = Ref(0.0)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - (hull_zc0 + z_h[]))
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
    deck_h = T_c / 2,        # T1: deck for proper Archimedes at z_h ≠ 0
    map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== R1: heave 1-DOF (Containership) ===\n"
@printf "  Grid     = %d × %d × %d   NSTEPS=%d\n" NX NY NZ NSTEPS
@printf "  L=%.1f B=%.1f T=%.1f   V₀=%.3f   M_ship=%.3f\n" L_c B_c T_c V0 M_ship
@printf "  Fr=%.2f  Re=%.0f  g_cell=%.4f\n" Fr Re G_c
@printf "  Spring K = ρ_w·g·2BL/3 = %.2f  →  T_h = %.2f cell-time = %.1f steps\n" (ρ_w*G_c*2*B_c*L_c/3) (2π*sqrt(M_ship/(ρ_w*G_c*2*B_c*L_c/3))) (2π*sqrt(M_ship/(ρ_w*G_c*2*B_c*L_c/3))/0.25)
flush(stdout)

vof = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float64)
function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
end
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0.0, 0.0), L_c;
    T=Float64, ν=vof.ν,
    g=(i, x, t) -> i == 3 ? -G_c : 0.0,
    Δt=0.25, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=U∞,
)

# Ramp gravity from 0 → G_c over the first quarter of the run so we don't
# get an initial slamming impulse from cold-starting the buoyancy/gravity
# balance against zero velocity field.
const GRAVITY_RAMP_STEPS = NSTEPS ÷ 4

z_hist     = Float64[]
zdot_hist  = Float64[]
Fbnet_hist = Float64[]
Fhz_hist   = Float64[]

t0 = time()
for step in 1:NSTEPS
    # Re-measure the BDIM kernel each step so the moving map function
    # (which reads z_h[]) propagates into the flow's μ₀, μ₁ coefficients.
    WaterLily.measure!(sim, sum(sim.flow.Δt))
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1e-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    # WaterLily.pressure_force / viscous_force return the force the
    # BODY exerts ON the FLUID (drag scripts use −Fp[1] for drag in
    # the inflow direction). So the force ON the body in z is −Fp[3].
    # At hydrostatic equilibrium this gives the Archimedes upthrust.
    g_scale = min(1.0, step / GRAVITY_RAMP_STEPS)
    g_now = G_c * g_scale
    Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
    F_buoy = -(Fp[3] + Fv[3])         # +z = up; positive at equilibrium
    F_grav = -M_ship * g_now           # negative
    F_net  = F_buoy + F_grav
    F_hyz  = F_buoy                    # keep variable for the print below
    # Strong damping: the Containership SDF has sharp transitions
    # that make explicit Euler unstable at the natural heave
    # frequency. β=0.5/dt makes the system overdamped but stable.
    β_heave = 0.5 / sim.flow.Δt[end-1]
    z_ddot = F_net / M_ship - β_heave * zdot_h[]
    zdot_h[] += z_ddot * sim.flow.Δt[end-1]
    z_h[]    += zdot_h[] * sim.flow.Δt[end-1]
    # Hard clamp so even if damping fails the body stays in domain.
    z_h[] = clamp(z_h[], -T_c, T_c)

    push!(z_hist, z_h[])
    push!(zdot_hist, zdot_h[])
    push!(Fbnet_hist, F_grav)
    push!(Fhz_hist, F_hyz)
    if step ≤ 5 || step % 30 == 0
        @printf "  step=%3d  z_h=%+.4f  ż=%+.4f  F_hyz=%+g  F_grav=%+g  F_net=%+g  Fp_z=%+g  Fv_z=%+g  |u|=%.3f\n" step z_h[] zdot_h[] F_hyz F_grav F_net Fp[3] Fv[3] maximum(abs, sim.flow.u)
        flush(stdout)
    end
end

# Equilibrium estimate: average over last 25 %
tail = (NSTEPS - NSTEPS÷4 + 1):NSTEPS
z_eq = mean(z_hist[tail])
zdot_rms = sqrt(mean(zdot_hist[tail].^2))
@printf "\nFinal: ⟨z_h⟩_tail = %+.4f  zdot_rms_tail = %.4f\n" z_eq zdot_rms

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "heave_1dof_containership"))
mkpath(OUTDIR)
xs = collect(1:NSTEPS)
fig = Figure(size=(950, 720))
ax1 = Axis(fig[1, 1]; ylabel="z_h (cells)",
    title="Heave 1-DOF response (gravity ramped over first $(GRAVITY_RAMP_STEPS) steps)")
lines!(ax1, xs, z_hist; color=:steelblue, linewidth=2)
hlines!(ax1, [0]; color=:grey, linestyle=:dash)
vlines!(ax1, [GRAVITY_RAMP_STEPS]; color=:tomato, linestyle=:dash)

ax2 = Axis(fig[2, 1]; ylabel="ż_h (cells/step)")
lines!(ax2, xs, zdot_hist; color=:tomato, linewidth=2)
hlines!(ax2, [0]; color=:grey, linestyle=:dash)

ax3 = Axis(fig[3, 1]; xlabel="step", ylabel="F (cell units)")
lines!(ax3, xs, Fbnet_hist; color=:steelblue, linewidth=2, label="−M·g")
lines!(ax3, xs, Fhz_hist;   color=:tomato, linewidth=2, label="F_hyz (≈ buoyancy + dynamic)")
hlines!(ax3, [0]; color=:grey, linestyle=:dash)
axislegend(ax3; position=:rt)

save(joinpath(OUTDIR, "heave.png"), fig)
println("Wrote $(joinpath(OUTDIR, "heave.png"))")
