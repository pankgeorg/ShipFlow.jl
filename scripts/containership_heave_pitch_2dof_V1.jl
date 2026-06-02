#!/usr/bin/env julia
#
# K1: heave + pitch 2-DOF. Extends J1 with a pitch DOF about the
# y-axis through the centre of mass. Pitch moment from
# `WaterLily.pressure_moment(x_CG, sim)[2]` (y component, negated
# to match the body→fluid sign convention).
#
# Equations of motion:
#   M · z̈ = −F_p,z − M·g                 (heave)
#   I · θ̈ = −M_p,y                       (pitch, no gravity moment
#                                          for symmetric body around CG)
#
# `I = M · L²/12` is the slender-body pitch inertia approximation.

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

const NX = 128; const NY = 64; const NZ = 32
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

const V0      = 0.75 * L_c * B_c * T_c
const M_ship  = ρ_w * V0
const I_pitch = M_ship * L_c^2 / 12.0

z_h    = Ref(0.0); zdot_h = Ref(0.0)
θ      = Ref(0.0); θdot   = Ref(0.0)

# Map: subtract CG translation (z_h), then unrotate by θ around y.
hull_map = (x, t) -> begin
    dx = x[1] - hull_xc
    dz = x[3] - (hull_zc0 + z_h[])
    cθ, sθ = cos(θ[]), sin(θ[])
    SVector(cθ * dx + sθ * dz, x[2] - hull_yc, -sθ * dx + cθ * dz)
end
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
    deck_h = T_c / 2,        # L1: deck for metacentric stability
    map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== W2: heave + pitch 2-DOF (Wigley) ===\n"
@printf "  Grid     = %d × %d × %d   NSTEPS=%d\n" NX NY NZ NSTEPS
@printf "  M=%.0f  I_pitch=%.0f  g_cell=%.4f\n" M_ship I_pitch G_c
@printf "  Estimated periods: T_heave≈17 ct, T_pitch≈14 ct\n"
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
    T=Float64, ν = VoF.viscosity(vof),
    g=(i, x, t) -> i == 3 ? -G_c : 0.0,
    Δt=0.25, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=U∞,
)

const GRAVITY_RAMP_STEPS = NSTEPS ÷ 4

z_hist = Float64[]; zdot_hist = Float64[]
θ_hist = Float64[]; θdot_hist = Float64[]
Fhz_hist = Float64[]; My_hist = Float64[]

x_CG = [hull_xc, hull_yc, hull_zc0]

t0 = time()
for step in 1:NSTEPS
    WaterLily.measure!(sim, sum(sim.flow.Δt))
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=1e-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    g_scale = min(1.0, step / GRAVITY_RAMP_STEPS)
    g_now = G_c * g_scale
    Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
    Mp = WaterLily.pressure_moment(x_CG, sim)
    F_buoy = -(Fp[3] + Fv[3])
    F_grav = -M_ship * g_now
    F_net  = F_buoy + F_grav
    # Pitch moment about y (component 2). Negate for fluid→body.
    M_y_bdim = -Mp[2]
    # V1: explicit analytical metacentric restoring on top of the
    # BDIM-measured moment. BDIM under-captures the waterplane-
    # intersection moment at our grid resolution, leaving the wave-
    # forcing component without an in-kind restoring response.
    # For the Wigley waterplane (rectangular waterplane).  B·L³/12.
    K_pitch = ρ_w * g_now * B_c * L_c^3 / 12
    M_y_restore = -K_pitch * θ[]
    M_y = M_y_bdim + M_y_restore

    # Artificial damping: needs to be strong because Wigley's
    # parabolic waterplane gives weaker K_pitch (B·L³/20) than the
    # Containership (B·L³/12). Per U1, β=0.2/dt was enough for
    # Containership but Wigley needs β=0.5/dt to converge.
    β_pitch = 0.5 / sim.flow.Δt[end-1]
    β_heave = 0.1 / sim.flow.Δt[end-1]
    z_ddot = F_net / M_ship - β_heave * zdot_h[]
    θ_ddot = M_y / I_pitch - β_pitch * θdot[]
    zdot_h[] += z_ddot * sim.flow.Δt[end-1]
    z_h[]    += zdot_h[] * sim.flow.Δt[end-1]
    θdot[]   += θ_ddot * sim.flow.Δt[end-1]
    θ[]      += θdot[] * sim.flow.Δt[end-1]
    # Hard clamp to keep the body inside the domain even if damping
    # fails — the hull should oscillate at most ±5° in steady state.
    θ[] = clamp(θ[], -deg2rad(20), deg2rad(20))
    z_h[] = clamp(z_h[], -T_c, T_c)

    push!(z_hist, z_h[]); push!(zdot_hist, zdot_h[])
    push!(θ_hist, θ[]);   push!(θdot_hist, θdot[])
    push!(Fhz_hist, F_buoy); push!(My_hist, M_y)

    if step % 40 == 0
        @printf "  step=%3d  z=%+.3f ż=%+.3f  θ=%+.4f rad (%+.2f°)  θ̇=%+g  F_buoy=%+.1f M_y=%+.1f  |u|=%.3f  elapsed=%.1fs\n" step z_h[] zdot_h[] θ[] (θ[]*180/π) θdot[] F_buoy M_y maximum(abs, sim.flow.u) (time()-t0)
        flush(stdout)
    end
end

tail = (NSTEPS - NSTEPS÷4 + 1):NSTEPS
z_eq = mean(z_hist[tail])
θ_eq = mean(θ_hist[tail])
@printf "\nFinal: ⟨z_h⟩_tail = %+.4f cells   ⟨θ⟩_tail = %+.4f rad (%+.2f°)\n" z_eq θ_eq (θ_eq*180/π)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "heave_pitch_2dof_containership_V1"))
mkpath(OUTDIR)
xs = collect(1:NSTEPS)
fig = Figure(size=(950, 900))
ax1 = Axis(fig[1, 1]; ylabel="z_h (cells)", title="Heave + pitch 2-DOF")
lines!(ax1, xs, z_hist; color=:steelblue, linewidth=2)
hlines!(ax1, [0]; color=:grey, linestyle=:dash)
vlines!(ax1, [GRAVITY_RAMP_STEPS]; color=:tomato, linestyle=:dash)
ax2 = Axis(fig[2, 1]; ylabel="θ (degrees)")
lines!(ax2, xs, θ_hist .* (180/π); color=:purple, linewidth=2)
hlines!(ax2, [0]; color=:grey, linestyle=:dash)
ax3 = Axis(fig[3, 1]; ylabel="F_buoy (cell units)")
lines!(ax3, xs, Fhz_hist; color=:tomato, linewidth=2)
hlines!(ax3, [M_ship * G_c]; color=:steelblue, linestyle=:dash, label="M·g")
axislegend(ax3)
ax4 = Axis(fig[4, 1]; xlabel="step", ylabel="M_y (cell units)")
lines!(ax4, xs, My_hist; color=:darkgreen, linewidth=2)
hlines!(ax4, [0]; color=:grey, linestyle=:dash)
save(joinpath(OUTDIR, "heave_pitch.png"), fig)
println("Wrote $(joinpath(OUTDIR, "heave_pitch.png"))")
