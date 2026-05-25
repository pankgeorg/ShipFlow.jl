#!/usr/bin/env julia
#
# Rudder-effectiveness trial: drive δ(t) in a tri-state ramp profile
# (0 → 15° hold → 0) and record CY_rudder vs F_hull_y over time. This
# is the canonical manoeuvring experiment — the smallest demonstration
# that the lifting-surface tier produces a meaningful side-force
# history a manoeuvring-model coupler would consume.
#
# G2 of NEW_PLAN.md.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add(["VortexLattice", "CairoMakie"]; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, CairoMakie, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "180"))
const δ_MAX_DEG = parse(Float32, get(ENV, "WL_DELTA_MAX", "15"))
const RAMP_UP = NSTEPS ÷ 3
const HOLD    = NSTEPS ÷ 3
# rest is ramp-down

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

# Use the G1 self-prop operating point.
const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J = 0.32          # ≈ G1's self-prop J
const Ω = Float64(π) * U∞ / (J * R_prop)

const rud_xc = Float32(prop_xc + 2 * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

@printf "=== Rudder-effectiveness trial ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Steps      = %d  (ramp_up=%d, hold=%d, ramp_down=%d)\n" NSTEPS RAMP_UP HOLD (NSTEPS - RAMP_UP - HOLD)
@printf "  δ_max      = %.1f°\n" δ_MAX_DEG
@printf "  Rotor J=%.2f (G1 self-prop point)\n\n" J
flush(stdout)

vof = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)

# Compute rotor CT, CQ once.
r_rot  = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop
@printf "  rotor |CT|=%.3f, thrust=%.3f, |CQ|=%.3f, torque=%.3f\n" abs(r_rot.CT) thrust abs(r_rot.CQ) torque

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
end

sim = WaterLily.Simulation((NX, NY, NZ),
    (U∞, 0f0, 0f0), L_c;
    T=Float32, ν=vof.ν,
    g=(i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=U∞,
)

# Tri-state ramp δ(step)
function δ_at(step)
    if step ≤ RAMP_UP
        return Float32(δ_MAX_DEG) * step / RAMP_UP
    elseif step ≤ RAMP_UP + HOLD
        return Float32(δ_MAX_DEG)
    else
        s = step - RAMP_UP - HOLD
        return Float32(δ_MAX_DEG) * (1 - s / (NSTEPS - RAMP_UP - HOLD))
    end
end

# Histories
step_hist   = Int[]
δ_hist      = Float32[]
CL_rud_hist = Float32[]
Fhy_hist    = Float32[]   # hull side-force in body-y

@info "Running trial…"
t0 = time()
for step in 1:NSTEPS
    δ_deg = δ_at(step)
    # Rotor: thrust + torque
    function combo_udf(flow, t; kwargs...)
        smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                     SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
        smear_torque!(flow.f, Float32(torque),
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                      N=8, ε=2.0f0)
        # Rudder
        r_rud = rudder_forces(rudder, deg2rad(δ_deg), U∞)
        side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
        smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                     SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
        # Stash rudder CL for the time-series (just the freestream
        # response; doesn't include flow.u sampling).
        push!(CL_rud_hist, Float32(r_rud.CL))
        return nothing
    end

    WaterLily.mom_step!(sim.flow, sim.pois; udf=combo_udf, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    # Pop the most recent CL — combo_udf pushes once per step.
    # Hull y-side-force = pressure_force[2] + viscous_force[2]
    Fp = WaterLily.pressure_force(sim)
    Fv = WaterLily.viscous_force(sim)
    F_hy = Float32(Fp[2] + Fv[2])

    push!(step_hist, step)
    push!(δ_hist, δ_deg)
    push!(Fhy_hist, F_hy)

    if step % 20 == 0
        u_max = maximum(abs, sim.flow.u)
        @info @sprintf("step=%3d  δ=%5.2f°  |u|=%.3f  F_hy=%+.3f  elapsed=%.1fs",
            step, δ_deg, u_max, F_hy, time()-t0)
        flush(stdout)
    end
end
# combo_udf is called once per step but the runtime may call it twice
# (predictor + corrector). Sub-sample the CL stream to one per step.
n = length(CL_rud_hist) ÷ NSTEPS
CL_rud_hist = CL_rud_hist[n:n:end]

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder_trial"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "trial.png")

fig = Figure(size=(1100, 700))
ax_δ = Axis(fig[1,1]; ylabel="δ (°)", title="Rudder-effectiveness trial")
lines!(ax_δ, step_hist, δ_hist; color=:black, linewidth=2)
hidexdecorations!(ax_δ, grid=false)
ax_F = Axis(fig[2,1]; ylabel="F_hull,y  (cell-units)")
lines!(ax_F, step_hist, Fhy_hist; color=:steelblue, linewidth=2)
hlines!(ax_F, [0]; color=:grey, linestyle=:dash)
hidexdecorations!(ax_F, grid=false)
ax_C = Axis(fig[3,1]; xlabel="step", ylabel="CL_rudder")
lines!(ax_C, step_hist, CL_rud_hist; color=:tomato, linewidth=2)
hlines!(ax_C, [0]; color=:grey, linestyle=:dash)
save(out, fig)
@printf "\nWrote %s\n" out

# Quick numerical summary
F_hy_max  = maximum(Fhy_hist[RAMP_UP:RAMP_UP+HOLD])
F_hy_min  = minimum(Fhy_hist[RAMP_UP:RAMP_UP+HOLD])
CL_hold   = mean(CL_rud_hist[RAMP_UP:RAMP_UP+HOLD])
@printf "Summary over hold phase (steps %d–%d):\n" RAMP_UP RAMP_UP+HOLD
@printf "  Rudder CL          = %.3f\n" CL_hold
@printf "  Hull F_y range     = [%+.3f, %+.3f]\n" F_hy_min F_hy_max
@printf "  Ratio |F_hull|/CL ≈ %.2f (a maneuvering-coupler input)\n" max(abs(F_hy_max), abs(F_hy_min)) / abs(CL_hold)
