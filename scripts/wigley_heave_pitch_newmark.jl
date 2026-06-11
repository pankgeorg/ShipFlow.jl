#!/usr/bin/env julia
#
# Newmark-β 2-DOF heave+pitch for the Wigley hull — the "proper
# integrator + added-mass model" follow-up to K1/V1
# (wigley_heave_pitch_2dof.jl). Setup (grid, hull, VoF, force
# measurement) matches that script; only the rigid-body update differs.
#
# Integrator: Newmark average-acceleration (β=1/4, γ=1/2,
# unconditionally stable for the linear part), per-DOF scalar solve:
#
#   (M_eff + γ·dt·C + β·dt²·K) a⁺ =
#         F_meas − C·(v + (1−γ)·dt·a) − K·(dt·v + (½−β)·dt²·a)
#
# with
#   * M_eff = M + A — added mass A = Ca·ρ∇ (heave) / Ca·ρ∇·L²/12 (pitch),
#     Ca env-tunable (defaults 1.0 — slender-body O(1) estimates);
#   * K — analytic hydrostatic stiffness acting on the *increment* only
#     (K·Δz, not K·z), so equilibrium is still set by the measured BDIM
#     force (no double-counting of the partially-captured hydrostatics)
#     while the staggered coupling gets an implicit restoring term at
#     the frequencies it under-resolves. Wigley waterplane:
#     A_wp = (2/3)·L·B  ⇒  K_z = ρg·A_wp ;  I_wp = B·L³/20 ⇒ K_θ = ρg·I_wp.
#   * C — radiation-damping stand-in as a fraction ζ of critical,
#     C = 2ζ√(K·M_eff)  (default ζ=0.2 — c.f. the explicit script's
#     effectively-overdamped β=0.5/dt ad-hoc term).
#
# ENV: WL_NSTEPS (240), WL_ZETA (0.2), WL_CA33 (1.0), WL_CA55 (1.0)

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

const V0      = 4 * L_c * B_c * T_c / 9
const M_ship  = ρ_w * V0
const I_pitch = M_ship * L_c^2 / 12.0

# Newmark / structural parameters
const βN = 0.25; const γN = 0.5
const ζ    = parse(Float64, get(ENV, "WL_ZETA", "0.2"))
const Ca33 = parse(Float64, get(ENV, "WL_CA33", "1.0"))
const Ca55 = parse(Float64, get(ENV, "WL_CA55", "1.0"))
const A33  = Ca33 * ρ_w * V0
const A55  = Ca55 * ρ_w * V0 * L_c^2 / 12
const Mz_eff = M_ship + A33
const Iθ_eff = I_pitch + A55
# hydrostatic stiffness (analytic, Wigley): full-load values at z=0
const A_wp = 2/3 * L_c * B_c          # ∫ B(1-(2x/L)²) dx
const I_wp = B_c * L_c^3 / 20         # ∫ B(1-(2x/L)²) x² dx

z_h    = Ref(0.0); zdot_h = Ref(0.0); zacc = Ref(0.0)
θ      = Ref(0.0); θdot   = Ref(0.0); θacc = Ref(0.0)

hull_map = (x, t) -> begin
    dx = x[1] - hull_xc
    dz = x[3] - (hull_zc0 + z_h[])
    cθ, sθ = cos(θ[]), sin(θ[])
    SVector(cθ * dx + sθ * dz, x[2] - hull_yc, -sθ * dx + cθ * dz)
end
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c,
    deck_h = T_c / 2, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== Newmark-β heave+pitch 2-DOF (Wigley) ===\n"
@printf "  Grid %d×%d×%d  NSTEPS=%d   β=%.2f γ=%.2f ζ=%.2f  Ca33=%.1f Ca55=%.1f\n" NX NY NZ NSTEPS βN γN ζ Ca33 Ca55
@printf "  M=%.0f (+A33 %.0f)   I=%.0f (+A55 %.0f)\n" M_ship A33 I_pitch A55
flush(stdout)

vof = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float64)
function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    # tol/itmx at construction (struct defaults from the stacked
    # WaterLily PR) — the old `pois_tol=` mom_step! kwargs are silently
    # swallowed on foam-integration (see damBreak script note).
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,),
                                tol=1e-6, itmx=50)
end
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0.0, 0.0), L_c;
    T=Float64, ν = VoF.viscosity(vof),
    g=(i, x, t) -> i == 3 ? -G_c : 0.0,
    Δt=0.25, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor, U=U∞,
)

const GRAVITY_RAMP_STEPS = NSTEPS ÷ 4
# Hold the body fixed until the impulsive-start transient (freestream
# + gravity ramp) has passed, then release from rest — the numerical
# analogue of the model-basin release. Without this the body responds
# to the unphysical startup splash and the staggered coupling spikes
# (both lightly-damped attempts diverged at step ~44, right as the
# ramp ends; V1 only survived it by being effectively over-damped).
const RELEASE_STEP = parse(Int, get(ENV, "WL_RELEASE", string(GRAVITY_RAMP_STEPS + 20)))

# One Newmark scalar DOF update. Returns (x⁺, v⁺, a⁺).
function newmark(x, v, a, F, M, C, K, dt)
    x_pred = dt * v + (0.5 - βN) * dt^2 * a       # increment predictor (Δx terms)
    v_pred = v + (1 - γN) * dt * a
    a⁺ = (F - C * v_pred - K * x_pred) / (M + γN * dt * C + βN * dt^2 * K)
    v⁺ = v_pred + γN * dt * a⁺
    x⁺ = x + x_pred + βN * dt^2 * a⁺
    return x⁺, v⁺, a⁺
end

z_hist = Float64[]; θ_hist = Float64[]
Fhz_hist = Float64[]; My_hist = Float64[]
x_CG = [hull_xc, hull_yc, hull_zc0]

t0 = time()
for step in 1:NSTEPS
    WaterLily.measure!(sim, sum(sim.flow.Δt))
    WaterLily.mom_step!(sim.flow, sim.pois)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))

    g_scale = min(1.0, step / GRAVITY_RAMP_STEPS)
    g_now = G_c * g_scale
    Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
    Mp = WaterLily.pressure_moment(x_CG, sim)
    F_buoy = -(Fp[3] + Fv[3])
    F_net  = F_buoy - M_ship * g_now
    # Absolute analytic metacentric restoring on top of the measured
    # moment (the V1 finding, reconfirmed here: at Fr=0.25 the BDIM
    # moment slope is net DESTABILIZING — Munk moment plus
    # under-captured waterplane restoring — so without an absolute
    # −K_θ·θ the pitch diverges slowly no matter how the integrator is
    # damped; the increment-implicit K below only stabilizes the
    # staggered coupling, it cannot flip an absolute slope sign).
    M_y    = -Mp[2] - ρ_w * g_now * I_wp * θ[]

    dt = sim.flow.Δt[end-1]
    # K and C from FULL gravity, not the ramped g_now: with ramped
    # values the first ~ramp steps have neither damping nor implicit
    # stiffness and the startup imbalance pumps velocity unchecked
    # (diverged at step 39 in the first attempt). Only the gravity
    # *force* ramps.
    K_z = ρ_w * G_c * A_wp
    K_θ = ρ_w * G_c * I_wp
    C_z = 2ζ * sqrt(K_z * Mz_eff)
    C_θ = 2ζ * sqrt(K_θ * Iθ_eff)

    if step > RELEASE_STEP
        z_h[], zdot_h[], zacc[] = newmark(z_h[], zdot_h[], zacc[], F_net, Mz_eff, C_z, K_z, dt)
        θ[],   θdot[],   θacc[] = newmark(θ[],   θdot[],   θacc[], M_y,  Iθ_eff, C_θ, K_θ, dt)
    end

    push!(z_hist, z_h[]); push!(θ_hist, θ[])
    push!(Fhz_hist, F_buoy); push!(My_hist, M_y)

    if step % 40 == 0
        @printf "  step=%3d  z=%+.3f ż=%+.3f  θ=%+.4f rad (%+.2f°)  F_buoy=%+.1f M_y=%+.1f  |u|=%.3f  %.1fs\n" step z_h[] zdot_h[] θ[] (θ[]*180/π) F_buoy M_y maximum(abs, sim.flow.u) (time()-t0)
        flush(stdout)
    end
    if abs(θ[]) > deg2rad(25) || abs(z_h[]) > 1.5 * T_c
        @warn "diverged at step $step: z=$(z_h[]) θ=$(θ[])"
        break
    end
end

tail = (length(z_hist) - length(z_hist)÷4 + 1):length(z_hist)
@printf "\nFinal: ⟨z⟩_tail=%+.4f cells  ⟨θ⟩_tail=%+.4f rad (%+.2f°)   max|θ|=%.2f°\n" mean(z_hist[tail]) mean(θ_hist[tail]) (mean(θ_hist[tail])*180/π) maximum(abs, θ_hist)*180/π

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "heave_pitch_newmark"))
mkpath(OUTDIR)
xs = collect(1:length(z_hist))
fig = Figure(size=(950, 700))
ax1 = Axis(fig[1, 1]; ylabel="z_h (cells)", title="Newmark-β heave+pitch (ζ=$(ζ), Ca33=$(Ca33), Ca55=$(Ca55))")
lines!(ax1, xs, z_hist; color=:steelblue, linewidth=2)
hlines!(ax1, [0]; color=:grey, linestyle=:dash)
vlines!(ax1, [GRAVITY_RAMP_STEPS]; color=:tomato, linestyle=:dash)
ax2 = Axis(fig[2, 1]; ylabel="θ (degrees)", xlabel="step")
lines!(ax2, xs, θ_hist .* (180/π); color=:purple, linewidth=2)
hlines!(ax2, [0]; color=:grey, linestyle=:dash)
save(joinpath(OUTDIR, "newmark.png"), fig)
println("Wrote $(joinpath(OUTDIR, "newmark.png"))")
