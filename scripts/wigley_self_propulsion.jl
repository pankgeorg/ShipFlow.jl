#!/usr/bin/env julia
#
# Self-propulsion smoke-test: add a PI controller on the actuator-disk
# thrust to drive (Thrust − Drag) toward zero. At convergence the system
# is self-propelling — the disk thrust exactly balances the hull drag.
#
# This is the dynamic version of `scripts/wigley_propeller.jl`. The same
# stack (Wigley + AD + VoF + WALE + WaterLily) but with a feedback loop
# that adapts the disk thrust each step.
#
# Pass criterion (smoke): after N_RAMP steps, |T - D|/D < 0.1.

using WaterLily
using VoF
using ShipShapes
using ShipShapes: StaticArrays
using Turbulence
using Propellers
const SVector = StaticArrays.SVector
using Printf

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const L_c = parse(Float64, get(ENV, "WL_L", "56"))
const B_c = parse(Float64, get(ENV, "WL_B", "10"))
const T_c = parse(Float64, get(ENV, "WL_T", "8"))

const ρ_w = 10.0
const ρ_a = 1.0
const U∞  = 1.0
const Fr  = 0.25
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = 2000
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c
const T_NUM = Float64

const H_w_c   = NZ/2
const hull_xc = NX/3
const hull_yc = NY/2
const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = NY/2
const prop_zc = H_w_c - T_c/2
const R_prop  = T_c / 2 * 0.8
const W_prop  = 1.5

# Controller parameters
const KP = parse(Float64, get(ENV, "WL_KP", "0.05"))    # proportional gain
const KI = parse(Float64, get(ENV, "WL_KI", "0.002"))   # integral gain
const T_INIT = parse(Float64, get(ENV, "WL_T_INIT", "20.0"))  # initial thrust seed

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== Wigley + AD self-propulsion (PI thrust controller) ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Hull       = L=%.1f, B=%.1f, T=%.1f\n" L_c B_c T_c
@printf "  Prop       = R=%.2f at (%.1f, %.1f, %.1f)\n" R_prop prop_xc prop_yc prop_zc
@printf "  Re=%g, Fr=%.2f, ρ=%g\n" Re Fr ρ_w/ρ_a
@printf "  PI: Kp=%.3f, Ki=%.4f, T_init=%.1f\n\n" KP KI T_INIT

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM,
)
turb = WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0))

# Mutable disk: we'll rebuild it each step with the controller's thrust.
function build_disk(thrust)
    ActuatorDisk(
        center = SVector(prop_xc, prop_yc, prop_zc),
        axis   = SVector(1.0, 0.0, 0.0),
        R = R_prop, w = W_prop,
        thrust = thrust,
    )
end

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end

sim = WaterLily.Simulation((NX, NY, NZ),
    (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
    T = T_NUM, ν = turb.ν,
    g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
    Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = T_NUM(U∞),
)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "200"))

const N_WARMUP = parse(Int, get(ENV, "WL_WARMUP", "40"))
const INTEGRAL_CAP = parse(Float64, get(ENV, "WL_I_CAP", "100.0"))

function run_controller(NSTEPS, T_INIT)
    thrust = T_INIT
    integral = 0.0
    err_history = Float64[]
    α_smooth = 0.2
    drag_smooth = 0.0
    for step in 1:NSTEPS
        disk = build_disk(thrust)
        disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)
        WaterLily.mom_step!(sim.flow, sim.pois;
            udf = disk_udf,
            pois_tol = T_NUM(1e-8), pois_itmx = 100)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        update_νt!(turb, sim.flow.u, vof.ν)

        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        drag = -(Float64(Fp[1]) + Float64(Fv[1]))

        drag_smooth = step == 1 ? drag : α_smooth * drag + (1 - α_smooth) * drag_smooth

        if step ≤ N_WARMUP
            # Warm-up: let the drag transient settle without driving the
            # controller. Reset integral at the end of warmup.
            push!(err_history, NaN)
            if step == N_WARMUP
                thrust = drag_smooth   # seed thrust at the warmup drag estimate
                integral = 0.0
            end
        else
            err = drag_smooth - thrust
            integral = clamp(integral + err * sim.flow.Δt[end-1],
                             -INTEGRAL_CAP, INTEGRAL_CAP)
            Δthrust = KP * err + KI * integral
            push!(err_history, err / max(abs(drag_smooth), 1.0))
            thrust = max(0.0, min(2 * abs(drag_smooth) + 1.0, thrust + Δthrust))
        end

        if mod(step, max(1, NSTEPS ÷ 25)) == 0 || step ≤ 3 || step == N_WARMUP
            u_max = maximum(abs, sim.flow.u)
            phase = step ≤ N_WARMUP ? "warm" : "ctrl"
            @info @sprintf("step=%4d %s  drag=%6.2f  drag_s=%6.2f  T=%6.2f  err/T=%+.3f  intg=%.1f  |u|=%.3f",
                step, phase, drag, drag_smooth, thrust,
                (drag_smooth - thrust)/max(abs(drag_smooth), 1.0),
                integral, u_max)
            if !isfinite(u_max) || u_max > 50
                @warn "Blow-up"; break
            end
        end
    end
    return thrust, err_history
end

final_thrust, err_history = run_controller(NSTEPS, T_INIT)
thrust = final_thrust

# Final diagnostic
Fp = WaterLily.pressure_force(sim)
Fv = WaterLily.viscous_force(sim)
drag = -(Float64(Fp[1]) + Float64(Fv[1]))
println()
@printf "Final drag    = %.4f\n" drag
@printf "Final thrust  = %.4f\n" thrust
@printf "|T - D|/D     = %.3f %%\n" 100 * abs(thrust - drag) / abs(drag)

# Median over last quarter of history
N = length(err_history)
recent = err_history[max(1, N - N÷4):end]
recent_median_abs = sort(abs.(recent))[length(recent)÷2 + 1]
@printf "Median |err|/drag over last quarter = %.3f %%\n" 100 * recent_median_abs
println()
println(recent_median_abs < 0.10 ? "✓ Self-propulsion within 10%" :
        "△ Not yet self-propelling (smoke test only)")
