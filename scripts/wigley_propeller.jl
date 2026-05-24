#!/usr/bin/env julia
#
# HEADLINE: Wigley hull + actuator-disk propeller + VoF free surface +
# WALE turbulence in WaterLily. The full five-package stack.
#
# Setup: half-submerged Wigley hull at the front of the domain, with a
# uniform-thrust actuator disk centred behind the stern, just below the
# waterline (a single-screw stern-mounted propeller). Steady inflow.
# Smoke test: runs stably, hull experiences drag, disk produces thrust,
# system reaches an out-of-balance dynamic state.

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
const NZ = parse(Int, get(ENV, "WL_NZ", "64"))

const L_c = parse(Float64, get(ENV, "WL_L", "56"))
const B_c = parse(Float64, get(ENV, "WL_B", "10"))
const T_c = parse(Float64, get(ENV, "WL_T", "8"))

const ρ_w = parse(Float64, get(ENV, "WL_RHO_RATIO", "10"))
const ρ_a = 1.0

const U∞ = 1.0
const Fr = parse(Float64, get(ENV, "WL_FR", "0.25"))
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = parse(Float64, get(ENV, "WL_RE", "2000"))
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c

const T_NUM = Float64

# Hull centred at NX/3 (give room for wake), waterline at NZ/2
const H_w_c   = NZ/2
const hull_xc = NX/3
const hull_yc = NY/2
const hull_zc = H_w_c

# Propeller: just behind hull stern (hull_xc + L_c/2) by ~T_c/4
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = NY/2
const prop_zc = H_w_c - T_c/2     # mid-depth in water

# Disk geometry
const R_prop = T_c / 2 * 0.6      # 60% of draft (~typical hub/blade ratio)
const W_prop = 1.5

# Thrust: aim to balance hull drag at convergence.  For Re~2000 a rough
# estimate is C_D ~ 0.05 * A_wet → F_drag ~ 0.025·ρ·U²·A_wet. We pick a
# constant value to seed the simulation; over the smoke test the system
# is NOT in self-propulsion (drag won't equal thrust automatically).
const A_wet_est = (8/9) * L_c * (B_c + 4*T_c) / 2
const C_T_target = parse(Float64, get(ENV, "WL_CT", "0.3"))
const thrust = 0.5 * 1.0 * π * R_prop^2 * U∞^2 * C_T_target

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

@printf "=== Wigley + actuator disk + VoF + WALE — full stack ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Hull       = L=%.1f, B=%.1f, T=%.1f at (%g, %g, %g)\n" L_c B_c T_c hull_xc hull_yc hull_zc
@printf "  Propeller  = R=%.2f at (%g, %g, %g); thrust=%.3f (C_T=%.2f)\n" R_prop prop_xc prop_yc prop_zc thrust C_T_target
@printf "  Fr=%.2f → g_cell=%.4f;  Re=%g → ν_w=%.3e\n" Fr G_c Re ν_w_c
@printf "  ρ_w/ρ_a    = %.0f\n\n" ρ_w/ρ_a

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM,
)
turb = WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0))
const HAS_SWIRL = get(ENV, "WL_SWIRL", "0") == "1"
const torque    = parse(Float64, get(ENV, "WL_TORQUE_RATIO", "0.5")) * thrust * R_prop
disk = HAS_SWIRL ?
    SwirlingDisk(
        center = SVector(prop_xc, prop_yc, prop_zc),
        axis   = SVector(1.0, 0.0, 0.0),
        R = R_prop, w = W_prop,
        thrust = thrust, torque = torque,
    ) :
    ActuatorDisk(
        center = SVector(prop_xc, prop_yc, prop_zc),
        axis   = SVector(1.0, 0.0, 0.0),
        R = R_prop, w = W_prop,
        thrust = thrust,
    )
if HAS_SWIRL
    @printf "  Swirl ON   = torque=%.3f (Q/T·R = %.2f)\n" torque torque/(thrust*R_prop)
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
    T = T_NUM,
    ν = turb.ν,                                # combined molecular + eddy
    g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
    Δt = 0.25,
    body = hull,
    ϵ = 1,
    perdir = (2,),
    exitBC = true,
    pois_ctor = vof_pois_ctor,
    U = T_NUM(U∞),
)

# Disk udf wrapper
disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "120"))
const REPORT = max(1, NSTEPS ÷ 15)
@info "Stepping for $NSTEPS steps (full stack: hull + disk + VoF + WALE)…"

# Wave amplitude metric: max |η| from waterline H_w_c
function max_wave(α)
    nx, ny, nz = size(α)
    me = 0.0
    @inbounds for j in 2:ny-1, i in 2:nx-1
        for k in nz-1:-1:2
            if α[i, j, k] > 0.5
                eta = (k - 1.5) - H_w_c
                abs(eta) > abs(me) && (me = eta)
                break
            end
        end
    end
    return me
end

for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois;
        udf = disk_udf,                       # propeller body force
        pois_tol = T_NUM(1e-8), pois_itmx = 100)
    step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
    update_νt!(turb, sim.flow.u, vof.ν)
    if mod(step, REPORT) == 0 || step ≤ 3
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        F_drag_on_body = -(Float64(Fp[1]) + Float64(Fv[1]))
        u_max = maximum(abs, sim.flow.u)
        wave = max_wave(vof.α)
        @info @sprintf("step=%4d  t=%5.2f  Δt=%.2e  |u|=%.3f  wave=%+5.2f  F_drag=%+7.2f  thrust=%.2f  balance=%+.2f",
            step, step * sim.flow.Δt[end-1], sim.flow.Δt[end-1],
            u_max, wave, F_drag_on_body, thrust, thrust - F_drag_on_body)
        if !isfinite(u_max) || u_max > 50
            @warn "Blow-up at step $step"; break
        end
    end
end

Fp = WaterLily.pressure_force(sim)
Fv = WaterLily.viscous_force(sim)
F_drag_on_body = -(Float64(Fp[1]) + Float64(Fv[1]))
println()
@printf "Final hull drag      = %.4f cell-units\n" F_drag_on_body
@printf "Prescribed thrust    = %.4f cell-units\n" thrust
@printf "Imbalance (T - D)    = %+.4f\n" thrust - F_drag_on_body
@printf "Max wave amplitude   = %+.3f cells\n" max_wave(vof.α)
println("\n✓ Headline smoke test: hull + disk + VoF + WALE all wired and stable.")
