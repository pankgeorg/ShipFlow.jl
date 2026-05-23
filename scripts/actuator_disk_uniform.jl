#!/usr/bin/env julia
#
# Actuator-disk momentum-theory cross-validation (Propellers.jl Layer 1).
#
# Setup: uniform inflow at x=0, periodic in y,z (no walls); uniform-thrust
# actuator disk in the middle of the domain.  Run to steady state and
# compare the velocity at the disk plane and far wake to 1D actuator-disk
# momentum theory.
#
# Theory (Froude–Rankine, propeller convention — force in flow direction
# adds energy; this is OPPOSITE-sign of the wind-turbine convention):
#   - thrust T          = 2 ρ A U∞² a (1 + a)
#   - C_T = T / (½ρAU²) = 4 a (1 + a)
#   - induction         a = ½(-1 + √(1 + C_T))
#   - U at disk         = U∞ (1 + a)
#   - U in far wake     = U∞ (1 + 2a)
#
# Pass: induction within ±5% and far-wake velocity within ±10% of theory.

using WaterLily
using Propellers
using Printf
using Propellers: StaticArrays
const SVector = Propellers.StaticArrays.SVector

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "64"))
const R     = parse(Float64, get(ENV, "WL_R",   "8"))     # disk radius (cells)
const R_HUB = parse(Float64, get(ENV, "WL_RHB", "0"))
const Wd    = parse(Float64, get(ENV, "WL_W",   "2"))     # disk thickness (cells)
const U∞    = 1.0
const C_T_target = parse(Float64, get(ENV, "WL_CT", "0.5"))

const A_disk = π * (R^2 - R_HUB^2)
const a_theory = 0.5 * (-1 + sqrt(1 + C_T_target))   # propeller convention
const U_disk_theory = U∞ * (1 + a_theory)
const U_wake_theory = U∞ * (1 + 2*a_theory)
const thrust_target = 0.5 * 1.0 * A_disk * U∞^2 * C_T_target  # ρ=1 cell-units

@printf "=== Actuator disk in uniform stream — momentum-theory check ===\n"
@printf "  Grid          = %d × %d × %d\n" NX NY NZ
@printf "  Disk          = R=%.1f, R_hub=%.1f, w=%.1f\n" R R_HUB Wd
@printf "  U_inflow      = %.3f\n" U∞
@printf "  C_T (target)  = %.3f → a=%.4f, U_disk=%.4f, U_wake=%.4f\n" C_T_target a_theory U_disk_theory U_wake_theory
@printf "  thrust        = (1/2)·ρ·A·U²·C_T = %.4f\n\n" thrust_target

# Disk centred at x = 1/3 of domain length so we have room downstream
const cx = NX / 3
const cy = (NY + 1) / 2 - 1
const cz = (NZ + 1) / 2 - 1

disk = ActuatorDisk(
    center = SVector(cx, cy, cz),
    axis   = SVector(1.0, 0.0, 0.0),
    R      = R, R_hub = R_HUB, w = Wd,
    thrust = thrust_target,
)

flow = WaterLily.Flow((NX, NY, NZ), (Float32(U∞), 0f0, 0f0);
    T = Float32,
    ν = Float32(U∞ * 2R / 5000),     # Re_D = U·D/ν = 5000
    perdir = (2, 3),                  # periodic in y,z
    exitBC = true,                    # convective outlet in x
    Δt = 0.25,
)
pois = WaterLily.MultiLevelPoisson(flow.p, flow.μ₀, flow.σ; perdir=(2,3))

# Run to (approximate) steady state — wait T_END · NX cells for the wake
# to be flushed through.
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "600"))
const REPORT_EVERY = 50

function disk_plane_velocity(flow, cx_idx, R_cells, R_hub_cells)
    # Average u_x over the disk plane (annulus).
    u = flow.u
    s = 0.0; cnt = 0
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        x = (Float32(j) - 1.5f0 - (size(u,2)+1)/2 + 1)
        y = (Float32(k) - 1.5f0 - (size(u,3)+1)/2 + 1)
        r² = x^2 + y^2
        if R_hub_cells^2 ≤ r² ≤ R_cells^2
            s += Float64(u[cx_idx, j, k, 1])
            cnt += 1
        end
    end
    cnt == 0 ? NaN : s / cnt
end

function downstream_velocity(flow, x_idx, R_cells, R_hub_cells)
    return disk_plane_velocity(flow, x_idx, R_cells, R_hub_cells)
end

disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)
@info "Running $(NSTEPS) steps…"
for step in 1:NSTEPS
    WaterLily.mom_step!(flow, pois; udf = disk_udf)
    if mod(step, REPORT_EVERY) == 0 || step ≤ 3
        u_dk = disk_plane_velocity(flow, Int(round(cx)) + 1, R, R_HUB)
        u_5R = downstream_velocity(flow, Int(round(cx + 5R)) + 1, R, R_HUB)
        u_max = maximum(abs, flow.u)
        @info @sprintf("step=%4d  Δt=%.2e  u_disk=%.4f  u_wake_5R=%.4f  |u|_max=%.3f",
                       step, flow.Δt[end-1], u_dk, u_5R, u_max)
    end
end

# Final measurement
u_disk = disk_plane_velocity(flow, Int(round(cx)) + 1, R, R_HUB)
u_wake_5R = downstream_velocity(flow, Int(round(cx + 5R)) + 1, R, R_HUB)
u_wake_10R = downstream_velocity(flow, min(NX - 1, Int(round(cx + 10R))) + 1, R, R_HUB)

println()
@printf "  U_at_disk   measured = %.4f  vs theory %.4f  → δ = %+5.2f %%\n" (
    u_disk) U_disk_theory 100*(u_disk - U_disk_theory)/U_disk_theory
@printf "  U_at_5R     measured = %.4f  vs theory %.4f  → δ = %+5.2f %%\n" (
    u_wake_5R) U_wake_theory 100*(u_wake_5R - U_wake_theory)/U_wake_theory
@printf "  U_at_10R    measured = %.4f  (wake recovery starts beyond)\n" u_wake_10R

# Induction at disk (propeller convention): a = U_disk/U∞ - 1
a_measured = u_disk / U∞ - 1
@printf "\n  a (theory)   = %.4f\n" a_theory
@printf "  a (measured) = %.4f  → δ = %+5.2f %%\n" a_measured 100*(a_measured - a_theory)/a_theory

# Pass criteria
ok = abs(u_disk - U_disk_theory)/U_disk_theory < 0.05 &&
     abs(u_wake_5R - U_wake_theory)/U_wake_theory < 0.10
println()
println(ok ? "✓ PASS (induction within 5%, wake within 10%)" : "✗ FAIL")
exit(ok ? 0 : 1)
