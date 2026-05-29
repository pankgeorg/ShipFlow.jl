#!/usr/bin/env julia
#
# Backward-facing step (Layer-3 RANS validation). Measures the
# reattachment length x_r/H on the bottom wall downstream of the step;
# compares to OpenFOAM kOmegaSST (run separately on the arm64-native
# opencfd/openfoam-default:2406 image) and the experimental x_r/H ≈ 6–7.
#
# Stabilised vs the first attempt:
#   - function inflow BC injecting U only ABOVE the step, tanh-ramped at
#     the step lip (the uniform inflow forced U into the solid block,
#     making a shear singularity that blew up);
#   - model-selectable (SA default — one equation, robust at the sharp
#     re-entrant corner; SST optional);
#   - SA uses the BDIM wall function on the bottom wall.
#
# TURB_MODEL=sa|sst   RE_H, WL_H, WL_XSTEP, WL_XOUT, WL_NSTEPS

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
]; io=devnull)

using WaterLily, Turbulence, Printf, Statistics

const MODEL  = get(ENV, "TURB_MODEL", "sa")
const H      = parse(Int,     get(ENV, "WL_H",      "20"))
const X_STEP = parse(Int,     get(ENV, "WL_XSTEP",  "40"))
const X_OUT  = parse(Int,     get(ENV, "WL_XOUT",   "260"))
const RE_H   = parse(Float64, get(ENV, "RE_H",      "5000"))
const NSTEPS = parse(Int,     get(ENV, "WL_NSTEPS", "16000"))
const U_IN   = 1.0
const T_NUM  = Float64

const N_Y = 2H
const N_X = X_OUT
const NU  = U_IN * H / RE_H

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "bfs_$(MODEL)"))
mkpath(OUTDIR)

@printf("BFS-%s: %d×%d  H=%d  x_step=%d  ν=%.3e  Re_H=%.0f\n",
        uppercase(MODEL), N_X, N_Y, H, X_STEP, NU, RE_H)

# Floor profile and body (step block + top wall). Fluid: floor(x) < y < N_Y.
floor_y(x) = x < X_STEP ? T_NUM(H) : zero(T_NUM)
body = AutoBody((x, t) -> min(x[2] - floor_y(x[1]), N_Y - x[2]))

# Smooth inflow profile: U above the step, 0 below, tanh-ramped over ~2
# cells at the step lip so there is no discontinuity feeding the solver.
ramp(y) = (one(T_NUM) + tanh((y - H) / T_NUM(2))) / 2
inflow_u(y) = T_NUM(U_IN) * ramp(y)
uBC = (i, x, t) -> i == 1 ? inflow_u(x[2]) : zero(T_NUM)
uλ  = (i, x)    -> i == 1 ? inflow_u(x[2]) : zero(T_NUM)

# Advection limiter (WL_LAMBDA=quick|vanLeer|cds) for k/ω + momentum.
const LAM = let s = lowercase(get(ENV, "WL_LAMBDA", "quick"))
    s == "vanleer" ? WaterLily.vanLeer : s == "cds" ? WaterLily.cds : WaterLily.quick
end

if MODEL == "sst"
    k_in = 1.5 * (0.05*U_IN)^2
    ω_in = sqrt(k_in) / (0.09^0.25 * 0.07*H)
    model = KOmegaSST((N_X, N_Y), body; ν=NU, k∞=k_in, ω∞=ω_in, T=T_NUM)
    stepturb!(u, dt) = step_sst!(model, u, dt; λ=LAM)          # native ω-wall
else
    model = SpalartAllmaras((N_X, N_Y), body; ν=NU, ν̃∞=3NU, T=T_NUM)
    stepturb!(u, dt) = step_sa!(model, u, dt; wallfn=true, band=(T_NUM(1), T_NUM(4)))
end

sim = Simulation((N_X, N_Y), uBC, Float64(H);
    U=U_IN, uλ=uλ, ν=model.ν, body=body, exitBC=true, T=T_NUM)

# Reattachment length = the rightmost near-wall flow reversal of the
# main bubble: the largest x (within the bubble window) where the first
# interior cell above the bottom wall has u_x < 0. Taking the *last*
# reversal (not the first) skips the small counter-rotating corner eddy.
function reattachment(u)
    jw = 2; last_neg = 0
    for i in (X_STEP+2):min(X_STEP + 15H, N_X-2)
        u[i, jw, 1] < 0 && (last_neg = i)
    end
    return last_neg == 0 ? NaN : (last_neg - X_STEP)
end

for step in 1:NSTEPS
    stepturb!(sim.flow.u, sim.flow.Δt[end])
    sim_step!(sim; remeasure=false, λ=LAM)
    if step % 2000 == 0
        um = maximum(abs, sim.flow.u); xr = reattachment(sim.flow.u)
        @printf("  step %5d  |u|max=%.3f  νt/ν_max=%.1f  x_r/H=%s\n",
                step, um, (maximum(model.ν)-NU)/NU,
                isnan(xr) ? "—" : @sprintf("%.2f", xr/H))
        flush(stdout)
        isfinite(um) || (println("  DIVERGED"); break)
    end
end

xr = reattachment(sim.flow.u)
@printf("\n=== BFS-%s summary ===\n", uppercase(MODEL))
@printf("  x_r/H = %s   (exp ≈ 6–7)\n", isnan(xr) ? "not found" : @sprintf("%.2f", xr/H))
open(joinpath(OUTDIR, "wall_u.csv"), "w") do io
    println(io, "x_over_H,u_wall")
    for i in (X_STEP+2):(N_X-2)
        @printf(io, "%.3f,%.5f\n", (i-X_STEP)/H, sim.flow.u[i, 2, 1])
    end
end
@printf("  wrote %s\n=======================\n", joinpath(OUTDIR, "wall_u.csv"))
