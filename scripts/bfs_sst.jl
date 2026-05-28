#!/usr/bin/env julia
#
# Backward-facing step (Driver–Seegmiller 1985) — Turbulence.jl k–ω SST
# Layer-3 validation. Measures the reattachment length x_r/H on the
# bottom wall downstream of the step and compares to the experimental
# x_r/H ≈ 6.0–6.3 (kOmegaSST literature typically 6–7).
#
# OpenFOAM can't run on this (arm64) host — the amd64 OF image fails
# exec — so we validate against the published experiment, not a fresh
# OF run.
#
# Geometry (expansion ratio 2): inlet channel of height H above a step
# of height H; the floor drops from y=H to y=0 at x=x_step. Top is a
# wall. Uniform inflow U; BDIM walls; SST with its native ω-wall
# treatment (no Spalding override — channel showed it double-counts).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
]; io=devnull)

using WaterLily, Turbulence, Printf, Statistics

const H       = parse(Int,     get(ENV, "WL_H",      "24"))    # step height (cells)
const X_STEP  = parse(Int,     get(ENV, "WL_XSTEP",  "48"))    # step location
const X_OUT   = parse(Int,     get(ENV, "WL_XOUT",   "336"))   # domain length
const RE_H    = parse(Float64, get(ENV, "RE_H",      "37500")) # Re on step height
const NSTEPS  = parse(Int,     get(ENV, "WL_NSTEPS", "20000"))
const U_IN    = 1.0
const T_NUM   = Float64

const N_Y = 2H                     # full channel height (expansion ratio 2)
const N_X = X_OUT
const NU  = U_IN * H / RE_H

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "bfs_sst"))
mkpath(OUTDIR)

@printf("BFS-SST: %d×%d  H=%d  x_step=%d  ν=%.3e  Re_H=%.0f\n",
        N_X, N_Y, H, X_STEP, NU, RE_H)

# Floor profile: y=H upstream of the step, y=0 downstream. Body is the
# step block + the top wall: fluid is floor(x) < y < N_Y.
floor_y(x) = x < X_STEP ? T_NUM(H) : zero(T_NUM)
body = AutoBody((x, t) -> min(x[2] - floor_y(x[1]), N_Y - x[2]))

# SST: inlet turbulence intensity ~5%, ω from a mixing-length estimate.
k_in = 1.5 * (0.05*U_IN)^2
ω_in = sqrt(k_in) / (0.09^0.25 * 0.07*H)
sst = KOmegaSST((N_X, N_Y), body; ν=NU, k∞=k_in, ω∞=ω_in, T=T_NUM)

# IC: uniform U above the step, 0 inside it.
uλ = (i, x) -> (i == 1 && x[2] > floor_y(x[1])) ? T_NUM(U_IN) : zero(T_NUM)

sim = Simulation((N_X, N_Y), (U_IN, 0.0), Float64(H);
    U=U_IN, uλ=uλ, ν=sst.ν, body=body, exitBC=true, T=T_NUM)

# Near-wall streamwise velocity along the bottom wall (first fluid cell
# above y=0), downstream of the step. Sign change neg→pos = reattachment.
function reattachment_x(u)
    jwall = 2                       # first interior cell above y=0
    prev_neg = false; xr = NaN
    for i in (X_STEP+2):(N_X-2)
        uw = u[i, jwall, 1]
        if uw < 0
            prev_neg = true
        elseif prev_neg && uw ≥ 0
            xr = i - X_STEP          # cells past the step
            break
        end
    end
    return xr
end

for step in 1:NSTEPS
    step_sst!(sst, sim.flow.u, sim.flow.Δt[end])   # no wall fn (native ω-wall)
    sim_step!(sim; remeasure=false)
    if step % 2000 == 0
        xr = reattachment_x(sim.flow.u)
        @printf("  step %5d  |u|max=%.3f  x_r/H=%s\n",
                step, maximum(abs, sim.flow.u),
                isnan(xr) ? "—" : @sprintf("%.2f", xr/H))
        flush(stdout)
    end
end

xr = reattachment_x(sim.flow.u)
@printf("\n=== BFS-SST summary ===\n")
@printf("  reattachment x_r = %s cells\n", isnan(xr) ? "not found" : @sprintf("%.1f", xr))
@printf("  x_r/H = %s   (Driver–Seegmiller exp ≈ 6.0–6.3; kΩSST lit 6–7)\n",
        isnan(xr) ? "—" : @sprintf("%.2f", xr/H))
@printf("=======================\n")

# Dump the bottom-wall u profile for inspection.
open(joinpath(OUTDIR, "wall_u.csv"), "w") do io
    println(io, "x_over_H,u_wall")
    for i in (X_STEP+2):(N_X-2)
        @printf(io, "%.3f,%.5f\n", (i-X_STEP)/H, sim.flow.u[i, 2, 1])
    end
end
