#!/usr/bin/env julia
#
# Spalart–Allmaras RANS channel at Re_τ ≈ 395 — Turbulence.jl validation
# against the law of the wall (the canonical SA / DNS target; OpenFOAM's
# SA reproduces the same log law).
#
# RANS, unlike LES, models *all* turbulence via ν_t — so no transition
# trick is needed. A steady 2D channel (periodic streamwise, BDIM walls
# in y) develops ν_t(y) and the turbulent mean profile u(y) directly.
#
# Segregated loop:  step_sa!(model, u, dt)  then  sim_step!.
#
# Pass: u⁺(y⁺) matches u⁺ = (1/0.41)·ln y⁺ + 5.2 in the log layer
# (30 < y⁺ < 0.3·Re_τ) to within ±10%.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
]; io=devnull)

using WaterLily, Turbulence, Printf, Statistics

const N_HC   = parse(Int,     get(ENV, "WL_N_HC",  "64"))   # cells per half-height δ
const N_X    = parse(Int,     get(ENV, "N_X",      "32"))   # streamwise (periodic)
const RE_TAU = parse(Float64, get(ENV, "RE_TAU",   "395"))
const NSTEPS = parse(Int,     get(ENV, "WL_NSTEPS","6000"))
const WALLFN = parse(Bool,    get(ENV, "WL_WALLFN","false"))  # BDIM wall function
const BAND_LO = parse(Float64, get(ENV, "WL_BAND_LO", "1.0"))
const BAND_HI = parse(Float64, get(ENV, "WL_BAND_HI", "3.0"))
const U_TAU  = 1.0                                          # work in wall units
const T_NUM  = Float64

const N_Y   = 2 * N_HC
const DELTA = Float64(N_HC)
const NU    = U_TAU * DELTA / RE_TAU          # ν from Re_τ definition
const G_X   = U_TAU^2 / DELTA                 # constant streamwise driver −dP/dx

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "channel_sa"))
mkpath(OUTDIR)

@printf("SA channel: N=%d×%d  δ=%.0f  ν=%.4e  Re_τ=%.0f  g_x=%.4e\n",
        N_X, N_Y, DELTA, NU, RE_TAU, G_X)

# Two-plate channel SDF (positive inside fluid).
body = AutoBody((x, t) -> min(x[2], N_Y - x[2]))

# SA model: ν̃∞ = 3ν as a turbulent seed.
sa = SpalartAllmaras((N_X, N_Y), body; ν=NU, ν̃∞=3NU, perdir=(1,), T=T_NUM)

# IC: a modest plug-ish profile; RANS evolves it to the turbulent mean.
uλ = (i, x) -> begin
    if i == 1
        η = (x[2] - N_HC) / N_HC          # −1..1 across channel
        T_NUM(18 * U_TAU * (1 - η^2))     # ~bulk/u_τ ≈ 16.9 scale
    else
        zero(T_NUM)
    end
end

g = (i, x, t) -> i == 1 ? T_NUM(G_X) : zero(T_NUM)

sim = Simulation((N_X, N_Y), (0.0, 0.0), DELTA;
    U=U_TAU, uλ=uλ, ν=sa.ν, body=body, g=g, perdir=(1,), T=T_NUM)

# ── time loop: SA substep first (sets ν), then momentum ──
for step in 1:NSTEPS
    step_sa!(sa, sim.flow.u, sim.flow.Δt[end]; wallfn=WALLFN, band=(BAND_LO, BAND_HI))
    sim_step!(sim; remeasure=false)
    if step % 500 == 0
        umax = maximum(@view sim.flow.u[:, :, 1])
        νtmax = (maximum(sa.ν) - NU) / NU
        @printf("  step %5d  u_max=%.3f  νt/ν_max=%.2f  Δt=%.3e\n",
                step, umax, νtmax, sim.flow.Δt[end])
        flush(stdout)
    end
end

# ── sample mean profile u(y), averaged over x ──
function ux_of_y(u)
    nx, ny, _ = size(u)
    [mean(@view u[2:nx-1, j, 1]) for j in 2:ny-1]
end
prof = ux_of_y(sim.flow.u)

# u_τ from the exact channel force balance: at steady state the net
# body force equals the wall drag, τ_w = g_x·δ ⇒ u_τ = √(g_x·δ). This is
# the physically correct friction velocity (the near-wall slope estimate
# is unreliable under the smeared BDIM wall).
uτ_fb = sqrt(G_X * DELTA)
# Near-wall-slope estimate, for reference / to expose the BDIM gap.
uτ_slope = sqrt(NU * prof[1] / 0.5)

@printf("\n=== SA channel summary ===\n")
@printf("  u_τ (force balance) = %.4f  (target %.4f)\n", uτ_fb, U_TAU)
@printf("  u_τ (wall slope)    = %.4f  (BDIM-smeared, underestimates)\n", uτ_slope)
const uτ_meas = uτ_fb
const Reτ_meas = uτ_meas * DELTA / NU
@printf("  Re_τ (force bal.)   = %.1f  (target %.1f)\n", Reτ_meas, RE_TAU)
@printf("  centreline u⁺       = %.2f  (DNS ≈ 20.5)\n", maximum(prof)/uτ_meas)

# u⁺(y⁺) and log-law comparison.
open(joinpath(OUTDIR, "uplus.csv"), "w") do io
    println(io, "y_plus,u_plus,loglaw")
    err_acc = Float64[]
    for (jj, u) in pairs(prof)
        y = jj - 0.5                         # distance from wall (cells)
        y > DELTA && continue                # lower half only
        yplus = y * uτ_meas / NU
        uplus = u / uτ_meas
        loglaw = (1/0.41) * log(yplus) + 5.2
        @printf(io, "%.4f,%.4f,%.4f\n", yplus, uplus, loglaw)
        if 30 < yplus < 0.3*RE_TAU
            push!(err_acc, abs(uplus - loglaw) / loglaw)
        end
    end
    if !isempty(err_acc)
        @printf("  log-layer mean rel. error vs law of the wall = %.1f%%\n",
                100*mean(err_acc))
        @printf("  log-layer max  rel. error = %.1f%%\n", 100*maximum(err_acc))
    else
        println("  (no points in 30 < y⁺ < 0.3 Re_τ — check resolution)")
    end
end
@printf("  wrote %s\n", joinpath(OUTDIR, "uplus.csv"))
@printf("==========================\n")
