#!/usr/bin/env julia
#
# k–ω SST RANS channel at Re_τ ≈ 395 — Turbulence.jl validation against
# the law of the wall (OpenFOAM's kOmegaSST reproduces the same log law).
#
# Mirror of channel_sa.jl: steady 2D channel (periodic streamwise, BDIM
# walls), segregated loop step_sst!(model,u,dt) → sim_step!.
#
# Two modes (WL_WALLFN): without the BDIM wall function we expect the
# same B-offset as SA (confirming the wall, not the closure, is the
# limit); with it, the outer log layer should recover.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
]; io=devnull)

using WaterLily, Turbulence, Printf, Statistics

const N_HC   = parse(Int,     get(ENV, "WL_N_HC",  "64"))
const N_X    = parse(Int,     get(ENV, "N_X",      "32"))
const RE_TAU = parse(Float64, get(ENV, "RE_TAU",   "395"))
const NSTEPS = parse(Int,     get(ENV, "WL_NSTEPS","25000"))
const WALLFN = parse(Bool,    get(ENV, "WL_WALLFN","false"))
const BAND_LO= parse(Float64, get(ENV, "WL_BAND_LO","4.0"))
const BAND_HI= parse(Float64, get(ENV, "WL_BAND_HI","12.0"))
const U_TAU  = 1.0
const T_NUM  = Float64

const N_Y   = 2 * N_HC
const DELTA = Float64(N_HC)
const NU    = U_TAU * DELTA / RE_TAU
const G_X   = U_TAU^2 / DELTA

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "channel_sst"))
mkpath(OUTDIR)

@printf("SST channel: N=%d×%d  δ=%.0f  ν=%.4e  Re_τ=%.0f  wallfn=%s band=(%.0f,%.0f)\n",
        N_X, N_Y, DELTA, NU, RE_TAU, WALLFN, BAND_LO, BAND_HI)

body = AutoBody((x, t) -> min(x[2], N_Y - x[2]))

# SST seeds: k∞ ≈ (turb intensity·U)², ω∞ ≈ U/δ scale.
sst = KOmegaSST((N_X, N_Y), body; ν=NU, k∞=1e-3, ω∞=1.0, perdir=(1,), T=T_NUM)

uλ = (i, x) -> i == 1 ? T_NUM(18 * U_TAU * (1 - ((x[2]-N_HC)/N_HC)^2)) : zero(T_NUM)
g  = (i, x, t) -> i == 1 ? T_NUM(G_X) : zero(T_NUM)

sim = Simulation((N_X, N_Y), (0.0, 0.0), DELTA;
    U=U_TAU, uλ=uλ, ν=sst.ν, body=body, g=g, perdir=(1,), T=T_NUM)

for step in 1:NSTEPS
    step_sst!(sst, sim.flow.u, sim.flow.Δt[end]; wallfn=WALLFN, band=(BAND_LO, BAND_HI))
    sim_step!(sim; remeasure=false)
    if step % 2500 == 0
        umax = maximum(@view sim.flow.u[:, :, 1])
        νtmax = (maximum(sst.ν) - NU) / NU
        @printf("  step %5d  u_max=%.3f  νt/ν_max=%.2f  k_max=%.3e  ω_max=%.2e\n",
                step, umax, νtmax, maximum(sst.k), maximum(sst.ω))
        flush(stdout)
    end
end

ux_of_y(u) = (nx=size(u,1); ny=size(u,2); [mean(@view u[2:nx-1, j, 1]) for j in 2:ny-1])
prof = ux_of_y(sim.flow.u)

uτ_fb = sqrt(G_X * DELTA)
@printf("\n=== SST channel summary ===\n")
@printf("  u_τ (force balance) = %.4f\n", uτ_fb)
@printf("  Re_τ (force bal.)   = %.1f  (target %.1f)\n", uτ_fb*DELTA/NU, RE_TAU)
@printf("  centreline u⁺       = %.2f  (DNS ≈ 20.5)\n", maximum(prof)/uτ_fb)

open(joinpath(OUTDIR, "uplus.csv"), "w") do io
    println(io, "y_plus,u_plus,loglaw")
    err = Float64[]
    for (jj, u) in pairs(prof)
        y = jj - 0.5
        y > DELTA && continue
        yplus = y * uτ_fb / NU
        uplus = u / uτ_fb
        loglaw = (1/0.41)*log(yplus) + 5.2
        @printf(io, "%.4f,%.4f,%.4f\n", yplus, uplus, loglaw)
        30 < yplus < 0.3*RE_TAU && push!(err, abs(uplus-loglaw)/loglaw)
    end
    isempty(err) || @printf("  log-layer mean rel. error = %.1f%%   max = %.1f%%\n",
                            100*mean(err), 100*maximum(err))
end
@printf("  wrote %s\n==========================\n", joinpath(OUTDIR, "uplus.csv"))
