#!/usr/bin/env julia
#
# Phase-3 propeller Layer-3: BLADE-RESOLVED rotating propeller in open water.
# The DTMB 4381 (unskewed 5-blade parent) is immersed as an analytic SDF
# (NavalArchitectToolbox.blade_sdf + a capped-cylinder hub) wrapped in a
# time-dependent AutoBody map that spins the shaft; WaterLily's moving-body
# BDIM (`sim_step!` remeasure) does the rest. Uniform inflow U∞ along +x,
# shaft along +x, advance ratio J = U∞/(nD).
#
# Conventions (cell units, RESULTS-propeller-layer2 lineage):
#   U∞ = 1, Δx = 1, ρ = 1;   n = U∞/(J·D)  [rev/time],  Ω = 2πn.
#   pressure_force/viscous_force return the force the BODY exerts on the
#   FLUID, so with the wake pushed downstream (+x) the propeller thrust is
#   T = +(Fp+Fv)[1] and the shaft torque is Q = +(Mp+Mv)[1].
#   KT = T/(ρ n² D⁴), KQ = Q/(ρ n² D⁵), η = J·KT/(2π·KQ).
#
# Modes:
#   WL_STATIC=1  -> measure the SDF on the grid, report immersed volume vs
#                   analytic section integral + thin-blade cell coverage,
#                   run a few frozen-rotor steps. No spinning.
#   default      -> spin WL_REVS revolutions, log KT/KQ per step to CSV,
#                   report last-rev windowed mean ± std.
#
# Run (smoke):  WL_D=64 WL_REVS=3 julia --project=. --threads=32 scripts/prop_bladeresolved_openwater.jl
#      (prod):  WL_D=128 WL_REVS=6 ...

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using WaterLily
using NavalArchitectToolbox
using NavalArchitectToolbox: blade_sdf
using StaticArrays, Printf, Statistics

const Dc    = parse(Float64, get(ENV, "WL_D", "64"))     # prop diameter [cells]
const J     = parse(Float64, get(ENV, "WL_J", "0.889"))  # advance ratio
const ReD   = parse(Float64, get(ENV, "WL_RE", "5000"))  # U∞·D/ν
const REVS  = parse(Float64, get(ENV, "WL_REVS", "3"))
const STATIC= get(ENV, "WL_STATIC", "0") == "1"
const SPIN  = parse(Float64, get(ENV, "WL_SPIN", "-1"))  # handedness: ±1
# Correct spin is NEGATIVE about the downstream shaft axis: the section
# apparent wind U ẑ − Ωr θ̂ only aligns LE→TE (chord = +cosβ θ̂ + sinβ ẑ)
# for Ω<0, and both codes confirm it empirically — at Ω>0 WaterLily and
# OpenFOAM MRF agree on REVERSED thrust (KT ≈ −0.75 vs −0.79 at J=0.889).
# (An early frozen-rotor sign read said otherwise — that was the
#  impulsive-start transient, not thrust.)
const PROP  = get(ENV, "WL_PROP", "4381")
const tab   = PROP == "4382" ? dtmb4382 : dtmb4381

const U∞ = 1.0
const Dh = 0.2Dc                    # hub diameter (r/R=0.2 table root)
const HUBL = 0.16Dc                 # hub half-length
const n  = U∞ / (J * Dc)            # rev / time-unit
const Ω  = SPIN * 2π * n            # rad / time-unit (sign = handedness)
const NX = round(Int, 4Dc); const NY = round(Int, 2Dc); const NZ = NY
const cx = 1.5Dc                    # prop plane 1.5D from inlet
const cy = NY/2; const cz = NZ/2

# ── geometry: analytic SDF in the blade frame (shaft = blade-frame z)
sdfb = blade_sdf(tab, Dc, Dh)
@inline function hub(ξ)
    q1 = hypot(ξ[1], ξ[2]) - Dh/2
    q2 = abs(ξ[3]) - HUBL
    min(max(q1, q2), 0.0) + hypot(max(q1, 0.0), max(q2, 0.0))
end
prop_sdf(ξ, t) = min(sdfb(ξ), hub(ξ))

# grid → blade frame: translate to shaft centre, un-rotate by the shaft
# angle about grid-x, then permute so blade-frame z = grid x (shaft axis).
# The spin-up RAMP is essential: an impulsive full-Ω start on sub-cell-thick
# blades fires a spurious velocity spike that crushes the CFL Δt ~100×.
const T_RAMP = parse(Float64, get(ENV, "WL_TRAMP", "0.25")) / n   # revs → t.u.
@inline shaft_angle(t) = t < T_RAMP ? Ω * t^2 / (2T_RAMP) : Ω * (t - T_RAMP/2)
@inline function prop_map(x, t)
    φ = shaft_angle(t)
    s, c = sincos(φ)
    yr = x[2] - cy; zr = x[3] - cz
    SA[c*yr + s*zr, -s*yr + c*zr, x[1] - cx]
end
body = AutoBody(prop_sdf, prop_map)

@printf("DTMB %s | D=%g cells, J=%.3f, Re_D=%g, n=%.5f rev/t, Ω=%+.5f\n",
        PROP, Dc, J, ReD, n, Ω)
@printf("domain %d×%d×%d (%.1fM cells), 1 rev = %.1f t.u.\n",
        NX, NY, NZ, NX*NY*NZ/1e6, 1/n)

sim = Simulation((NX, NY, NZ), (U∞, 0., 0.), Dc;
                 ν = U∞*Dc/ReD, body, exitBC = true)

const x₀ = SA[cx, cy, cz]
thrust_torque(s) = begin
    F = WaterLily.pressure_force(s) .+ WaterLily.viscous_force(s)  # body→fluid
    M = WaterLily.pressure_moment(x₀, s) .+ WaterLily.viscous_moment(x₀, s)
    F[1], M[1]
end
kt(T) = T / (n^2 * Dc^4); kq(Q) = Q / (n^2 * Dc^5)

OUT = joinpath(@__DIR__, "..", "runs", "prop_bladeresolved"); mkpath(OUT)

if STATIC
    # ── SDF sanity: immersed volume vs analytic thin-section integral
    a = sim.flow
    vol = 0.0; thin = 0
    σ = similar(a.p); WaterLily.measure_sdf!(σ, body, 0.0)
    vol = sum(x -> x < 0 ? 1.0 : 0.0, σ)
    # analytic: Z ∫ c(r)·t̄(r) dr per blade (t̄ ≈ 0.685·tmax NACA66-ish area
    # factor) + hub cylinder
    rs = range(Dh/Dc/2 + 1e-3, 0.5 - 1e-3; length=200) .* 2   # r/R
    dr = (rs[2]-rs[1]) * Dc/2
    Vb = sum(rs) do rR
        g = NavalArchitectToolbox.dimensional(tab, rR, Dc)
        0.685 * g.c^2 * g.tmax
    end * dr * tab.Z
    Vh = π * (Dh/2)^2 * 2HUBL
    @printf("immersed volume: %.0f cells³ | analytic est: blades %.0f + hub %.0f = %.0f (ratio %.3f)\n",
            vol, Vb, Vh, Vb+Vh, vol/(Vb+Vh))
    # thin-blade coverage: min |sdf| gradient sanity + μ₀ stats
    μ = a.μ₀
    @printf("μ₀ ∈ [%.3f, %.3f], cells with μ₀<0.5 (solid-ish): %d\n",
            minimum(μ), maximum(μ), sum(x -> x < 0.5 ? 1 : 0, μ))
    # a few frozen-rotor steps to confirm stability
    for _ in 1:20; WaterLily.mom_step!(a, sim.pois); end
    T, Q = thrust_torque(sim)
    @printf("frozen rotor after 20 steps: T=%+.3f Q=%+.3f (signs only)\n", T, Q)
else
    t_rev = 1 / n
    t_end = REVS * t_rev
    csv = joinpath(OUT, "kt_kq_$(PROP)_D$(Int(Dc))_J$(J).csv")
    open(csv, "w") do io
        println(io, "t,rev,KT,10KQ"); flush(io)
        step = 0
        while sim_time(sim) < t_end
            sim_step!(sim; remeasure = true)
            t = sim_time(sim)
            T, Q = thrust_torque(sim)
            @printf(io, "%.4f,%.3f,%.6f,%.6f\n", t, t/t_rev, kt(T), 10kq(Q))
            step += 1
            if step % 25 == 0
                flush(io)
                @printf("rev %.3f/%.1f  Δt=%.4f  KT=%.4f 10KQ=%.4f\n",
                        t/t_rev, REVS, sim.flow.Δt[end], kt(T), 10kq(Q))
            end
        end
    end
    # last-revolution windowed stats
    rows = [parse.(Float64, split(l, ',')) for l in readlines(csv)[2:end]]
    last_rev = filter(r -> r[2] > REVS - 1, rows)
    KT = mean(r[3] for r in last_rev); sKT = std(r[3] for r in last_rev)
    KQ10 = mean(r[4] for r in last_rev); sKQ = std(r[4] for r in last_rev)
    η = J * KT / (2π * KQ10/10)
    @printf("\n=== DTMB %s D=%g J=%.3f Re=%g | last rev (n=%d samples) ===\n",
            PROP, Dc, J, ReD, length(last_rev))
    @printf("KT   = %.4f ± %.4f\n10KQ = %.4f ± %.4f\nη    = %.4f\n", KT, sKT, KQ10, sKQ, η)
    println("csv: $csv")
end
