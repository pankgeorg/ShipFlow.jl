#!/usr/bin/env julia
#
# Ladder 1 (Phase-3 propeller Layer-1.5): quantitative actuator-disk
# momentum-theory curve. Extend the single-point spot check
# (actuator_disk_uniform.jl / RESULTS-propeller.md) into a C_T sweep,
# comparing the in-grid disk's induced velocity and slipstream
# contraction to 1-D Froude-Rankine actuator-disk theory.
#
# Propeller (energy-adding) convention:
#   C_T = T/(½ρAU∞²) = 4 a (1 + a)            a = ½(−1 + √(1 + C_T))
#   U_disk = U∞(1 + a)                        U_wake = U∞(1 + 2a)
#   slipstream contraction:  A_wake/A_disk = (1+a)/(1+2a)
#                            R_wake/R_disk = √((1+a)/(1+2a))
#
# Run: a uniform stream, periodic y/z, convective exit, top-hat disk.
# Sweep C_T; for each, record U_disk, U_wake(5R), and the measured
# slipstream radius (from the axial-velocity excess contour), then
# compare to theory. Writes a CSV + a summary table.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using WaterLily
using Propellers
using Printf
using Statistics
using Propellers: StaticArrays
const SVector = Propellers.StaticArrays.SVector

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "96"))
const NZ = parse(Int, get(ENV, "WL_NZ", "96"))
const R     = parse(Float64, get(ENV, "WL_R",   "12"))   # disk radius (cells)
const Wd    = parse(Float64, get(ENV, "WL_W",   "3"))    # disk thickness (cells)
const U∞    = 1.0
const ReD   = parse(Float64, get(ENV, "WL_RE",  "5000"))
const NSTEPS= parse(Int, get(ENV, "WL_NSTEPS", "800"))

# Disk centred 1/3 along x so there is ≥ several D of downstream domain.
const cx = NX / 3
const cy = (NY + 1) / 2 - 1
const cz = (NZ + 1) / 2 - 1
const A_disk = π * R^2

const CT_LIST = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]

theory_a(CT)     = 0.5 * (-1 + sqrt(1 + CT))
theory_Udisk(CT) = U∞ * (1 + theory_a(CT))
theory_Uwake(CT) = U∞ * (1 + 2*theory_a(CT))
theory_Rwake(CT) = R * sqrt((1 + theory_a(CT)) / (1 + 2*theory_a(CT)))

# Mean u_x over the annulus (r ≤ R) at axial index x_idx.
function disk_plane_velocity(flow, x_idx, R_cells)
    u = flow.u
    s = 0.0; cnt = 0
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        y = Float64(j) - 1.5 - cy
        z = Float64(k) - 1.5 - cz
        if y^2 + z^2 ≤ R_cells^2
            s += Float64(u[x_idx, j, k, 1]); cnt += 1
        end
    end
    cnt == 0 ? NaN : s/cnt
end

# Slipstream radius at axial index x_idx: the radius enclosing the half-
# excess axial velocity. We bin u_x by radius, find the centreline excess
# (U_c − U∞), and report the radius where the excess falls to half.
function slipstream_radius(flow, x_idx)
    u = flow.u
    rmax = min(cy, cz) - 1
    nb = 40
    edges = range(0, rmax; length=nb+1)
    sums = zeros(nb); cnts = zeros(Int, nb)
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        y = Float64(j) - 1.5 - cy
        z = Float64(k) - 1.5 - cz
        r = sqrt(y^2 + z^2)
        b = searchsortedlast(edges, r)
        (1 ≤ b ≤ nb) || continue
        sums[b] += Float64(u[x_idx, j, k, 1]); cnts[b] += 1
    end
    centers = [0.5*(edges[i]+edges[i+1]) for i in 1:nb]
    prof = [cnts[i] > 0 ? sums[i]/cnts[i] : NaN for i in 1:nb]
    # centreline excess from the innermost valid bins
    valid = findall(!isnan, prof)
    isempty(valid) && return NaN
    Uc = mean(prof[valid[1:min(3,length(valid))]])
    half = U∞ + 0.5*(Uc - U∞)
    # walk outward to where the profile drops through `half`
    for i in valid
        if prof[i] ≤ half && centers[i] > 0.4R
            return centers[i]
        end
    end
    return NaN
end

run_one(CT) = begin
    thrust = 0.5 * 1.0 * A_disk * U∞^2 * CT      # ρ=1 cell-units
    disk = ActuatorDisk(center=SVector(cx,cy,cz), axis=SVector(1.0,0.0,0.0),
                        R=R, w=Wd, thrust=thrust)
    flow = WaterLily.Flow((NX,NY,NZ), (Float32(U∞),0f0,0f0);
        T=Float32, ν=Float32(U∞*2R/ReD), perdir=(2,3), exitBC=true, Δt=0.25)
    pois = WaterLily.MultiLevelPoisson(flow.p, flow.μ₀, flow.σ; perdir=(2,3))
    udf = (fl,t;kw...) -> disk(fl,t;kw...)
    for _ in 1:NSTEPS
        WaterLily.mom_step!(flow, pois; udf=udf)
    end
    xd = Int(round(cx)) + 1
    x5 = min(NX-1, Int(round(cx + 5R)) + 1)
    Ud = disk_plane_velocity(flow, xd, R)
    Uw = disk_plane_velocity(flow, x5, R)
    Rw = slipstream_radius(flow, x5)
    (; CT, Ud, Uw, Rw)
end

@printf "=== Ladder 1: actuator-disk C_T curve vs Froude-Rankine theory ===\n"
@printf "  Grid %d×%d×%d  R=%.0f w=%.0f  Re_D=%.0f  steps=%d\n\n" NX NY NZ R Wd ReD NSTEPS
flush(stdout)

OUT = joinpath(@__DIR__, "..", "runs", "actuator_disk_ct")
mkpath(OUT)
csv = joinpath(OUT, "ct_curve.csv")
open(csv, "w") do io
    println(io, "CT,a_theory,Udisk_theory,Udisk_meas,err_Udisk_pct,",
                "Uwake_theory,Uwake_meas,err_Uwake_pct,",
                "Rwake_theory,Rwake_meas,err_Rwake_pct")
end

@printf "%5s | %8s %8s %7s | %8s %8s %7s | %8s %8s %7s\n" "C_T" "Ud_th" "Ud_ms" "δ%%" "Uw_th" "Uw_ms" "δ%%" "Rw_th" "Rw_ms" "δ%%"
for CT in CT_LIST
    r = run_one(CT)
    a = theory_a(CT)
    Udt = theory_Udisk(CT); Uwt = theory_Uwake(CT); Rwt = theory_Rwake(CT)
    eUd = 100*(r.Ud-Udt)/Udt
    eUw = 100*(r.Uw-Uwt)/Uwt
    eRw = isnan(r.Rw) ? NaN : 100*(r.Rw-Rwt)/Rwt
    @printf "%5.2f | %8.4f %8.4f %+6.1f | %8.4f %8.4f %+6.1f | %8.4f %8.4f %+6.1f\n" CT Udt r.Ud eUd Uwt r.Uw eUw Rwt r.Rw eRw
    open(csv, "a") do io
        @printf io "%.3f,%.5f,%.5f,%.5f,%.3f,%.5f,%.5f,%.3f,%.5f,%.5f,%.3f\n" CT a Udt r.Ud eUd Uwt r.Uw eUw Rwt r.Rw eRw
    end
    flush(stdout)
end
@printf "\nWrote %s\n" csv
