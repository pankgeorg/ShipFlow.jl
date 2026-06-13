#!/usr/bin/env julia
#
# Ladder 2 (Phase-3 propeller Layer-2 core): VLM → in-grid disk →
# integrated K_T/K_Q/η round-trip for DTMB 4382.
#
# Pipeline:
#   1. NavalArchitectToolbox open-water VLM (validated vs experiment ~2%)
#      gives KT, KQ, η AND a radial loading shape dT/dr, dQ/dr vs r/R
#      (dumped to runs/radial_loading_J*.csv by
#       NavalArchitectToolbox.jl/examples/dump_radial_loading.jl).
#   2. Build a Propellers.GradedDisk carrying that radial shape, with the
#      total thrust/torque set from the VLM coefficients in cell units.
#   3. Run it in a uniform WaterLily stream (single phase, periodic y/z,
#      convective exit).
#   4. Integrate the RESOLVED thrust & torque the fluid carries away
#      (axial momentum-flux balance across the box; swirl angular-momentum
#      flux) and recompute KT/KQ/η. Check round-trip vs VLM within ±10%.
#
# Cell units. ρ=1, U∞=1 (cell/Δt). Choose the in-grid diameter D (cells)
# so R = D/2 ≥ ~10 cells; the rotation rate follows from the advance
# ratio: J = U∞/(n·D) ⇒ n = U∞/(J·D). Coefficients use this same (n,D).
#   KT = T/(ρ n² D⁴),  KQ = Q/(ρ n² D⁵),  η = J·KT/(2π·KQ).

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using WaterLily
using Propellers
using Printf
using Statistics
using Propellers: StaticArrays
const SVector = Propellers.StaticArrays.SVector

const J      = parse(Float64, get(ENV, "WL_J", "0.889"))
const Rcells = parse(Float64, get(ENV, "WL_R", "12"))
const Dcell  = 2*Rcells
const Wd     = parse(Float64, get(ENV, "WL_W", "3"))
const U∞     = 1.0
const ReD    = parse(Float64, get(ENV, "WL_RE", "5000"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "1000"))

# advance ratio → rotation rate (rev per Δt) in cell units
const nRPS = U∞ / (J * Dcell)

# domain: long in x for the slipstream, periodic in y/z
const NX = parse(Int, get(ENV, "WL_NX", "240"))
const NY = parse(Int, get(ENV, "WL_NY", "96"))
const NZ = parse(Int, get(ENV, "WL_NZ", "96"))
const cx = NX/4
const cy = (NY+1)/2 - 1
const cz = (NZ+1)/2 - 1

const LOADDIR = abspath(joinpath(@__DIR__, "..", "..", "NavalArchitectToolbox.jl", "runs"))

function read_loading(J)
    fn = joinpath(LOADDIR, @sprintf("radial_loading_J%.3f.csv", J))
    isfile(fn) || error("missing loading CSV $fn — run dump_radial_loading.jl first")
    hdr = readline(fn)
    # parse KT/KQ/eta from header "# ... KT=.. KQ=.. eta=.. ..."
    getf(key) = (m = match(Regex(key*"=([-0-9.eE]+)"), hdr); parse(Float64, m.captures[1]))
    KT = getf("KT"); KQ = getf("KQ"); eta = getf("eta")
    rR = Float64[]; dT = Float64[]; dQ = Float64[]
    for ln in eachline(fn)
        startswith(ln, "#") && continue
        startswith(ln, "rR") && continue
        t = split(ln, ",")
        length(t) < 5 && continue
        push!(rR, parse(Float64, t[1]))
        push!(dT, parse(Float64, t[4]))
        push!(dQ, parse(Float64, t[5]))
    end
    (; KT, KQ, eta, rR, dT, dQ)
end

ld = read_loading(J)
# total thrust/torque the in-grid disk must carry, in cell units, from
# the VLM coefficients at this (n, D)
const Tt = ld.KT * 1.0 * nRPS^2 * Dcell^4
const Qt = ld.KQ * 1.0 * nRPS^2 * Dcell^5

@printf "=== Ladder 2: VLM → in-grid GradedDisk → KT/KQ/η round-trip ===\n"
@printf "  DTMB 4382  J=%.3f   VLM: KT=%.4f  10KQ=%.4f  η=%.4f\n" J ld.KT 10*ld.KQ ld.eta
@printf "  Grid %d×%d×%d  D=%.0f (R=%.0f) w=%.0f  Re_D=%.0f  steps=%d\n" NX NY NZ Dcell Rcells Wd ReD NSTEPS
@printf "  cell units: U∞=%.3f n=%.5f → prescribed T=%.4f Q=%.4f\n\n" U∞ nRPS Tt Qt
flush(stdout)

disk = GradedDisk(
    center = SVector(cx, cy, cz), axis = SVector(1.0, 0.0, 0.0),
    R = Rcells, w = Wd, thrust = Tt, torque = Qt,
    rR = ld.rR, w_thrust = ld.dT, w_torque = ld.dQ,
)

flow = WaterLily.Flow((NX,NY,NZ), (Float32(U∞),0f0,0f0);
    T=Float32, ν=Float32(U∞*Dcell/ReD), perdir=(2,3), exitBC=true, Δt=0.25)
pois = WaterLily.MultiLevelPoisson(flow.p, flow.μ₀, flow.σ; perdir=(2,3))
udf = (fl,t;kw...) -> disk(fl,t;kw...)

# Axial momentum flux  ∫ ρ u_x² dA  over a y-z plane at axial index i.
function axial_mom_flux(flow, i)
    u = flow.u; s = 0.0
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        ux = Float64(u[i,j,k,1]); s += ux*ux
    end
    s   # ρ=1, cell area =1
end
# Net axial pressure force on a plane. In the confined periodic duct the
# slipstream cannot expand, so a large part of the thrust shows up as a
# streamwise pressure rise (ducted-fan limit) — this term is NOT small and
# is part of the control-volume momentum balance, not a correction.
function plane_pressure(flow, i)
    p = flow.p; s = 0.0
    for k in 2:size(p,3)-1, j in 2:size(p,2)-1
        s += Float64(p[i,j,k])
    end
    s
end
# Swirl angular-momentum flux  ∫ ρ u_x · u_θ · r dA  about the axis.
function swirl_flux(flow, i)
    u = flow.u; s = 0.0
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        y = Float64(j)-1.5-cy; z = Float64(k)-1.5-cz
        r = sqrt(y*y+z*z); r < 1e-6 && continue
        ux = Float64(u[i,j,k,1])
        uy = 0.5*(Float64(u[i,j,k,2]) + Float64(u[i,j+1,k,2]))
        uz = 0.5*(Float64(u[i,j,k,3]) + Float64(u[i,j,k+1,3]))
        uθ = (-z*uy + y*uz)/r
        s += ux*uθ*r
    end
    s
end

@info "running…"; flush(stdout)
for step in 1:NSTEPS
    WaterLily.mom_step!(flow, pois; udf=udf)
    if step % 100 == 0
        @printf "  step=%4d |u|max=%.3f\n" step maximum(abs, flow.u)
        flush(stdout)
    end
end

# Momentum-flux balance. Inlet plane just upstream of the disk, outlet
# plane several R downstream (slipstream developed, before the exit).
const i_in  = max(2, Int(round(cx - 3Rcells)))
const i_out = min(NX-1, Int(round(cx + 6Rcells)))
Mout = axial_mom_flux(flow, i_out)
Min  = axial_mom_flux(flow, i_in)
Pin  = plane_pressure(flow, i_in); Pout = plane_pressure(flow, i_out)
# Control-volume axial-momentum theorem (steady, periodic sides ⇒ no net
# side pressure/flux): T_body = (M_out − M_in) + (P_out − P_in). In the
# confined periodic duct the slipstream cannot contract, so the thrust is
# carried mostly by the streamwise pressure rise (ducted-fan limit) rather
# than a momentum-flux gain — both terms are reported.
T_res = (Mout - Min) + (Pout - Pin)
# resolved torque = swirl angular-momentum flux gained across the disk
Q_res = swirl_flux(flow, i_out) - swirl_flux(flow, i_in)

# round-trip coefficients from the RESOLVED fluxes
KT_res = T_res / (1.0 * nRPS^2 * Dcell^4)
KQ_res = Q_res / (1.0 * nRPS^2 * Dcell^5)
eta_res = KQ_res != 0 ? J*KT_res/(2π*KQ_res) : NaN

# also the prescribed coefficients (exact bookkeeping closure)
KT_pre = Tt / (1.0*nRPS^2*Dcell^4)
KQ_pre = Qt / (1.0*nRPS^2*Dcell^5)

errT = 100*(KT_res - ld.KT)/ld.KT
errQ = 100*(KQ_res - ld.KQ)/ld.KQ
errE = 100*(eta_res - ld.eta)/ld.eta

println()
@printf "  --- prescribed (cell-unit bookkeeping closure) ---\n"
@printf "  KT_prescribed = %.4f  (VLM %.4f, δ=%+.2f%%)\n" KT_pre ld.KT 100*(KT_pre-ld.KT)/ld.KT
@printf "  KQ_prescribed = %.5f (VLM %.5f, δ=%+.2f%%)\n" KQ_pre ld.KQ 100*(KQ_pre-ld.KQ)/ld.KQ
@printf "  --- resolved (in-grid momentum-flux balance) ---\n"
@printf "  T_resolved=%.3f  prescribed=%.3f  (Δmom %.3f + Δp %.3f)\n" T_res Tt (Mout-Min) (Pout-Pin)
@printf "  KT_resolved = %.4f  vs VLM %.4f → δ = %+6.2f %%\n" KT_res ld.KT errT
@printf "  KQ_resolved = %.5f vs VLM %.5f → δ = %+6.2f %%\n" KQ_res ld.KQ errQ
@printf "  η_resolved  = %.4f  vs VLM %.4f → δ = %+6.2f %%\n" eta_res ld.eta errE

# Radial profiles just behind the disk (axial-velocity excess and swirl)
# — the in-grid signature of the imposed radial loading. Compare the
# SHAPE to the VLM dT/dr, dQ/dr that drove the disk.
function radial_profiles(flow, i; nb=24)
    u = flow.u
    edges = range(0, Rcells*1.05; length=nb+1)
    sax = zeros(nb); sth = zeros(nb); cnt = zeros(Int, nb)
    for k in 2:size(u,3)-1, j in 2:size(u,2)-1
        y = Float64(j)-1.5-cy; z = Float64(k)-1.5-cz
        r = sqrt(y*y+z*z)
        b = searchsortedlast(edges, r); (1≤b≤nb) || continue
        ux = Float64(u[i,j,k,1])
        uy = 0.5*(Float64(u[i,j,k,2])+Float64(u[i,j+1,k,2]))
        uz = 0.5*(Float64(u[i,j,k,3])+Float64(u[i,j,k+1,3]))
        uθ = r<1e-6 ? 0.0 : (-z*uy+y*uz)/r
        sax[b]+=ux-U∞; sth[b]+=uθ; cnt[b]+=1
    end
    rc = [0.5*(edges[i]+edges[i+1])/Rcells for i in 1:nb]   # r/R
    ax = [cnt[i]>0 ? sax[i]/cnt[i] : NaN for i in 1:nb]
    th = [cnt[i]>0 ? sth[i]/cnt[i] : NaN for i in 1:nb]
    rc, ax, th
end
i_prof = min(NX-1, Int(round(cx + 1.5Rcells)))
rprof, axexc, swprof = radial_profiles(flow, i_prof)

OUT = joinpath(@__DIR__, "..", "runs", "vlm_disk_roundtrip"); mkpath(OUT)
prof_csv = joinpath(OUT, @sprintf("profile_J%.3f.csv", J))
open(prof_csv, "w") do io
    println(io, "rR,axial_excess,swirl")
    for i in eachindex(rprof)
        @printf io "%.4f,%.5e,%.5e\n" rprof[i] axexc[i] swprof[i]
    end
end
# correlation of the in-grid radial shape with the VLM loading shape,
# over the bins that fall inside the loaded annulus
function shape_corr(rR_grid, vals, rR_vlm, vlm)
    interp(x) = begin
        x ≤ rR_vlm[1] && return vlm[1]; x ≥ rR_vlm[end] && return vlm[end]
        k = findlast(<=(x), rR_vlm); t=(x-rR_vlm[k])/(rR_vlm[k+1]-rR_vlm[k]); vlm[k]+t*(vlm[k+1]-vlm[k])
    end
    vi = Float64[]; vv = Float64[]
    for i in eachindex(rR_grid)
        r = rR_grid[i]
        (r ≥ rR_vlm[1] && r ≤ rR_vlm[end] && !isnan(vals[i])) || continue
        push!(vi, vals[i]); push!(vv, interp(r))
    end
    length(vi) < 3 ? NaN : cor(vi, vv)
end
cT = shape_corr(rprof, axexc, ld.rR, ld.dT)
cQ = shape_corr(rprof, swprof, ld.rR, ld.dQ)
@printf "  radial-shape correlation (in-grid vs VLM): axial r=%.3f  swirl r=%.3f\n" cT cQ

csv = joinpath(OUT, "roundtrip.csv")
newfile = !isfile(csv)
open(csv, newfile ? "w" : "a") do io
    newfile && println(io, "J,Rcells,NX,KT_vlm,KQ_vlm,eta_vlm,KT_res,KQ_res,eta_res,errKT_pct,errKQ_pct,errEta_pct")
    @printf io "%.3f,%.0f,%d,%.5f,%.6f,%.5f,%.5f,%.6f,%.5f,%.2f,%.2f,%.2f\n" J Rcells NX ld.KT ld.KQ ld.eta KT_res KQ_res eta_res errT errQ errE
end
@printf "\nwrote %s\n" csv

pass = abs(errT) ≤ 10 && abs(errQ) ≤ 10
@printf "\n  ±10%% gate (KT & KQ): %s\n" (pass ? "PASS" : "see deltas")
