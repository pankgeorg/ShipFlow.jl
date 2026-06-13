#!/usr/bin/env julia
#
# Analysis for the DTC bare-hull resistance grid ladder.
# Reads runs/dtc_resistance/forces_<TAG>.csv for each grid in WL_GRIDS,
# averages C_T over the settled window, applies the Froude decomposition,
# does a Richardson/GCI extrapolation, and compares against the published
# DTC reference (el Moctar, Shigunov & Zorn 2012, Table 4, Fr=0.218).
#
# ENV:
#   WL_GRIDS   comma-separated "NL:TAG" pairs, coarsest first
#              (default "96:96,144:144,192:192")
#   WL_WIN_LO  settled-window start in L/U (default 2.0)
#
# Output: runs/dtc_resistance/ladder_summary.csv  +  printed table.

using Printf, Statistics, DelimitedFiles

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "dtc_resistance"))
const WIN_LO = parse(Float64, get(ENV, "WL_WIN_LO", "2.0"))
const LPP_M  = 5.976

# --- Published reference (el Moctar/Shigunov/Zorn 2012, Tab. 4, model scale) -
# FULL Table 4 (verified against the source PDF, round-2): the SVA Potsdam
# towing tank tested SIX speeds, Fr = 0.174 … 0.218 ONLY. There is NO
# experimental point at Fr = 0.28 or 0.33 (those are above the tested
# envelope), so any sim at those Fr can only be reported with the
# ITTC-friction-subtracted C_R and an explicit "no experimental reference"
# caveat — never invent a number. C_T, C_F are ×1e-3; C_W is ×1e-4.
# Columns: Fr => (Re_M, CT, CF_ITTC, CW)  (CW = CT − (1+k)·CF, k = 0.094)
const K_FORM   = 0.094
const TAB4 = Dict(
    0.174 => (7.319e6, 3.661e-3, 3.170e-3, 1.932e-4),
    0.183 => (7.681e6, 3.605e-3, 3.142e-3, 1.672e-4),
    0.192 => (8.054e6, 3.588e-3, 3.116e-3, 1.791e-4),
    0.200 => (8.415e6, 3.602e-3, 3.092e-3, 2.194e-4),
    0.209 => (8.783e6, 3.623e-3, 3.069e-3, 2.660e-4),
    0.218 => (9.145e6, 3.670e-3, 3.047e-3, 3.360e-4),
)
# Default reference point (highest tested = largest wave fraction).
const FR_REF   = 0.218
const RE_M     = TAB4[FR_REF][1]
const CT_EXP   = TAB4[FR_REF][2]
const CF_EXP   = TAB4[FR_REF][3]      # ITTC-57 at Re_m (matches the paper)
const CW_EXP   = TAB4[FR_REF][4]      # paper, with form factor k=0.094
# Simple Froude residuary (no form factor) — the PLAN's C_R definition:
const CR_REF   = CT_EXP - CF_EXP            # = 0.623e-3

CF_ittc(Re) = 0.075 / (log10(Re) - 2)^2

# Nearest tabulated reference for an arbitrary simulated Fr (within 0.005),
# else nothing. Returns (Fr_ref, CT, CF, CW, CR=CT−CF).
function tab4_ref(fr)
    best = nothing; bd = 0.005
    for (f, v) in TAB4
        d = abs(f - fr)
        d ≤ bd && (bd = d; best = (f, v[2], v[3], v[4], v[2]-v[3]))
    end
    return best
end

# --- Per-grid window average ------------------------------------------------
function grid_stats(tag)
    path = joinpath(OUTDIR, "forces_$(tag).csv")
    isfile(path) || return nothing
    data, _ = readdlm(path, ','; header = true)
    t  = Float64.(data[:, 1]); CT = Float64.(data[:, 4])
    CP = Float64.(data[:, 5]); CF = Float64.(data[:, 6])
    tmax = maximum(t)
    lo = min(WIN_LO, 0.6 * tmax)        # if run is short, use the last 40%
    mask = t .≥ lo
    count(mask) < 3 && return (tag=tag, n=0, tmax=tmax, lo=lo)
    idx = findall(mask)
    h1 = idx[1:end÷2]; h2 = idx[end÷2+1:end]
    split = 100 * abs(mean(CT[h1]) - mean(CT[h2])) / max(abs(mean(CT[mask])), eps())
    return (tag=tag, n=count(mask), tmax=tmax, lo=lo,
            CT=mean(CT[mask]), CT_std=std(CT[mask]),
            CP=mean(CP[mask]), CF=mean(CF[mask]), split=split)
end

# --- Parse grid list --------------------------------------------------------
const GRIDS = let s = get(ENV, "WL_GRIDS", "96:96,144:144,192:192")
    [(parse(Int, split(p, ':')[1]), split(p, ':')[2]) for p in split(s, ',')]
end

rows = []
for (nl, tag) in GRIDS
    st = grid_stats(tag)
    st === nothing && (@printf("  [skip] NL=%d tag=%s — no CSV\n", nl, tag); continue)
    push!(rows, (nl=nl, st...))
end
isempty(rows) && error("no grid CSVs found in $OUTDIR")

println("="^90)
@printf "DTC bare-hull resistance — grid ladder (Fr=%.3f, Re_sim from each run)\n" FR_REF
println("="^90)
@printf "%-6s %-8s %6s %7s   %-12s %-12s %-12s %-10s %-7s\n" "NL" "ΔX[m]" "n_win" "t_end" "CT" "CP(press)" "CF(visc)" "CR=CT-CFsim" "split%"
println("-"^90)

CFsim = CF_ittc(1e5)   # all production runs are at Re=1e5

ladder = []
for r in rows
    haskey(r, :CT) || (@printf("%-6d %-8s   n=%d  (insufficient window, t_end=%.2f, lo=%.2f)\n", r.nl, "", r.n, r.tmax, r.lo); continue)
    ΔX = LPP_M / r.nl
    CR = r.CT - CFsim
    @printf "%-6d %-8.4f %6d %7.2f   %-12.5e %-12.5e %-12.5e %-10.4e %-7.2f\n" (
        r.nl) ΔX r.n r.tmax r.CT r.CP r.CF CR r.split
    push!(ladder, (nl=r.nl, CT=r.CT, CP=r.CP, CF=r.CF, CR=CR))
end
println("-"^90)

# --- Richardson / GCI on CT (and CR) ----------------------------------------
function richardson(vals, nls)
    # vals coarsest→finest, refinement ratio from cell counts.
    length(vals) < 3 && return nothing
    f3, f2, f1 = vals[end-2], vals[end-1], vals[end]   # coarse, med, fine
    n3, n2, n1 = nls[end-2], nls[end-1], nls[end]
    r21 = n1 / n2; r32 = n2 / n3
    ε21 = f1 - f2; ε32 = f2 - f3
    R = ε21 / ε32
    if !(0 < R < 1)
        return (monotone=false, R=R)              # oscillatory or divergent
    end
    p = abs(log(abs(ε32 / ε21)) / log(r21))
    f_ext = f1 + ε21 / (r21^p - 1)                 # Richardson extrapolate
    e_ext = abs((f_ext - f1) / f_ext)
    GCI_fine = 1.25 * abs(ε21 / f1) / (r21^p - 1)  # Roache GCI (Fs=1.25)
    return (monotone=true, p=p, f_ext=f_ext, e_ext=e_ext, GCI=GCI_fine)
end

println()
if length(ladder) ≥ 3
    nls = [l.nl for l in ladder]
    for (name, vals) in (("CT", [l.CT for l in ladder]), ("CR", [l.CR for l in ladder]))
        rc = richardson(vals, nls)
        if rc === nothing
            @printf "  %s: need ≥3 grids for Richardson\n" name
        elseif rc.monotone
            @printf "  %s: monotone, p=%.2f, extrap=%.5e, GCI_fine=%.1f%% (band ±%.1f%%)\n" (
                name) rc.p rc.f_ext 100*rc.GCI 100*rc.GCI
        else
            @printf "  %s: NON-monotone (R=%.2f) — report spread, not GCI\n" name rc.R
        end
    end
else
    @printf "  Only %d grid(s) settled — Richardson needs 3. Reporting trend only.\n" length(ladder)
end

# --- Reference comparison ---------------------------------------------------
println("\n" * "="^90)
println("Reference comparison (el Moctar/Shigunov/Zorn 2012, Table 4, Fr=0.218, model scale)")
println("="^90)
@printf "  Experiment:  CT=%.4e  CF_ITTC(Re_m=%.2e)=%.4e  CR(=CT-CF)=%.4e  CW(form-factor)=%.4e\n" (
    CT_EXP) RE_M CF_EXP CR_REF CW_EXP
@printf "  Sim Re=1e5:  CF_ITTC(Re_sim)=%.4e\n" CFsim
println("-"^90)
if !isempty(ladder)
    fine = ladder[end]
    # The classical Froude subtraction C_T,sim − C_F,ITTC(Re_sim) is INVALID
    # here: BDIM at this resolution/Re develops ~no turbulent wall friction
    # (CF_meas ≈ 1e-4 « ITTC 8.3e-3), so the subtraction overcorrects into a
    # negative number. The simulation's resolved force is essentially ALL
    # pressure drag (form + wave + BDIM-staircase + tank sloshing). The
    # physically meaningful comparison is therefore the PRESSURE coefficient
    # against the reference RESIDUARY coefficient (which is itself the
    # pressure-origin resistance):  C_P,sim  vs  C_R,ref.
    err_CP = 100 * (fine.CP - CR_REF) / CR_REF
    @printf "  [Froude C_R = CT − CF_ITTC(Re_sim) = %.4e — NEGATIVE/invalid, see note]\n" fine.CR
    @printf "  Primary:  C_P,sim=%.4e (pressure/residuary)  vs  C_R,ref=%.4e  →  Δ = %+.0f %% (gate ±15%%)\n" (
        fine.CP) CR_REF err_CP
    verdict = abs(err_CP) ≤ 15 ? "PASS" : "FAIL"
    @printf "  C_R gate (±15%%) on C_P: %s\n" verdict
    @printf "  NOTE: C_T,sim is sloshing-dominated (std≈56%% of mean); the true\n"
    @printf "        wave signal at Fr=0.218 (C_R_ref=0.6e-3) is near the method's noise floor.\n"
end
println("="^90)

# --- Write summary CSV ------------------------------------------------------
open(joinpath(OUTDIR, "ladder_summary.csv"), "w") do io
    println(io, "NL,dX_m,CT,CP_press,CF_visc,CR,CR_ref,err_CR_pct")
    for l in ladder
        @printf io "%d,%.5f,%.6e,%.6e,%.6e,%.6e,%.6e,%.2f\n" (
            l.nl) (LPP_M/l.nl) l.CT l.CP l.CF l.CR CR_REF 100*(l.CR-CR_REF)/CR_REF
    end
end
@printf "wrote %s\n" joinpath(OUTDIR, "ladder_summary.csv")
