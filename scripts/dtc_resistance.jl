#!/usr/bin/env julia
#
# DTC bare-hull calm-water resistance at Fr = 0.218 (Phase-3 step 1).
# ============================================================================
# First resistance computation on the *actual* DTC benchmark geometry
# (every prior resistance number was Wigley or the Containership stand-in).
# Full WaterLily five-package stack: WaterLily (flow + BDIM) + VoF (free
# surface, variable density) + Turbulence (WALE/SST eddy viscosity) +
# ShipShapes (tabulated DTC SDF).
#
# Method (see PLAN-dtc-resistance.md §4):
#   - Hull FIXED at design draft, zero trim (matches the OpenFOAM DTCHull
#     tutorial; free sinkage/trim is a stretch goal).
#   - Fr = 0.218 matched exactly via g_cell = U∞²/(Fr²·L_c), U∞ = 1.
#   - Re = 1e5 (NOT model-scale 6.5e6 — unreachable on a Cartesian cut with
#     no y+ wall function). Compared to reference via the Froude
#     decomposition: C_R = C_T − C_F,ITTC(Re). The ±15 % gate applies to
#     C_R; raw C_T at both Re is reported for transparency.
#   - Grid ladder over WL_NL (cells per Lpp). Run coarsest-first.
#
# Unit system (cell-units, per HANDOFF.md §3 / damBreak recipe):
#   length  = ΔX = Lpp / N_L   (one cell)
#   velocity = U∞ = 1
#   ⇒ ν_cell = U∞·L_c/Re   with L_c = N_L cells
#     g_cell = U∞²/(Fr²·L_c)
#   ρ is a dimensionless ratio (ρ_w = 100, ρ_a = 1 — see note below).
#
# ENV knobs:
#   WL_NL       cells per Lpp                          (default 96)
#   WL_FR       Froude number U/√(g·Lpp)               (default 0.218)
#   WL_RE       Reynolds U·Lpp/ν_water                 (default 1e5)
#   WL_TEND_LU  end time in convective L/U units       (default 5.0)
#   WL_TURB     wale | sst | none                      (default wale)
#   WL_VOF      mules | clamp                          (default clamp)
#   WL_RHO      ρ_w/ρ_a density ratio                  (default 100)
#   WL_RAMP_LU  gravity-ramp length in L/U units       (default 1.0)
#   WL_TAG      output filename suffix                 (default = NL value)
#   WL_HYDRO    if "1", run hydrostatic check only and exit
#
# Output: runs/dtc_resistance/forces_<TAG>.csv
#   columns: t_LU, Fx_pressure, Fx_viscous, CT, CP, CF_meas, umax, wave, niter
# plus a final SUMMARY line printed to stdout (parsed by the analysis step).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Turbulence.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
]; io=devnull)

using WaterLily, VoF, Turbulence, ShipShapes
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics

# --- DTC reference data (model scale 1:59.407) ------------------------------
const REFDATA = get(ENV, "WL_REFDATA",
    abspath(joinpath(@__DIR__, "..", "..", "cerulean-reference-data")))
const SDF_PATH = joinpath(REFDATA, "hulls", "dtc_sdf_h12.5mm.bin")
const LPP_M = 5.976          # Lpp [m], model scale
const B_M   = 0.859          # beam [m]
const T_M   = 0.244          # design draft [m]
const VOL_M = 0.827          # published model displacement [m³]
const S_M   = 6.243          # published wetted surface [m²] (verify in RESULTS)
const G_PHYS = 9.81

# --- knobs ------------------------------------------------------------------
const N_L      = parse(Int,     get(ENV, "WL_NL",      "96"))
const FR       = parse(Float64, get(ENV, "WL_FR",      "0.218"))
const RE       = parse(Float64, get(ENV, "WL_RE",      "1e5"))
const TEND_LU  = parse(Float64, get(ENV, "WL_TEND_LU", "5.0"))
const TURB     = lowercase(get(ENV, "WL_TURB", "wale"))
const VOFMODE  = lowercase(get(ENV, "WL_VOF",  "clamp"))
const RHO      = parse(Float64, get(ENV, "WL_RHO",     "100"))
const RAMP_LU  = parse(Float64, get(ENV, "WL_RAMP_LU", "1.0"))
const HYDRO    = get(ENV, "WL_HYDRO", "0") == "1"
const TAG      = get(ENV, "WL_TAG", string(N_L))

const T_NUM = Float64
const U∞ = 1.0

# --- Geometry: cell-units ---------------------------------------------------
const ΔX = LPP_M / N_L                        # metres per cell
const L_c = Float64(N_L)                      # Lpp in cells (reference length)

# Domain in Lpp (PLAN §4.1): x ∈ [-0.75, +1.75] about midship,
# width ≥ 1.2 Lpp, depth ≥ 0.75 below waterline + ≥ 0.25 above.
# The SDF frame has x≈0 at the bbox centre ≈ midship, the hull spanning
# ±0.589 Lpp incl. overhangs; place midship at x = 0.75 Lpp from inlet.
const DOM_X_LO = -0.75; const DOM_X_HI = 1.85   # Lpp upstream/downstream of midship
const DOM_Y    = 1.30                            # total width in Lpp
const DOM_Z_DN = 0.80                            # below design waterline, Lpp
const DOM_Z_UP = 0.30                            # above design waterline, Lpp

# MultiLevelPoisson requires each grid dimension = a·2ⁿ (n>2), i.e. a
# multiple of 8. Round each to a multiple of 16 so ≥4 multigrid levels exist,
# then derive the *actual* domain extent (in Lpp) and hull anchor cells from
# the rounded dims — the requested fractions above are nominal.
snap16(x) = max(16, 16 * round(Int, x / 16))
const NX = snap16((DOM_X_HI - DOM_X_LO) * N_L)
const NY = snap16(DOM_Y * N_L)
const NZ = snap16((DOM_Z_DN + DOM_Z_UP) * N_L)

# Cell index of midship (SDF x=0), symmetry plane (SDF y=0), waterline (SDF z=0).
# Keep the upstream fetch (0.75 Lpp) and below-water depth (0.80 Lpp) exact;
# the extra rounded cells extend the downstream / above-water regions.
const MIDSHIP_CX = -DOM_X_LO * N_L            # cells from inlet to midship
const CENTER_CY  = NY / 2                     # hull symmetry plane at box centre
const WL_CZ      = DOM_Z_DN * N_L             # design waterline height in cells

# --- Physical scales (cell-units) -------------------------------------------
const G_c   = U∞^2 / (FR^2 * L_c)             # gravity, cell-units
const ν_w_c = U∞ * L_c / RE                   # water kinematic ν, cell-units
const ν_a_c = ν_w_c * 18                      # air ~18× more viscous kinematically
const ρ_w   = RHO
const ρ_a   = 1.0
const μ_w_c = ρ_w * ν_w_c                     # VoFFlow consumes μ = ρ·ν
const μ_a_c = ρ_a * ν_a_c

# Reference wetted surface in cell-units: S_c = S_m / ΔX²
const S_c = S_M / ΔX^2

# Blockage: hull midship section (B·T) vs domain cross-section.
const blockage = (B_M * T_M) / ((DOM_Y * LPP_M) * ((DOM_Z_DN + DOM_Z_UP) * LPP_M))

@printf "=== DTC bare-hull resistance — Fr=%.3f, Re=%.1e ===\n" FR RE
@printf "  Grid          = %d × %d × %d   (N_L=%d cells/Lpp, ΔX=%.4f m)\n" NX NY NZ N_L ΔX
@printf "  L_c           = %.1f cells   U∞=%.1f\n" L_c U∞
@printf "  Domain (Lpp)  = x[%.2f,%.2f] × y%.2f × z[-%.2f,+%.2f]\n" DOM_X_LO DOM_X_HI DOM_Y DOM_Z_DN DOM_Z_UP
@printf "  Midship cx    = %.1f   waterline cz = %.1f   sym cy = %.1f\n" MIDSHIP_CX WL_CZ CENTER_CY
@printf "  g_cell        = %.4e   ν_w=%.3e  ν_a=%.3e\n" G_c ν_w_c ν_a_c
@printf "  ρ_w/ρ_a       = %.0f   S_ref = %.4f m² = %.1f cells²\n" ρ_w/ρ_a S_M S_c
@printf "  blockage      = %.2f %%   turb=%s  vof=%s\n" 100*blockage TURB VOFMODE
@printf "  t_end         = %.1f L/U   ramp=%.1f L/U\n" TEND_LU RAMP_LU
flush(stdout)

# --- Hull SDF ---------------------------------------------------------------
# Map: WaterLily world cell coords → DTC SDF frame (metres).
#   x_sdf_m = (x_cell - MIDSHIP_CX) * ΔX        (SDF x=0 at midship)
#   y_sdf_m = (x_cell[2] - CENTER_CY) * ΔX      (SDF y=0 at symmetry plane)
#   z_sdf_m = (x_cell[3] - WL_CZ) * ΔX          (SDF z=0 at design waterline)
const table = ShipShapes.load_tabulated(SDF_PATH)
hull_map = (x, t) -> SVector(
    (x[1] - MIDSHIP_CX) * ΔX,
    (x[2] - CENTER_CY)  * ΔX,
    (x[3] - WL_CZ)      * ΔX,
)
hull = ShipShapes.tabulated_sdf(table; map = hull_map)

# --- Initial α: water below the design waterline ----------------------------
α₀(_i, x_cell) = (x_cell[3] ≤ WL_CZ) ? 1.0 : 0.0

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "dtc_resistance"))
mkpath(OUTDIR)

# --- Build VoF + (optional) turbulence + Poisson ----------------------------
vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)

# Turbulence model. WALE composes cleanly with VoF: update_νt!(model,u,vof.ν)
# writes ν_mol(per-cell, from VoF) + ν_eddy into model.ν, which we feed to Flow.
# SST stores a *scalar* ν_mol (water value) and writes ν_mol+νt — it ignores
# the air/water ν contrast (acceptable: the hull sits in water, νt dominates).
turb_model = if TURB == "wale"
    WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0), T = T_NUM)
elseif TURB == "smagorinsky"
    Smagorinsky((NX, NY, NZ); Cs = T_NUM(0.17), ν₀ = T_NUM(0), T = T_NUM)
elseif TURB == "sst"
    KOmegaSST((NX, NY, NZ), hull; ν = ν_w_c, k∞ = 1e-4, ω∞ = 1.0,
              perdir = (2,), T = T_NUM)
elseif TURB == "none"
    nothing
else
    error("unknown WL_TURB=$TURB (wale|smagorinsky|sst|none)")
end

is_sst = turb_model isa KOmegaSST
is_les = (turb_model isa WALE) || (turb_model isa Smagorinsky)

# ν closure passed to Flow. The foam-integration WaterLily reads ν through a
# callable `ν(I)`, NOT a raw array (a raw array is rejected). VoF.viscosity
# wraps vof.ν this way; we wrap the turbulence model's ν field identically so
# every in-place refresh is seen by the next mom_step!.
#   none → vof.ν (molecular, refreshed each step by VoF)
#   les  → turb_model.ν (molecular-from-VoF + eddy, refreshed each step)
#   sst  → turb_model.ν (scalar molecular + eddy)
ν_for_flow = if isnothing(turb_model)
    VoF.viscosity(vof)
else
    let νarr = turb_model.ν
        I -> @inbounds νarr[I]
    end
end

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,),
        tol = T_NUM(1e-6), itmx = 50)
end

# Gravity ramp: only the gravity *force* ramps over the first RAMP_LU L/U,
# so the hydrostatic field is established gently (template pattern).
const RAMP_STEPS_TARGET = RAMP_LU * L_c / U∞    # in cell-time units (approx)
g_ramp = Ref(0.0)                                # current ramp scale ∈ [0,1]
g_fn = (i, x, t) -> i == 3 ? T_NUM(-G_c * g_ramp[]) : T_NUM(0)

sim = WaterLily.Simulation((NX, NY, NZ), (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
    T = T_NUM, ν = ν_for_flow, g = g_fn, Δt = 0.25, body = hull, ϵ = 1,
    perdir = (2,), exitBC = true, pois_ctor = vof_pois_ctor, U = T_NUM(U∞))

# Prime the LES eddy-viscosity field (initialised to ν₀=0) so the first
# mom_step! sees the molecular viscosity rather than zero.
is_les && update_νt!(turb_model, sim.flow.u, vof.ν)

# --- Diagnostics ------------------------------------------------------------
function max_wave(α)
    nx, ny, nz = size(α)
    max_eta = 0.0
    @inbounds for j in 2:ny-1, i in 2:nx-1
        for k in nz-1:-1:2
            if α[i, j, k] > 0.5
                eta = (k - 1.5) - WL_CZ
                abs(eta) > abs(max_eta) && (max_eta = eta)
                break
            end
        end
    end
    return max_eta
end

water_mass(α) = (s = 0.0; @inbounds for I in CartesianIndices(α); s += α[I]; end; s)

# --- Hydrostatic sanity check (WL_HYDRO=1): U=0, gravity on, ~50 steps. ------
# Vertical pressure force should ≈ ρ_w·g·V_displaced within ~15 % (BDIM bias
# ~12 %, known) before any production run.
if HYDRO
    @info "Hydrostatic check: quiescent box (zero inflow), gravity ramped on, 80 steps"
    # Dedicated zero-inflow sim so no inlet flux contaminates the static field.
    vofh = VoFFlow((NX, NY, NZ); α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a,
        μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
    gh = Ref(0.0)
    gh_fn = (i, x, t) -> i == 3 ? T_NUM(-G_c * gh[]) : T_NUM(0)
    poisctor_h = function (flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L); L[I] = flow.μ₀[I] * vofh.L[I]; end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,),
            tol = T_NUM(1e-6), itmx = 80)
    end
    simh = WaterLily.Simulation((NX, NY, NZ), (T_NUM(0), T_NUM(0), T_NUM(0)), L_c;
        T = T_NUM, ν = VoF.viscosity(vofh), g = gh_fn, Δt = 0.25, body = hull,
        ϵ = 1, perdir = (2,), exitBC = false, pois_ctor = poisctor_h, U = T_NUM(U∞))
    buoy = ρ_w * G_c * VOL_M / ΔX^3   # ρ·g·V_displaced in cell-units
    for step in 1:80
        gh[] = min(1.0, step / 40)     # ramp gravity over the first 40 steps
        WaterLily.mom_step!(simh.flow, simh.pois)
        dt = simh.flow.Δt[end-1]
        step_vof!(vofh, simh; dt = dt, mass_repair = true)
        if step % 10 == 0
            Fp = WaterLily.pressure_force(simh)
            Fz = Fp[3]                  # vertical pressure force, body→fluid
            @printf "  step=%2d  Fz_press=%.3e  ρgV=%.3e  ratio=%.3f  |u|=%.2e  g=%.2f  nit=%d\n" (
                step) Fz buoy (-Fz/buoy) maximum(abs, simh.flow.u) gh[] simh.pois.n[end]
            flush(stdout)
        end
    end
    @info "Hydrostatic check done — final ratio should be ≈1 (BDIM bias ~12% expected)."
    exit(0)
end

# --- Production time integration --------------------------------------------
const csv_path = joinpath(OUTDIR, "forces_$(TAG).csv")
io = open(csv_path, "w")
println(io, "t_LU,Fx_pressure,Fx_viscous,CT,CP,CF_meas,umax,wave,niter")

# Convective end-time in cell-time units: t_end_cell = TEND_LU · L_c / U∞.
const T_END_CELL = TEND_LU * L_c / U∞
const RAMP_END_CELL = RAMP_LU * L_c / U∞

t_cell = 0.0
step = 0
m0 = water_mass(vof.α)
blew_up = false

@info @sprintf("Stepping to t=%.1f L/U (= %.0f cell-time units)…", TEND_LU, T_END_CELL)
flush(stdout)

while t_cell < T_END_CELL
    global step += 1
    g_ramp[] = min(1.0, t_cell / max(RAMP_END_CELL, eps()))

    WaterLily.mom_step!(sim.flow, sim.pois)
    dt = sim.flow.Δt[end-1]

    # VoF advection (PLAN §2: plain MULES or clamp/mass_repair; never c_α=1).
    if VOFMODE == "clamp"
        step_vof!(vof, sim; dt = dt, mass_repair = true)
    else
        step_vof_mules!(vof, sim; dt = dt, perdir = (2,))
    end

    # Refresh eddy viscosity.
    if is_les
        update_νt!(turb_model, sim.flow.u, vof.ν)
    elseif is_sst
        step_sst!(turb_model, sim.flow.u, dt; production = :kato_launder)
    end

    global t_cell += dt
    t_LU = t_cell / L_c

    if step % 5 == 0 || step ≤ 3
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        Fxp = -Float64(Fp[1])             # drag on hull (negate body→fluid)
        Fxv = -Float64(Fv[1])
        denom = 0.5 * ρ_w * U∞^2 * S_c
        CT = (Fxp + Fxv) / denom
        CP = Fxp / denom
        CF = Fxv / denom
        umax = maximum(abs, sim.flow.u)
        wv = max_wave(vof.α)
        niter = isempty(sim.pois.n) ? -1 : sim.pois.n[end]
        @printf io "%.5f,%.6e,%.6e,%.6e,%.6e,%.6e,%.4f,%.3f,%d\n" (
            t_LU) Fxp Fxv CT CP CF umax wv niter
        flush(io)
        if step % 25 == 0 || step ≤ 3
            mass = water_mass(vof.α) / m0
            @printf "  t=%.3f L/U  CT=%.4e (CP=%.4e CF=%.4e)  |u|=%.2f  η=%+.1f  m/m0=%.4f  nit=%d  g=%.2f\n" (
                t_LU) CT CP CF umax wv mass niter g_ramp[]
            flush(stdout)
        end
        if !isfinite(umax) || umax > 25
            @warn "Blow-up at step $step (t=$t_LU L/U, |u|=$umax)"
            global blew_up = true
            break
        end
    end
end
close(io)

# --- Final averaged summary -------------------------------------------------
# Re-read the CSV, average CT over the settled window t_LU ∈ [2, end].
using DelimitedFiles
data, hdr = readdlm(csv_path, ','; header = true)
t_col   = Float64.(data[:, 1])
CT_col  = Float64.(data[:, 4])
CP_col  = Float64.(data[:, 5])
CF_col  = Float64.(data[:, 6])

const WIN_LO = 2.0
mask = t_col .≥ WIN_LO
n_win = count(mask)

# ITTC-57 friction line for the simulated Re.
CF_ittc(Re) = 0.075 / (log10(Re) - 2)^2
CFsim = CF_ittc(RE)

println("\n" * "="^78)
if n_win ≥ 3
    CT_mean = mean(CT_col[mask]); CT_std = std(CT_col[mask])
    CP_mean = mean(CP_col[mask]); CF_mean = mean(CF_col[mask])
    # Split-half convergence check.
    idx = findall(mask); half = idx[1:end÷2]; half2 = idx[end÷2+1:end]
    CT_h1 = mean(CT_col[half]); CT_h2 = mean(CT_col[half2])
    split_pct = 100 * abs(CT_h1 - CT_h2) / max(abs(CT_mean), eps())
    CR_meas = CT_mean - CFsim
    @printf "SUMMARY NL=%d Fr=%.3f Re=%.1e  NX=%d NY=%d NZ=%d  ΔX=%.4f\n" N_L FR RE NX NY NZ ΔX
    @printf "SUMMARY  CT=%.6e ± %.2e  (window t/[%.1f,%.2f], n=%d)\n" CT_mean CT_std WIN_LO maximum(t_col) n_win
    @printf "SUMMARY  CP=%.6e  CF_meas=%.6e  CF_ITTC(Re_sim)=%.6e\n" CP_mean CF_mean CFsim
    @printf "SUMMARY  C_R = CT - CF_ITTC(Re_sim) = %.6e\n" CR_meas
    @printf "SUMMARY  split-half CT agreement = %.2f %% (target < 2%%)\n" split_pct
    @printf "SUMMARY  blew_up=%s\n" blew_up
else
    @printf "SUMMARY NL=%d  INSUFFICIENT settled window (n=%d < 3). blew_up=%s\n" N_L n_win blew_up
end
println("="^78)
@printf "wrote %s\n" csv_path
