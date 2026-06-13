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

# --- Fix 1: sponge / wave-damping + gentle inflow ramp knobs ----------------
# Rayleigh damping layer near the outlet (and a thinner inlet layer) relaxes
# u → u_ref = (U∞,0,0) and α → α_still through a body force applied via the
# `udf` hook, killing the trapped longitudinal sloshing mode (round-1 §1).
const SPONGE     = get(ENV, "WL_SPONGE", "1") == "1"     # on by default round-2
const SPONGE_W   = parse(Float64, get(ENV, "WL_SPONGE_W",     "0.65"))  # outlet zone width [Lpp]
const SPONGE_WIN = parse(Float64, get(ENV, "WL_SPONGE_WIN",   "0.25"))  # inlet  zone width [Lpp]
const SPONGE_SIG = parse(Float64, get(ENV, "WL_SPONGE_SIGMA", "2.5"))   # σ_max · Lpp/U∞
const SPONGE_ATOP= parse(Float64, get(ENV, "WL_SPONGE_TOP",   "0.0"))   # top zone width [Lpp], 0=off
# The velocity sponge (-σ(u−u_ref) body force) and the α (surface) damping are
# now independently toggleable. The velocity body force enters the pressure
# projection and, if strong, stalls the multigrid (nit→itmx); the α-damping is
# a post-VoF relaxation that does NOT touch the Poisson and is the piece that
# actually flattens the standing wave. Default: α-damp ON, velocity sponge OFF.
const SPONGE_UVEL = get(ENV, "WL_SPONGE_UVEL", "0") == "1"   # velocity body-force sponge
# POST-projection velocity relaxation (Poisson-safe): after mom_step! returns a
# divergence-free field, relax u → u_ref in the sponge by an explicit factor
# (1 − w), w = SPONGE_UPOST·(σ/σ_max). Applied to the CONVERGED field (not
# inside the solver), so it never stalls the multigrid; the small divergence it
# introduces is removed by the next step's projection. This bleeds the standing
# wave's kinetic energy that α-damping (surface only) leaves behind.
const SPONGE_UPOST = parse(Float64, get(ENV, "WL_SPONGE_UPOST", "0.0"))  # 0=off
const SPONGE_ADMP= parse(Float64, get(ENV, "WL_SPONGE_ADAMP", "1.0"))   # α-damp blend factor ∈[0,1]
# Gentle combined ramp: scale the *inflow target* (uBC) over the first
# RAMP_LU L/U as well as gravity, so the freestream is not impulsive.
const INFLOW_RAMP = get(ENV, "WL_INFLOW_RAMP", "1") == "1"

# --- Fix 3: near-wall stress model (BDIM Spalding wall function) knobs -------
# Turns Turbulence.jl's existing wall function on so resolved friction becomes
# physical. SA preferred (its BDIM wall function hit the channel log-law within
# 7.4 %); SST has a native ω-wall and normally should NOT also use wallfn.
const WALLFN  = get(ENV, "WL_WALLFN", "0") == "1"
const WL_BAND = let s = get(ENV, "WL_BAND", "1,3")
    p = split(s, ','); (parse(Float64, p[1]), parse(Float64, p[2]))
end

const T_NUM = Float64
const U∞ = 1.0

# --- Geometry: cell-units ---------------------------------------------------
const ΔX = LPP_M / N_L                        # metres per cell
const L_c = Float64(N_L)                      # Lpp in cells (reference length)

# Domain in Lpp (PLAN §4.1): x ∈ [-0.75, +1.75] about midship,
# width ≥ 1.2 Lpp, depth ≥ 0.75 below waterline + ≥ 0.25 above.
# The SDF frame has x≈0 at the bbox centre ≈ midship, the hull spanning
# ±0.589 Lpp incl. overhangs; place midship at x = 0.75 Lpp from inlet.
# Round 2: the downstream extent is now a knob so the outlet sponge (Fix 1) has
# wake room AFT of the stern (hull bbox spans ±0.589 Lpp; stern ≈ +0.59 Lpp from
# midship = +1.34 Lpp from inlet). Lengthen to keep the sponge clear of the hull.
const DOM_X_LO = -0.75
const DOM_X_HI = parse(Float64, get(ENV, "WL_DOM_XHI", "2.50"))  # downstream extent [Lpp]
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

# --- Initial α --------------------------------------------------------------
# Free surface: water below the design waterline (default).
# Deep-water reference (WL_DEEP=1): α=1 everywhere → fully-submerged hull,
# no surface piercing. Its drag is the form+friction baseline; subtracting it
# from the free-surface run isolates the WAVE-MAKING component (the
# RESULTS-drag-decomposition.md method), which removes the hull's form drag
# and the common BDIM-staircase bias.
const DEEP = get(ENV, "WL_DEEP", "0") == "1"
α₀ = DEEP ? ((_i, x_cell) -> 1.0) : ((_i, x_cell) -> (x_cell[3] ≤ WL_CZ) ? 1.0 : 0.0)

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
elseif TURB == "sa"
    # Spalart–Allmaras: the model whose BDIM wall function hit the channel
    # log-law within 7.4 % (RESULTS-channel-sa.md). ν̃∞=3ν fully turbulent.
    SpalartAllmaras((NX, NY, NZ), hull; ν = ν_w_c, ν̃∞ = 3ν_w_c,
                    perdir = (2,), T = T_NUM)
elseif TURB == "none"
    nothing
else
    error("unknown WL_TURB=$TURB (wale|smagorinsky|sst|sa|none)")
end

is_sst = turb_model isa KOmegaSST
is_sa  = turb_model isa SpalartAllmaras
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
const RAMP_END_CELL_K = RAMP_LU * L_c / U∞       # ramp length in cell-time units
g_ramp = Ref(0.0)                                # current ramp scale ∈ [0,1]
g_fn = (i, x, t) -> i == 3 ? T_NUM(-G_c * g_ramp[]) : T_NUM(0)

# Smooth ramp profile (0→1 over [0,1], C¹ at both ends): smoothstep. Written
# generically so ForwardDiff Duals (used for the reference-frame accel) survive.
@inline smoothramp(s) = s ≤ zero(s) ? zero(s) : s ≥ one(s) ? one(s) : s*s*(3 - 2s)

# Fix 1: gentle inflow ramp. Instead of an impulsive freestream, ramp the
# inlet/reference velocity 0→U∞ over the first RAMP_LU L/U. WaterLily evaluates
# uBC(i,x,t) at the actual flow time t (cell-time units), so we ramp directly
# off t — self-consistent with the convective transient the sponge absorbs.
# NOTE: return type must follow `t` (WaterLily differentiates uBC w.r.t. t via
# ForwardDiff for the reference-frame acceleration term), so do NOT hard-cast to
# Float64 here — let Duals propagate. smoothramp is written generically.
uBC_fn = if INFLOW_RAMP
    (i, x, t) -> i == 1 ? U∞ * smoothramp(t / max(RAMP_END_CELL_K, eps())) : zero(t)
else
    (i, x, t) -> i == 1 ? oftype(t, U∞) : zero(t)
end

sim = WaterLily.Simulation((NX, NY, NZ), uBC_fn, L_c;
    T = T_NUM, ν = ν_for_flow, g = g_fn, Δt = 0.25, body = hull, ϵ = 1,
    perdir = (2,), exitBC = true, pois_ctor = vof_pois_ctor, U = T_NUM(U∞))

# Prime the LES eddy-viscosity field (initialised to ν₀=0) so the first
# mom_step! sees the molecular viscosity rather than zero.
is_les && update_νt!(turb_model, sim.flow.u, vof.ν)

# ============================================================================
# Fix 1 — Rayleigh damping (sponge) layer.
# ============================================================================
# Precompute a per-cell damping rate σ(x) [1/cell-time] that is 0 over the hull
# and near wake, rising smoothly (quadratic) to σ_max at the domain outlet, plus
# a thinner inlet band that absorbs the upstream-running reflected wave, plus an
# optional top band. The udf adds  f[I,i] += -σ(x)·(u[I,i] - u_ref_i)  to the
# body-force accumulator (BDIM's μ₀ mask zeros it inside the hull automatically).
# σ_max is set in physical (Lpp/U∞) units: σ_max_cell = SPONGE_SIG · U∞ / L_c.
# All arrays are ghost-padded to (NX+2,NY+2,NZ+2); cell-centre coordinate of
# index I along a dim is (I − 1.5) (VoF._ic_loc convention, matches WL_CZ and
# the hull_map). Interior cells span index 2..N+1 → coords 0.5..N−0.5.
const σ_max_cell = SPONGE_SIG * U∞ / L_c
const X_OUT      = (NX + 1) - 1.5                 # last interior cell-x = NX-0.5
const X_IN       = 2 - 1.5                         # first interior cell-x = 0.5
const SP_W_CELL   = SPONGE_W   * L_c              # outlet zone width [cells]
const SP_WIN_CELL = SPONGE_WIN * L_c              # inlet  zone width [cells]
const SP_TOP_CELL = SPONGE_ATOP * L_c             # top    zone width [cells]
const Z_TOP       = (NZ + 1) - 1.5                # last interior cell-z

# σ array, ghost-padded to match flow.u / flow.f / vof.α. Built once.
σ_sponge = SPONGE ? zeros(T_NUM, NX + 2, NY + 2, NZ + 2) : nothing
if SPONGE
    inner_out = X_OUT - SP_W_CELL                 # inner edge of outlet zone
    inner_in  = X_IN  + SP_WIN_CELL               # inner edge of inlet  zone
    inner_top = Z_TOP - SP_TOP_CELL
    @inbounds for I in CartesianIndices(σ_sponge)
        i, k = I[1], I[3]
        xc = (i - 1.5); zc = (k - 1.5)
        s = 0.0
        # Outlet band: quadratic ramp 0→1 from inner edge to outlet.
        if SP_W_CELL > 0 && xc > inner_out
            ξ = (xc - inner_out) / SP_W_CELL
            s = max(s, ξ*ξ)
        end
        # Inlet band: quadratic ramp 1→0 from inlet to inner edge.
        if SP_WIN_CELL > 0 && xc < inner_in
            ξ = (inner_in - xc) / SP_WIN_CELL
            s = max(s, ξ*ξ)
        end
        # Top band (above the air gap): quadratic ramp 0→1 toward the lid.
        if SP_TOP_CELL > 0 && zc > inner_top
            ξ = (zc - inner_top) / SP_TOP_CELL
            s = max(s, ξ*ξ)
        end
        σ_sponge[I] = σ_max_cell * min(s, 1.0)
    end
    # Hull bbox in cell-x (from the SDF table extent ±0.589 Lpp about midship).
    bow_cx  = (0.75 - 0.589) * L_c        # forward perpendicular, cells from inlet
    stern_cx = (0.75 + 0.590) * L_c       # aft  perpendicular
    @printf "  SPONGE on: σ_max=%.3e/ct (=%.1f U/L)  out=%.2f Lpp inlet=%.2f Lpp top=%.2f Lpp  α-damp=%.2f  uvel-sponge=%s\n" (
        σ_max_cell) SPONGE_SIG SPONGE_W SPONGE_WIN SPONGE_ATOP SPONGE_ADMP (SPONGE_UVEL ? "ON" : "off")
    @printf "  post-sponge (Poisson-safe u-relax) = %.2f\n" SPONGE_UPOST
    @printf "  sponge edges: inlet→%.1f cells | hull[%.1f,%.1f] | %.1f cells←outlet  (X_OUT=%.1f)\n" (
        inner_in) bow_cx stern_cx inner_out X_OUT
    if inner_in > bow_cx
        @warn "INLET sponge intrudes on hull bow (inner_in=$(round(inner_in,digits=1)) > bow=$(round(bow_cx,digits=1)) cells) — reduce WL_SPONGE_WIN"
    end
    if inner_out < stern_cx
        @warn "OUTLET sponge intrudes on hull stern (inner_out=$(round(inner_out,digits=1)) < stern=$(round(stern_cx,digits=1)) cells) — reduce WL_SPONGE_W or raise WL_DOM_XHI"
    end
    flush(stdout)
end

# ============================================================================
# Fix 3 — gate the BDIM wall function to the water.
# ============================================================================
# step_sa!(...; wallfn=true) writes the Spalding-law ν override into m.ν over
# the WHOLE near-wall band, including the air-side BDIM cells above the
# waterline. Friction must only load in the water, so after the step we restore
# the un-wall-functioned effective viscosity (ν_mol + ν_t,model) in any band
# cell whose VoF fraction α ≤ 0.5. We hold the pre-wall-function ν in a buffer.
const _has_sa = (TURB == "sa")
ν_nowf_buf = _has_sa ? similar(turb_model.ν) : nothing
function gate_wallfn_to_water!(ν, α, νmol)
    # ν currently holds the wall-functioned field; ν_nowf_buf holds the
    # pre-override (model-only) field captured just before apply_wall_function!.
    @inbounds for I in CartesianIndices(ν)
        # α is interior-sized (NX,NY,NZ); ν is (NX,NY,NZ) too here.
        if α[I] ≤ 0.5
            ν[I] = ν_nowf_buf[I]
        end
    end
    return nothing
end

# udf closure: relax velocity toward u_ref=(U∞,0,0) inside the sponge.
# Signature (flow, t; kwargs...) per WaterLily's udf! hook. The force is an
# acceleration (BDIM does u += dt·f), so σ has units 1/cell-time. f[I,i] is the
# i-th velocity-component face value; we damp toward the freestream component.
function sponge_force!(flow, t; kwargs...)
    σ = σ_sponge
    f = flow.f
    u = flow.u
    D = ndims(flow.p)
    @inbounds for i in 1:D
        uref = (i == 1) ? T_NUM(U∞) : zero(T_NUM)
        for I in CartesianIndices(σ)
            σI = σ[I]
            σI == 0 && continue
            f[I, i] -= σI * (u[I, i] - uref)
        end
    end
    return nothing
end

# α-damping: after the VoF step, blend α toward the still-water profile inside
# the sponge so the OUTGOING wave's surface deformation is flattened (this is
# the piece that actually stops the standing wave, not just its velocity).
# α_still = 1 below the design waterline, 0 above. Blend weight = SPONGE_ADMP·
# (σ/σ_max) so it tracks the same spatial ramp as the velocity sponge.
function damp_alpha!(α)
    (SPONGE && SPONGE_ADMP > 0 && σ_max_cell > 0) || return nothing
    σ = σ_sponge; b = SPONGE_ADMP; inv_smax = 1.0 / σ_max_cell
    R = CartesianIndices(α)
    Threads.@threads for I in R
        @inbounds begin
            σI = σ[I]
            if σI != 0
                w = b * (σI * inv_smax)
                # still-water profile: α=1 below the design waterline (same as
                # the α₀ IC: cell-centre z = k − 1.5 vs WL_CZ).
                α_still = ((I[3] - 1.5) ≤ WL_CZ) ? one(T_NUM) : zero(T_NUM)
                α[I] += T_NUM(w) * (α_still - α[I])
            end
        end
    end
    return nothing
end

# Post-projection velocity relaxation (Poisson-safe sponge). Call AFTER
# mom_step! on the converged, divergence-free flow.u. Relaxes every velocity
# component toward u_ref=(U∞,0,0) by factor w in the sponge; the next step's
# projection re-imposes incompressibility, so this never enters the solver.
function post_sponge!(u)
    (SPONGE && SPONGE_UPOST > 0 && σ_max_cell > 0) || return nothing
    σ = σ_sponge; b = SPONGE_UPOST; inv_smax = 1.0 / σ_max_cell
    D = ndims(u) - 1
    R = CartesianIndices(σ)
    for i in 1:D
        uref = (i == 1) ? T_NUM(U∞) : zero(T_NUM)
        Threads.@threads for I in R
            @inbounds begin
                σI = σ[I]
                if σI != 0
                    w = min(b * (σI * inv_smax), 1.0)
                    u[I, i] += T_NUM(w) * (uref - u[I, i])
                end
            end
        end
    end
    return nothing
end

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
    g_ramp[] = smoothramp(t_cell / max(RAMP_END_CELL, eps()))

    # Fix 1: forward the velocity sponge udf only if explicitly enabled (it
    # enters the Poisson and can stall the solver — default OFF; α-damp below
    # is the real slosh killer and is projection-free).
    if SPONGE && SPONGE_UVEL
        WaterLily.mom_step!(sim.flow, sim.pois; udf = sponge_force!)
    else
        WaterLily.mom_step!(sim.flow, sim.pois)
    end
    dt = sim.flow.Δt[end-1]

    # Fix 1 (Poisson-safe): relax velocity toward freestream in the sponge on the
    # converged field, before VoF advects with it. No-op if WL_SPONGE_UPOST=0.
    post_sponge!(sim.flow.u)

    # VoF advection (PLAN §2: plain MULES or clamp/mass_repair; never c_α=1).
    if VOFMODE == "clamp"
        step_vof!(vof, sim; dt = dt, mass_repair = true)
    else
        step_vof_mules!(vof, sim; dt = dt, perdir = (2,))
    end
    # Fix 1: damp α toward still water inside the sponge (flattens the outgoing
    # wave's surface deformation — the piece that actually kills the slosh).
    damp_alpha!(vof.α)

    # Refresh eddy viscosity.
    if is_les
        update_νt!(turb_model, sim.flow.u, vof.ν)
    elseif is_sst
        step_sst!(turb_model, sim.flow.u, dt; production = :kato_launder,
                  wallfn = WALLFN, band = (T_NUM(WL_BAND[1]), T_NUM(WL_BAND[2])))
    elseif is_sa
        # Step SA WITHOUT the internal wall function so we can capture the
        # model-only ν, then apply + gate the wall function to water cells.
        step_sa!(turb_model, sim.flow.u, dt)
        if WALLFN
            copyto!(ν_nowf_buf, turb_model.ν)          # model-only ν snapshot
            apply_wall_function!(turb_model.ν, sim.flow.u, turb_model.d, ν_w_c;
                band = (T_NUM(WL_BAND[1]), T_NUM(WL_BAND[2])), perdir = (2,))
            gate_wallfn_to_water!(turb_model.ν, vof.α, ν_w_c)  # restore air cells
        end
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
