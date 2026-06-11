#!/usr/bin/env julia
#
# Gentle sloshing validation for MULES + interface compression (c_α).
# Follow-up 2″ of RESULTS-damBreak.md: the damBreak wall-slam showed a
# sharp algebraic interface at ρ-ratio 1000 destabilizes the
# u-advecting coupling in the violent regime; this case probes the
# regime compression was built for — a smooth free surface that never
# breaks — where the payoff is sharpness retention over many periods.
#
# Setup: closed 2D tank (no-slip walls), water depth h = Ly/2, free
# surface tilted by the first standing mode, a·cos(πx/Lx). Gravity
# makes it slosh; linear theory gives ω² = g·k·tanh(k·h), k = π/Lx.
# For the default box (0.585 m, h = 0.2925 m): T ≈ 0.90 s, so 5 s ≈
# 5.5 periods.
#
# Metrics per sample:
#   η_left, η_right — mass-equivalent surface elevation at the walls
#                     (Σ_j α per column, in metres) → period check
#   m/m₀            — global mass conservation
#   width           — interface cells, count(0.05 < α < 0.95)
#                     (sharp ≈ 1–3 per column; diffusion grows it)
#
# ENV knobs:
#   WL_N       grid per side (default 128)
#   WL_TEND    physical end time [s] (default 5.0)
#   WL_AMP     initial tilt amplitude as fraction of Ly (default 0.05)
#   WL_MULES   use step_vof_mules! (default true; false = clamp+repair)
#   WL_CALPHA  compression strength (default 1)
#   WL_TAG     output filename suffix
#
# Output: runs/sloshing_mules/slosh_<TAG>.csv (t, eta_l, eta_r, m_rel, width)

using WaterLily
using VoF
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "sloshing_mules"))
mkpath(OUTDIR)

# --- Physical parameters ----------------------------------------------------
const L_BOX = 0.585                  # tank side [m]
const G_p   = 9.81
const ν_w_p = 1.0e-6
const ν_a_p = 1.8e-5
const AMP   = parse(Float64, get(ENV, "WL_AMP", "0.05"))   # fraction of Ly

const N_RESOL = parse(Int, get(ENV, "WL_N", "128"))
const NX = N_RESOL; const NY = N_RESOL
const ΔX = L_BOX / NX
const H_w   = L_BOX / 2                               # still-water depth [m]
const H_c   = H_w / ΔX                                # in cells

# --- Non-dim scales (damBreak conventions) ----------------------------------
const U_ref = sqrt(G_p * H_w)
const G_c   = G_p * ΔX / U_ref^2
const ν_w_c = ν_w_p / (U_ref * ΔX)
const ν_a_c = ν_a_p / (U_ref * ΔX)

const ρ_w = 1000.0; const ρ_a = 1.0
const μ_w_c = ν_w_c * ρ_w
const μ_a_c = ν_a_c * ρ_a

const T_END        = parse(Float64, get(ENV, "WL_TEND", "5.0"))
const SAMPLE_EVERY = parse(Float64, get(ENV, "WL_SAMPLE", "0.02"))
const POIS_TOL     = parse(Float64, get(ENV, "WL_POIS_TOL", "1e-8"))
const POIS_ITMX    = parse(Int, get(ENV, "WL_POIS_ITMX", "200"))
const TAG          = get(ENV, "WL_TAG", "")
const T_NUM        = Float32

# Linear first-mode prediction (for the log header)
const k1 = π / L_BOX
const ω1 = sqrt(G_p * k1 * tanh(k1 * H_w))
@info @sprintf("""Sloshing setup:
  N=%d  ΔX=%.4f m  h=%.4f m  amp=%.3f·Ly
  linear mode-1: ω=%.3f rad/s  T=%.3f s  (%.1f periods in %.1f s)""",
    N_RESOL, ΔX, H_w, AMP, ω1, 2π/ω1, T_END*ω1/2π, T_END)

# --- VoF + Flow + Poisson ----------------------------------------------------
# Free surface at y = H + a·cos(πx/Lx) in cell coordinates.
const A_c = AMP * NY
α₀ = (i, x) -> begin
    η = H_c + A_c * cospi(x[1] / NX)
    x[2] < η ? 1f0 : 0f0
end

vof = VoFFlow((NX, NY);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
L0 = copy(vof.L)

flow = WaterLily.Flow((NX, NY), (T_NUM(0), T_NUM(0));
    T = T_NUM, ν = VoF.viscosity(vof),
    g = (i, x, t) -> i == 2 ? T_NUM(-G_c) : T_NUM(0), Δt = 0.25)
pois = WaterLily.MultiLevelPoisson(flow.p, L0, flow.σ;
    tol = T_NUM(POIS_TOL), itmx = POIS_ITMX)
sim = (flow = flow, pois = pois)

# --- Diagnostics --------------------------------------------------------------
# Mass-equivalent elevation of column i, in metres above the tank floor.
column_eta(α, i) = sum(@view α[i, 2:NY+1]) * ΔX
water_mass(α) = sum(@view α[2:NX+1, 2:NY+1])
intf_width(α) = count(c -> 0.05f0 < c < 0.95f0, @view α[2:NX+1, 2:NY+1])

# --- Time integration ---------------------------------------------------------
const USE_MULES = parse(Bool, get(ENV, "WL_MULES", "true"))
const CALPHA    = parse(Float64, get(ENV, "WL_CALPHA", "1"))
const Δt_grav_cap = T_NUM(0.3 * sqrt(2 / G_c))

csv = joinpath(OUTDIR, isempty(TAG) ? "slosh.csv" : "slosh_$(TAG).csv")
io = open(csv, "w")
println(io, "t,eta_l,eta_r,m_rel,width")

t_phys = 0.0; t_next = 0.0; nstep = 0
m0 = water_mass(vof.α)
while t_phys < T_END
    sim.flow.Δt[end] = min(sim.flow.Δt[end], Δt_grav_cap)
    WaterLily.mom_step!(sim.flow, sim.pois)
    dt_cell = sim.flow.Δt[end-1]
    if USE_MULES
        step_vof_mules!(vof, sim; dt = dt_cell, c_α = CALPHA)
    else
        step_vof!(vof, sim; dt = dt_cell, mass_repair = true)
    end
    global t_phys += dt_cell * ΔX / U_ref
    global nstep += 1
    if t_phys >= t_next
        @printf io "%.4f,%.5f,%.5f,%.6f,%d\n" t_phys column_eta(vof.α, 2) column_eta(vof.α, NX + 1) water_mass(vof.α)/m0 intf_width(vof.α)
        flush(io)
        global t_next += SAMPLE_EVERY
    end
    if nstep % 200 == 0
        @info @sprintf("step=%5d  t=%.3fs  η_l=%.4fm  m/m0=%.5f  width=%d  α∈[%.3f,%.3f]  |u|=%.2f  niter=%d",
            nstep, t_phys, column_eta(vof.α, 2), water_mass(vof.α)/m0,
            intf_width(vof.α), minimum(vof.α), maximum(vof.α),
            maximum(abs, sim.flow.u), isempty(sim.pois.n) ? -1 : sim.pois.n[end])
    end
    if !isfinite(maximum(vof.α)) || maximum(abs, sim.flow.u) > 100
        @warn "Blew up at step $nstep, t=$t_phys"
        break
    end
end
close(io)
@info "wrote $csv  ($(nstep) steps)"
