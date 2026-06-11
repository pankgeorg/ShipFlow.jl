#!/usr/bin/env julia
#
# Hysing et al. 2009 rising-bubble benchmark, test case 1 — the Layer-1
# validation for VoF.jl's CSF surface tension (PLAN milestone 3).
#
# Physical setup (Hysing, Int. J. Numer. Meth. Fluids 60, 2009):
#   domain [0,1]×[0,2] m, bubble D=0.5 m at (0.5, 0.5)
#   surrounding fluid 1: ρ₁=1000, μ₁=10 ; bubble fluid 2: ρ₂=100, μ₂=1
#   g = 0.98, σ = 24.5  ⇒  Re = ρ₁√(g D³)/μ₁ = 35, Eo = ρ₁ g D²/σ = 10
#   run to t = 3 s.
#
# Published reference values (benchmark groups TP2D/FreeLEM/MooNMD):
#   min circularity  ≈ 0.9011–0.9013  (at t ≈ 1.9)
#   max rise velocity ≈ 0.2417–0.2421 (at t ≈ 0.92–0.93)
#   centroid y_c(3)   ≈ 1.0813–1.0817
#
# Known deviations of this driver from the benchmark BCs: WaterLily's
# domain BC is uBC-Dirichlet on all walls (benchmark uses free-slip on
# the vertical sides) — expected to bias the rise velocity slightly low.
# Bubble = LIGHT phase: α=0 inside the bubble, α=1 outside (VoF.jl
# water convention).
#
# ENV knobs: WL_N (cells across width, default 128), WL_TEND (default 3.0),
#   WL_CALPHA (default 1), WL_MULES (default true), WL_PASSES (default 4),
#   WL_TAG.
#
# Output: runs/hysing_bubble/bubble_<TAG>.csv
#   (t, y_c, v_c, circularity, area_rel, mass_rel)

using WaterLily
using VoF
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "hysing_bubble"))
mkpath(OUTDIR)

# --- Physical parameters (Hysing case 1) -------------------------------------
const W_p  = 1.0                 # domain width [m]
const H_p  = 2.0                 # domain height [m]
const D_p  = 0.5                 # bubble diameter [m]
const G_p  = 0.98
const σ_p  = 24.5
const ρ1   = 1000.0; const μ1 = 10.0    # surrounding (α=1)
const ρ2   = 100.0;  const μ2 = 1.0     # bubble (α=0)

const NX = parse(Int, get(ENV, "WL_N", "128"))
const NY = 2 * NX
const ΔX = W_p / NX

# --- Cell-unit scales ---------------------------------------------------------
const U_ref = sqrt(G_p * D_p)            # 0.7 m/s
const G_c   = G_p * ΔX / U_ref^2
const ν1_c  = (μ1 / ρ1) / (U_ref * ΔX)
const ν2_c  = (μ2 / ρ2) / (U_ref * ΔX)
const σ_c   = σ_p / U_ref^2              # see csf_force! docstring
const μ1_c  = ν1_c * ρ1                  # VoFFlow stores μ with ρ in
const μ2_c  = ν2_c * ρ2                  # physical numeric units

const T_END   = parse(Float64, get(ENV, "WL_TEND", "3.0"))
const SAMPLE  = parse(Float64, get(ENV, "WL_SAMPLE", "0.05"))
const PASSES  = parse(Int, get(ENV, "WL_PASSES", "4"))
const CALPHA  = parse(Float64, get(ENV, "WL_CALPHA", "1"))
const USE_MULES = parse(Bool, get(ENV, "WL_MULES", "true"))
const TAG     = get(ENV, "WL_TAG", "")
const T_NUM   = Float32

@info @sprintf("""Hysing case 1: N=%d×%d  ΔX=%.4f m  U_ref=%.2f m/s
  Re=%.0f  Eo=%.0f  σ_c=%.1f  g_c=%.5f  ν1_c=%.3f  ν2_c=%.3f
  reference: c_min≈0.901, v_max≈0.242 @ t≈0.92, y_c(3)≈1.081""",
    NX, NY, ΔX, U_ref, ρ1*sqrt(G_p*D_p^3)/μ1, ρ1*G_p*D_p^2/σ_p,
    σ_c, G_c, ν1_c, ν2_c)

# --- VoF + Flow + Poisson ------------------------------------------------------
# Bubble (α=0) of radius R_c at (NX/2, NX/2) cell coords.
const R_c  = D_p / 2 / ΔX
const cx_c = 0.5 / ΔX
const cy_c = 0.5 / ΔX
α₀ = (i, x) -> (x[1]-cx_c)^2 + (x[2]-cy_c)^2 < R_c^2 ? 0f0 : 1f0

vof = VoFFlow((NX, NY);
    α₀ = α₀, ρ_w = ρ1, ρ_a = ρ2, μ_w = μ1_c, μ_a = μ2_c, T = T_NUM)
L0 = copy(vof.L)

flow = WaterLily.Flow((NX, NY), (T_NUM(0), T_NUM(0));
    T = T_NUM, ν = VoF.viscosity(vof),
    g = (i, x, t) -> i == 2 ? T_NUM(-G_c) : T_NUM(0), Δt = 0.1)
pois = WaterLily.MultiLevelPoisson(flow.p, L0, flow.σ;
    tol = T_NUM(1e-8), itmx = 200)
sim = (flow = flow, pois = pois)
st! = VoF.surface_tension(vof, T_NUM(σ_c); passes = PASSES)

# --- Bubble metrics (bubble fraction B = 1-α) ----------------------------------
function metrics(vof, flow)
    α = vof.α; u = flow.u
    A = 0.0; Sy = 0.0; Sv = 0.0; P = 0.0
    NXg, NYg = size(α)
    @inbounds for I in CartesianIndices((2:NXg-1, 2:NYg-1))
        B = 1 - α[I]
        if B > 0
            A  += B
            Sy += B * (I.I[2] - 1.5)
            vc = (u[I, 2] + u[I + WaterLily.δ(2, I), 2]) / 2
            Sv += B * vc
        end
        gx = (α[I + WaterLily.δ(1, I)] - α[I - WaterLily.δ(1, I)]) / 2
        gy = (α[I + WaterLily.δ(2, I)] - α[I - WaterLily.δ(2, I)]) / 2
        P += sqrt(gx^2 + gy^2)
    end
    y_c  = Sy / A * ΔX                       # [m]
    v_c  = Sv / A * U_ref                    # [m/s]
    circ = 2 * sqrt(π * A) / P               # π·d_a / P
    return y_c, v_c, circ, A
end

# --- Time integration -----------------------------------------------------------
const Δt_grav_cap = T_NUM(0.3 * sqrt(2 / G_c))
# Brackbill capillary timestep limit (cell units, Δx=1)
const Δt_cap = T_NUM(sqrt((ρ1 + ρ2) / 2 / (4π * σ_c)))

csv = joinpath(OUTDIR, isempty(TAG) ? "bubble.csv" : "bubble_$(TAG).csv")
io = open(csv, "w")
println(io, "t,y_c,v_c,circ,area_rel,mass_rel")

t_phys = 0.0; t_next = 0.0; nstep = 0
_, _, _, A0 = metrics(vof, flow)
m0 = sum(@view vof.α[2:NX+1, 2:NY+1])
v_max = 0.0; t_vmax = 0.0; c_min = 1.0
while t_phys < T_END
    sim.flow.Δt[end] = min(sim.flow.Δt[end], Δt_grav_cap, Δt_cap)
    WaterLily.mom_step!(sim.flow, sim.pois; udf = st!)
    dt_cell = sim.flow.Δt[end-1]
    if USE_MULES
        step_vof_mules!(vof, sim; dt = dt_cell, c_α = CALPHA)
    else
        step_vof!(vof, sim; dt = dt_cell, mass_repair = true)
    end
    global t_phys += dt_cell * ΔX / U_ref
    global nstep += 1
    if t_phys >= t_next
        y_c, v_c, circ, A = metrics(vof, flow)
        m = sum(@view vof.α[2:NX+1, 2:NY+1])
        if v_c > v_max; global v_max = v_c; global t_vmax = t_phys; end
        circ < c_min && (global c_min = circ)
        @printf io "%.4f,%.5f,%.5f,%.5f,%.6f,%.6f\n" t_phys y_c v_c circ A/A0 m/m0
        flush(io)
        global t_next += SAMPLE
    end
    if nstep % 250 == 0
        y_c, v_c, circ, _ = metrics(vof, flow)
        @info @sprintf("step=%5d t=%.3fs  y_c=%.4f  v_c=%.4f  circ=%.4f  α∈[%.2f,%.2f]  |u|=%.2f  niter=%d",
            nstep, t_phys, y_c, v_c, circ, minimum(vof.α), maximum(vof.α),
            maximum(abs, sim.flow.u), isempty(sim.pois.n) ? -1 : sim.pois.n[end])
    end
    if !isfinite(maximum(vof.α)) || maximum(abs, sim.flow.u) > 100
        @warn "Blew up at step $nstep, t=$t_phys"
        break
    end
end
close(io)
y_c, v_c, circ, _ = metrics(vof, flow)
@info @sprintf("""DONE (%d steps): y_c(%.2f)=%.4f  [ref 1.081]
  v_max=%.4f @ t=%.2f  [ref 0.242 @ 0.92]   c_min=%.4f  [ref 0.901]""",
    nstep, t_phys, y_c, v_max, t_vmax, c_min)
@info "wrote $csv"
