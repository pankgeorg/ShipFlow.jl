#!/usr/bin/env julia
#
# WaterLily damBreak using VoF.jl variable-density coupling.
#
# Setup matches OpenFOAM tutorials/incompressibleVoF/damBreak:
#   - 2D box, L_x = 0.585, L_y = 0.585 m
#   - Initial water column: 0 ≤ x ≤ 0.1461, 0 ≤ y ≤ 0.292
#   - ρ_water = 1000, ρ_air = 1, μ_water = 1e-3, μ_air = 1.8e-5
#   - Gravity (0, -9.81, 0) acting uniformly
#   - No surface tension
#
# Unit system (per HANDOFF.md §3):
#   length unit  = ΔX               (one cell)
#   velocity     = U_ref = √(g·H_w) (free-fall over column height)
#   time         = ΔX / U_ref
#   ⇒ g_cell    = g_phys · ΔX / U_ref²
#     ν_cell    = ν_phys / (U_ref · ΔX)
#   ρ is a dimensionless ratio (use ρ_w=1000, ρ_a=1).
#
# ENV knobs:
#   WL_N         grid cells per box side (default 64)
#   WL_TEND      physical end time in seconds (default 1.0)
#   WL_SAMPLE    sampling interval in seconds (default 0.05)
#   WL_RHO_RATIO ρ_w/ρ_a — try 10 first, then 1000 (default 10)
#   WL_POIS_TOL  Poisson L2 tolerance (default 1e-8)
#   WL_POIS_ITMX max V-cycles per projection (default 200)
#
# Output: runs/damBreak_waterlily/front_vs_t.csv  (t, x_front, τ, X)

using WaterLily
using VoF
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "damBreak_waterlily"))
mkpath(OUTDIR)

# --- Physical parameters ----------------------------------------------------
const L_BOX   = 0.585           # box side [m]
const COL_W_p = 0.1461          # water column width [m]
const COL_H_p = 0.292           # water column height [m]
const G_p     = 9.81            # gravity [m/s²]
const ν_w_p   = 1.0e-6          # water kinematic viscosity [m²/s]
const ν_a_p   = 1.8e-5          # air kinematic viscosity [m²/s]
const RHO_RATIO = parse(Float64, get(ENV, "WL_RHO_RATIO", "10"))

# --- Grid (cell units) ------------------------------------------------------
const N_RESOL = parse(Int, get(ENV, "WL_N", "64"))
const NX = N_RESOL
const NY = N_RESOL
const ΔX = L_BOX / NX                                # m / cell
const COL_W_c = COL_W_p / ΔX                         # column width  in cells
const COL_H_c = COL_H_p / ΔX                         # column height in cells

# --- Non-dim scales ---------------------------------------------------------
const U_ref = sqrt(G_p * COL_H_p)                    # m/s
const G_c   = G_p * ΔX / U_ref^2                     # gravity in cell-units
const ν_w_c = ν_w_p / (U_ref * ΔX)                   # cell-units ν_water
const ν_a_c = ν_a_p / (U_ref * ΔX)                   # cell-units ν_air

const ρ_w   = RHO_RATIO                              # dimensionless
const ρ_a   = 1.0
const μ_w_c = ρ_w * ν_w_c                            # VoFFlow uses μ/ρ internally
const μ_a_c = ρ_a * ν_a_c

const T_END        = parse(Float64, get(ENV, "WL_TEND",   "1.0"))
const SAMPLE_EVERY = parse(Float64, get(ENV, "WL_SAMPLE", "0.05"))
const POIS_TOL     = parse(Float64, get(ENV, "WL_POIS_TOL",  "1e-8"))
const POIS_ITMX    = parse(Int,     get(ENV, "WL_POIS_ITMX", "200"))
const TAG          = get(ENV, "WL_TAG", "")    # suffix for output filename

println("WaterLily damBreak — variable-density VoF (cell-units)")
@printf("  Grid          = %d × %d  (Δx = %.4e m)\n", NX, NY, ΔX)
@printf("  Column        = %.3f × %.3f m  (= %.1f × %.1f cells)\n",
        COL_W_p, COL_H_p, COL_W_c, COL_H_c)
@printf("  ρ_w/ρ_a       = %.1f  (ρ_w=%.1f, ρ_a=%.1f dimensionless)\n",
        ρ_w/ρ_a, ρ_w, ρ_a)
@printf("  U_ref         = %.4f m/s    g_cell = %.4f\n", U_ref, G_c)
@printf("  ν_w_cell      = %.3e        ν_a_cell = %.3e\n", ν_w_c, ν_a_c)
@printf("  Re_water      = %.0f          (= U·H/ν_water)\n",
        U_ref * COL_H_p / ν_w_p)
@printf("  Poisson tol   = %.1e        itmx = %d\n", POIS_TOL, POIS_ITMX)
@printf("  t_end         = %.2f s        sample = %.3f s\n", T_END, SAMPLE_EVERY)

# --- Initial α (water indicator) -------------------------------------------
# WaterLily uses cell-centred coords at x_cell = I.I .- 1.5  (so the centre
# of cell I=(2,2) is (0.5, 0.5) cell-units).  We're told the column is
# 0 ≤ x_phys ≤ COL_W_p, 0 ≤ y_phys ≤ COL_H_p; convert to cell-units.
α₀(_i, x_cell) = let
    (x_cell[1] ≤ COL_W_c && x_cell[2] ≤ COL_H_c) ? 1.0 : 0.0
end

const T_NUM = Float64

@info "Building VoFFlow + Flow + Poisson…"
vof = VoFFlow((NX, NY);
    α₀ = α₀,
    ρ_w = ρ_w, ρ_a = ρ_a,
    μ_w = μ_w_c, μ_a = μ_a_c,
    T = T_NUM,
)
L0 = copy(vof.L)

flow = WaterLily.Flow((NX, NY), (T_NUM(0), T_NUM(0));
    T  = T_NUM,
    ν  = vof.ν,                       # Hook 1: per-cell ν
    g  = (i, x, t) -> i == 2 ? T_NUM(-G_c) : T_NUM(0),
    Δt = 0.5,                         # initial Δt_cell — bootstrapped by CFL
)
pois = WaterLily.MultiLevelPoisson(flow.p, L0, flow.σ)
sim = (flow = flow, pois = pois)

# --- Diagnostics ------------------------------------------------------------

# Rightmost x (in metres) where α > threshold in the bottom y_band_cells rows
function front_position(α::AbstractArray{T,2}; threshold=0.5,
                        y_band_cells=3) where T
    nx, ny = size(α)
    best = 0.0
    @inbounds for j in 2:min(2 + y_band_cells, ny - 1)
        for i in 2:nx-1
            if α[i, j] > threshold
                x_phys = (i - 1.5) * ΔX
                x_phys > best && (best = x_phys)
            end
        end
    end
    return best
end

function water_mass_cells(α::AbstractArray{T,2}) where T
    s = zero(T)
    nx, ny = size(α)
    @inbounds for j in 2:ny-1, i in 2:nx-1
        s += α[i, j]
    end
    return s
end

# --- Time integration ------------------------------------------------------

function run!(sim, vof)
    @info "Stepping until t = $T_END s (= $(T_END * U_ref / ΔX) cell-time units)"
    csv = joinpath(OUTDIR,
        isempty(TAG) ? "front_vs_t.csv" : "front_vs_t_$(TAG).csv")
    io  = open(csv, "w")
    println(io, "t,x_front,tau,X")

    # Gravity-wave CFL in cell units: dt < CFL · √(2/g_cell). g_cell ≈ 0.06,
    # so dt < ~5.7. WaterLily's own CFL (convective) gives a different bound;
    # take the min.
    Δt_grav_cap = T_NUM(0.3 * sqrt(2 / G_c))

    t_cell  = 0.0
    t_phys  = 0.0
    t_next_sample = 0.0
    step_count = 0
    m0 = water_mass_cells(vof.α)

    while t_phys < T_END
        sim.flow.Δt[end] = min(sim.flow.Δt[end], Δt_grav_cap)

        WaterLily.mom_step!(sim.flow, sim.pois;
            pois_tol = T_NUM(POIS_TOL),
            pois_itmx = POIS_ITMX,
        )
        dt_cell = sim.flow.Δt[end-1]
        dt_phys = dt_cell * ΔX / U_ref
        if parse(Bool, get(ENV, "WL_MULES", "false"))
            step_vof_mules!(vof, sim; dt = dt_cell)
        else
            step_vof!(vof, sim; dt = dt_cell,
                mass_repair = parse(Bool, get(ENV, "WL_MASS_REPAIR", "false")))
        end
        t_cell += dt_cell
        t_phys += dt_phys
        step_count += 1

        if t_phys >= t_next_sample
            x_front = front_position(vof.α)
            τ = t_phys * sqrt(2 * G_p / COL_W_p)
            X = x_front / COL_W_p
            @printf io "%.4f,%.6f,%.4f,%.4f\n" t_phys x_front τ X
            flush(io)
            t_next_sample += SAMPLE_EVERY
        end

        if mod(step_count, 25) == 0 || step_count < 4
            u_max_cell = maximum(abs, sim.flow.u)
            mass = water_mass_cells(vof.α)
            @info @sprintf("step=%4d  t=%.4fs  Δt_phys=%.3e  x_f=%.3fm  α∈[%.3f,%.3f]  |u|_c=%.3f  m/m0=%.4f  niter=%d",
                step_count, t_phys, dt_phys,
                front_position(vof.α), minimum(vof.α), maximum(vof.α),
                u_max_cell, mass / m0,
                isempty(sim.pois.n) ? -1 : sim.pois.n[end])
        end
        if !isfinite(maximum(vof.α)) || maximum(abs, sim.flow.u) > 100
            @warn "Simulation blew up at step $step_count"
            break
        end
    end
    close(io)
    @info "wrote $csv"
end

run!(sim, vof)
