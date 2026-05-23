#!/usr/bin/env julia
#
# WaterLily damBreak using VoF.jl variable-density coupling.
#
# Setup matches OpenFOAM tutorials/incompressibleVoF/damBreak:
#   - 2D box, L_x = 0.585, L_y = 0.585 (3:1 cells)
#   - Initial water column: 0 ≤ x ≤ 0.1461, 0 ≤ y ≤ 0.292
#   - ρ_water = 1000, ρ_air = 1, μ_water = 1e-3, μ_air = 1.8e-5
#   - Gravity (0, -9.81, 0) acting uniformly
#   - No surface tension yet
#
# Output: runs/damBreak_waterlily/front_vs_t.csv

using WaterLily
using VoF
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "damBreak_waterlily"))
mkpath(OUTDIR)

# --- Parameters -------------------------------------------------------------
const N_RESOL = parse(Int,     get(ENV, "WL_N",    "64"))     # cells per L_box
const L_BOX   = 0.585                                          # box width [m]
const COL_W   = 0.1461                                         # water column width
const COL_H   = 0.292                                          # water column height
const ρ_W = 1000.0
const ρ_A = 1.0
const μ_W = 1e-3
const μ_A = 1.8e-5
const G   = 9.81
const T_END = parse(Float64, get(ENV, "WL_TEND", "1.0"))       # seconds
const SAMPLE_EVERY = parse(Float64, get(ENV, "WL_SAMPLE", "0.05"))

# --- Grid -------------------------------------------------------------------
# Match OF mesh proportions. OF blockMesh has 46 cells in x over 0.585 m
# (some refinement around the obstacle). We use a uniform grid: N_RESOL
# cells per L_BOX.
const NX = N_RESOL                                              # x cells
const NY = N_RESOL                                              # y cells
const ΔX = L_BOX / NX                                           # m / cell

# Length scale for WL nondim: L = N_RESOL (i.e., the whole box width in cells).
# Velocity scale: U = sqrt(g · COL_H) (free-fall over column height).
const U_SCALE = sqrt(G * COL_H)

# Bring physical (ρ, μ, g) into the same non-dim scaling.
# Pressure ~ ρ U², viscous ~ μ U / L. Use water as reference for ν.
# In WL the momentum is unit-density: Du/Dt = -∇p + g + ν∇²u (per unit mass).
# Variable density enters via the Poisson L = 1/ρ_local.
# To keep absolute magnitudes sane we feed:
#   ρ → ρ_phys (keeps L = 1/ρ in physical units)
#   μ → μ_phys
#   g → g_phys
# WL's `ν` is then μ/ρ_local. All in physical SI units. WL works in
# arbitrary units; consistency is the only requirement.

println("WaterLily damBreak — variable-density VoF")
@printf("  Grid          = %d x %d  (Δx = %.4e m)\n", NX, NY, ΔX)
@printf("  ρ_w, ρ_a      = %.1f, %.1f\n", ρ_W, ρ_A)
@printf("  μ_w, μ_a      = %.1e, %.1e\n", μ_W, μ_A)
@printf("  g             = %.2f m/s²\n", G)
@printf("  Column        = %.3f × %.3f m\n", COL_W, COL_H)
@printf("  t_end         = %.2f s\n", T_END)

# --- Build VoF + flow ------------------------------------------------------

# Cell-centred y position is (j - 1.5) Δx. Cells with both x and y inside
# the column are water.
α₀(_i, x_cell) = let
    x_m = x_cell[1] * ΔX
    y_m = x_cell[2] * ΔX
    (x_m ≤ COL_W && y_m ≤ COL_H) ? 1f0 : 0f0
end

const T_NUM = Float64

@info "Building VoFFlow + Flow + Poisson…"
vof = VoFFlow((NX, NY);
    α₀ = α₀,
    ρ_w = ρ_W, ρ_a = ρ_A, μ_w = μ_W, μ_a = μ_A,
    T = T_NUM,
)

# Initial L from current α (face-density weighted by VoFFlow ctor)
L0 = copy(vof.L)

flow = WaterLily.Flow((NX, NY), (T_NUM(0), T_NUM(0));
    T  = T_NUM,
    ν  = vof.ν,                      # PLAN 1 Hook 1
    g  = (i, x, t) -> i == 2 ? T_NUM(-G) : T_NUM(0),
    Δt = 1e-5,                       # tiny initial dt to bootstrap
)
pois = WaterLily.MultiLevelPoisson(flow.p, L0, flow.σ)

# Pack into a NamedTuple sim — step_vof! only needs sim.flow and sim.pois.
sim = (flow = flow, pois = pois)

# --- Time integration ------------------------------------------------------

function front_position(α::AbstractArray{T,2}; threshold=0.5f0,
                        y_band_cells=3) where T
    # Rightmost x (in physical metres) where α > threshold in the bottom y_band_cells.
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

function run!(sim, vof)
    @info "Stepping until t = $T_END s"
    csv = joinpath(OUTDIR, "front_vs_t.csv")
    io  = open(csv, "w")
    println(io, "t,x_front,tau,X")

    # WaterLily's CFL only considers convective velocity, not body forces.
    # For damBreak gravity-driven flow we must impose a gravity-wave CFL:
    #   Δt < sqrt(2 ΔX / g)  in physical s (CFL=0.5 → Δt_phys ≈ 0.022 s here).
    # In WaterLily cell-units (Δt_cell = Δt_phys / ΔX): Δt_cell ≈ 2.4.
    Δt_cell_cap = 0.1 * sqrt(2 * ΔX / G) / ΔX   # conservative gravity-wave CFL

    t_phys = 0.0
    t_next_sample = 0.0
    step_count = 0
    while t_phys < T_END
        # Cap the WaterLily-computed Δt before it's used by mom_step!
        sim.flow.Δt[end] = min(sim.flow.Δt[end], T_NUM(Δt_cell_cap))
        WaterLily.mom_step!(sim.flow, sim.pois)
        dt = sim.flow.Δt[end-1]
        dt_phys = dt * ΔX
        step_vof!(vof, sim; dt = dt)
        t_phys += dt_phys
        step_count += 1

        if t_phys >= t_next_sample
            x_front = front_position(vof.α)
            τ = t_phys * sqrt(2G / COL_W)
            X = x_front / COL_W
            @printf io "%.4f,%.6f,%.4f,%.4f\n" t_phys x_front τ X
            flush(io)
            t_next_sample += SAMPLE_EVERY
        end
        if mod(step_count, 5) == 0 || step_count < 6
            x_front = front_position(vof.α)
            u_max = maximum(abs, sim.flow.u)
            @info @sprintf("step=%d  t=%.4fs  Δt_phys=%.3e  x_f=%.4fm  α_max=%.3f  |u|_max=%.4f m/s",
                            step_count, t_phys, dt_phys, x_front,
                            maximum(vof.α), u_max)
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
