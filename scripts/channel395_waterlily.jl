#!/usr/bin/env julia
#
# WaterLily channel-flow LES at Re_τ ≈ 395 — Turbulence.jl Smagorinsky
# validation against OpenFOAM channel395.
#
# Setup:
#   - 3D periodic streamwise/spanwise channel
#   - Walls in y modelled via BDIM (two flat plates as an AutoBody SDF)
#   - Constant-mass-flow driver: udf adds a streamwise body force tuned
#     each step to hold ⟨u_x⟩ = U_bar
#   - Smagorinsky eddy viscosity refreshed each step into flow.ν
#   - MeanFlow time-average over t ∈ [t_avg_start, t_end]
#   - Output: u_x(y) profile averaged over (x, z), plus u_τ, u⁺(y⁺).

using WaterLily
using Turbulence
using Printf

# ───────────────────────────── parameters ─────────────────────────────
const N_HC     = parse(Int,    get(ENV, "WL_N_HC",   "32"))   # cells per channel half-height δ
const RE_TAU   = parse(Float64, get(ENV, "RE_TAU",   "395"))
const RE_BULK  = parse(Float64, get(ENV, "RE_BULK",  "6675"))  # bulk Re — matches OF ν=2e-5, Ubar=0.1335, δ=1
const U_BAR    = parse(Float64, get(ENV, "U_BAR",    "1.0"))   # bulk mean streamwise velocity (non-dim)
# Fixed body force g_x = u_τ²/δ that drives the channel at the target Re_τ.
# u_τ = U_bulk/16.9 (channel scaling); δ = N_HC cells.  In WL non-dim with
# L=δ, U=U_bar, the resulting g_x has dimensions of accel.
const U_TAU    = U_BAR / 16.9
const G_X      = U_TAU^2 / N_HC
const T_END    = parse(Float64, get(ENV, "WL_TEND",  "400.0")) # tU_bar/δ
const T_AVG    = parse(Float64, get(ENV, "WL_TAVG",  "100.0"))
const N_X      = parse(Int,    get(ENV, "N_X",      "128"))    # streamwise cells
const N_Z      = parse(Int,    get(ENV, "N_Z",      "64"))     # spanwise cells
const C_S      = parse(Float64, get(ENV, "C_S",     "0.17"))
const OUT_SUFFIX = get(ENV, "WL_SUFFIX", "")

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs",
                                 "channel395_waterlily" * OUT_SUFFIX))
mkpath(OUTDIR)

# Derived:
const N_Y    = 2 * N_HC         # full channel height in cells
const DELTA  = N_HC             # δ in cells (used as L scale)
# Set ν so that Re_bulk = U_bar · δ / ν = RE_BULK. For Re_τ=395 channel:
# Re_bulk ≈ 6675 (Moser-Kim-Mansour 1999 scaling; U_bulk/u_τ ≈ 16.9).
const NU_GRID = DELTA / RE_BULK

println("Channel parameters:")
@printf("  N_HC          = %d (half-channel cells)\n", N_HC)
@printf("  N_Y           = %d (full channel cells)\n", N_Y)
@printf("  N_X, N_Z      = %d, %d\n", N_X, N_Z)
@printf("  ν (grid units) = %.4e\n", NU_GRID)
@printf("  Smagorinsky Cs = %.3f\n", C_S)
@printf("  u_τ (target)   = %.4f\n", U_TAU)
@printf("  g_x (body force) = %.4e\n", G_X)
@printf("  t_end / tU/δ   = %.1f\n", T_END)
@printf("  t_avg_start    = %.1f\n", T_AVG)
@printf("  Output         = %s\n", OUTDIR)

# ─────────────────────── channel body (two walls) ─────────────────────
# Walls at y = 0 (mesh cell-index 1.5) and y = N_Y (mesh cell-index N_Y + 0.5).
# Fluid is 1.5 < y_cell < N_Y + 0.5. The body lies *outside* this band.
# SDF: signed distance to the nearest wall, positive inside fluid.

function channel_body()
    function sdf(x, t)
        # x[2] is cell-coord; channel inside is 0 ≤ x[2] ≤ N_Y.
        y = x[2]
        return min(y, N_Y - y)
    end
    AutoBody(sdf)
end

# ─────────────────────── mass-flow driver ─────────────────────────────
# Maintain ⟨u_x⟩ = U_BAR by adding a uniform streamwise body force
# proportional to the current mean-velocity deficit. Tuned each step
# with a relaxation factor so we don't ring.

mutable struct MassFlowDriver
    target   :: Float64
    g_x      :: Float64    # current body-force estimate
    relax    :: Float64
end

MassFlowDriver(target) = MassFlowDriver(target, 0.0, 0.05)   # small relax — slower integral control

function (drv::MassFlowDriver)(flow, t; kwargs...)
    # Compute current bulk u_x over the interior (skip ghost cells)
    u_mean = mean_interior_ux(flow.u)
    # Update body force estimate (proportional control)
    deficit = drv.target - u_mean
    drv.g_x += drv.relax * deficit
    # Add to flow.f everywhere (uniform body force)
    @inbounds for I in CartesianIndices(@view flow.f[:,:,:,1])
        flow.f[I, 1] += drv.g_x
    end
    return nothing
end

function mean_interior_ux(u)
    nx, ny, nz, _ = size(u)
    s = 0.0; n = 0
    @inbounds for k in 2:nz-1, j in 2:ny-1, i in 2:nx-1
        s += u[i, j, k, 1]
        n += 1
    end
    return s / n
end

# ──────────────────────────── simulation ──────────────────────────────

function build()
    body  = channel_body()
    model = Smagorinsky((N_X, N_Y, N_Z); Cs=Float32(C_S), ν₀=Float32(NU_GRID))

    # IC: parabolic mean + Schoppa-Hussain-style streamwise vortex pair
    # perturbation at 20 % of U_BAR. The pair injects O(1) wall-normal
    # vorticity which is the key inflectional driver for transition.
    # Uses deterministic seeding (no rand()) so the run is reproducible.
    uλ = (i, x) -> begin
        η = (x[2] - N_HC) / N_HC          # -1 at walls, 0 at centre
        # streamwise mean (parabolic)
        u_mean = 1.5 * U_BAR * (1 - η^2)
        # broadband 3D wave packet — multiple modes excited together
        ξx = 2π * x[1] / N_X
        ξy = π  * x[2] / N_Y               # half-wavelength in y (no-slip at walls)
        ξz = 2π * x[3] / N_Z
        # Streamwise vortex pair (u, v, w consistent ish — not exactly
        # divergence-free, but the projection step fixes that)
        amp = 0.20 * U_BAR
        if i == 1     # u-component: streamwise streaks
            u_pert = amp * (1 - η^2) * (
                0.6 * cos(ξx) * cos(2ξz) +
                0.4 * cos(2ξx + 1.7) * cos(ξz + 0.5) +
                0.3 * sin(3ξx) * cos(3ξz)
            )
            Float32(u_mean + u_pert)
        elseif i == 2  # v-component: wall-normal kicks
            v_pert = 0.5 * amp * sin(ξy) * (
                sin(ξx) * sin(2ξz) +
                0.5 * sin(2ξx) * sin(ξz)
            )
            Float32(v_pert)
        else           # w-component
            w_pert = 0.5 * amp * sin(ξy) * (
                cos(ξx) * sin(ξz) +
                0.4 * sin(3ξx + 0.8) * cos(ξz)
            )
            Float32(w_pert)
        end
    end

    # Constant streamwise body force — the standard channel-flow driver.
    g = (i, x, t) -> i == 1 ? Float32(G_X) : 0f0

    sim = Simulation(
        (N_X, N_Y, N_Z),
        (0f0, 0f0, 0f0),
        DELTA;
        U    = Float32(U_BAR),
        uλ   = uλ,
        ν    = model.ν,
        body = body,
        g    = g,
        perdir = (1, 3),
        T = Float32,
    )
    return sim, model
end

# ──────────────────────────── output ──────────────────────────────────

"Average u_x over (x, z) at every y-cell, returns vector of length N_Y."
function ux_profile_y(U::AbstractArray{T,4}) where T
    nx, ny, nz, _ = size(U)
    out = zeros(Float64, ny - 2)   # interior only
    @inbounds for j in 2:ny-1
        s = 0.0; n = 0
        for k in 2:nz-1, i in 2:nx-1
            s += U[i, j, k, 1]
            n += 1
        end
        out[j-1] = s / n
    end
    return out
end

# ─────────────────────────────── main ─────────────────────────────────

function main()
    @info "Building WaterLily channel sim…"
    sim, model = build()
    mean = WaterLily.MeanFlow(sim.flow; t_init=0.0)

    @info "Stepping; SGS via Smagorinsky; fixed body force g_x=$G_X"
    step_count = 0
    while sim_time(sim) < T_END
        sim_step!(sim;
            remeasure = false,
            udf = (flow, t; kwargs...) -> model(flow, t),
        )
        step_count += 1
        if sim_time(sim) >= T_AVG
            WaterLily.update!(mean, sim.flow)
        end
        if mod(step_count, 50) == 0
            u_mean = mean_interior_ux(sim.flow.u)
            @info @sprintf("t=%6.1f  ⟨u⟩=%.4f  Δt=%.3e  νₜ_max=%.3e",
                           sim_time(sim), u_mean, sim.flow.Δt[end],
                           maximum(model.ν) - NU_GRID)
        end
    end

    @info "Sampling u_x(y) profile averaged over (x, z) and time…"
    profile = ux_profile_y(mean.U)

    csv = joinpath(OUTDIR, "ux_profile.csv")
    open(csv, "w") do io
        println(io, "y_over_delta,u_over_Ubar")
        for (j, u) in pairs(profile)
            y_cell = j         # 1..N_Y in cell index after the offset
            y_over_delta = (y_cell - 0.5) / DELTA - 1.0
            # y/δ measured from channel center; centerline at j ≈ N_HC.
            @printf io "%.4f,%.6f\n" y_over_delta (u / U_BAR)
        end
    end
    @info "wrote $csv"

    # u_τ estimate: from wall-derivative
    Δy = 1.0
    u_first = profile[1]
    u_τ = sqrt(NU_GRID * u_first / Δy)
    Re_τ_est = u_τ * DELTA / NU_GRID
    @printf("\n=== WaterLily channel summary ===\n")
    @printf("  ⟨u⟩ (final)  = %.4f\n", mean_interior_ux(sim.flow.u))
    @printf("  u_τ estimate = %.4f\n", u_τ)
    @printf("  Re_τ est.    = %.1f  (target %.1f)\n", Re_τ_est, RE_TAU)
    @printf("=================================\n")

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
