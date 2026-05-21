#!/usr/bin/env julia
#
# WaterLily cylinder Re=100 — produce spatial validation data:
#   * time-averaged u_x(y) profiles at x = 2D, 5D, 10D downstream
#   * time-averaged wake centerline u_x(x)
#   * surface pressure coefficient Cp(θ) on the cylinder
#
# Outputs columnar CSVs into runs/cylinder_waterlily_profiles/.

using WaterLily
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_waterlily_profiles"))
mkpath(OUTDIR)

const N_RESOL = 32      # cells per D
const RE      = 100
const N_DIAM_X = 16
const N_DIAM_Y = 8
const T_AVG_START = 60.0      # start averaging after wake has developed
const T_END       = 200.0     # ~30 shedding cycles in the averaging window

# Cylinder placement (matches scripts/cylinder_Re100_waterlily.jl)
const CX = 4 * N_RESOL                # x_index of cylinder centre
const CY = N_DIAM_Y * N_RESOL / 2     # y_index
const R  = N_RESOL / 2                # radius in cells

function cylinder_sim()
    nx, ny = N_DIAM_X * N_RESOL, N_DIAM_Y * N_RESOL
    sdf(x, t) = √sum(abs2, x .- (CX, CY)) - R
    U = 1.0
    Simulation((nx, ny), (U, 0), N_RESOL;
        ν = U * N_RESOL / RE,
        body = AutoBody(sdf),
        T = Float32,
    )
end

# ──────────────────────────────────────────────────────────────────────────
# Spatial samplers (on the cell-centered, staggered Cartesian grid)
# ──────────────────────────────────────────────────────────────────────────

# u_x at a "physical" location (xpos, ypos) measured in cell units from
# the inlet/bottom corner. Returns u at the x-face (since u[I,1] sits at
# the lower-x face of cell I).
@inline function ux_at(u::AbstractArray{T,3}, xpos::Real, ypos::Real) where T
    # u[i, j, 1] is u_x at face x = i - 1.5 in cell-centre coords. Sample
    # by linear interpolation in (xpos, ypos).
    nx, ny, _ = size(u)
    fx = xpos + 1.5
    fy = ypos + 1.5
    i0 = clamp(floor(Int, fx), 1, nx-1)
    j0 = clamp(floor(Int, fy), 1, ny-1)
    a = fx - i0
    b = fy - j0
    @inbounds (
        (1-a)*(1-b)*u[i0,j0,1] + a*(1-b)*u[i0+1,j0,1] +
        (1-a)*b*u[i0,j0+1,1]   + a*b*u[i0+1,j0+1,1]
    )
end

@inline function p_at(p::AbstractArray{T,2}, xpos::Real, ypos::Real) where T
    nx, ny = size(p)
    fx = xpos + 1.5
    fy = ypos + 1.5
    i0 = clamp(floor(Int, fx), 1, nx-1)
    j0 = clamp(floor(Int, fy), 1, ny-1)
    a = fx - i0
    b = fy - j0
    @inbounds (
        (1-a)*(1-b)*p[i0,j0] + a*(1-b)*p[i0+1,j0] +
        (1-a)*b*p[i0,j0+1]   + a*b*p[i0+1,j0+1]
    )
end

# ──────────────────────────────────────────────────────────────────────────
# Run + collect MeanFlow + emit profile CSVs
# ──────────────────────────────────────────────────────────────────────────

function main()
    @info "Setting up WaterLily cylinder for profile sampling" N_RESOL RE
    sim = cylinder_sim()
    mean = WaterLily.MeanFlow(sim.flow; t_init=0.0)

    @info "Stepping to t = $T_END tU/D; averaging from t = $T_AVG_START"
    while sim_time(sim) < T_END
        sim_step!(sim; remeasure=false)
        if sim_time(sim) >= T_AVG_START
            WaterLily.update!(mean, sim.flow)
        end
    end

    @info "Done stepping. Sampling profiles…"

    Umean = mean.U
    Pmean = mean.P

    # Output coordinates in physical D-units, measured from the cylinder centre.
    function dump_line_y(filepath, x_phys, ys)
        open(filepath, "w") do io
            println(io, "y_over_D,u_over_U")
            for y_phys in ys
                xpos = CX + x_phys * N_RESOL
                ypos = CY + y_phys * N_RESOL
                ux = ux_at(Umean, xpos, ypos)
                @printf io "%.6f,%.6f\n" y_phys ux
            end
        end
    end

    function dump_line_x(filepath, xs, y_phys)
        open(filepath, "w") do io
            println(io, "x_over_D,u_over_U")
            for x_phys in xs
                xpos = CX + x_phys * N_RESOL
                ypos = CY + y_phys * N_RESOL
                ux = ux_at(Umean, xpos, ypos)
                @printf io "%.6f,%.6f\n" x_phys ux
            end
        end
    end

    function dump_cp(filepath; n=72)
        # Cp(θ) on a thin ring just outside the cylinder. Use r = R+1 (one
        # cell off the surface) to stay outside the BDIM kernel.
        open(filepath, "w") do io
            println(io, "theta_deg,Cp")
            r_phys = (R + 1) / N_RESOL   # in D-units; r ≈ 0.5D + 1 cell
            U2 = 1.0      # 0.5 * ρ * U^2 = 0.5 in non-dim (we drop the 0.5)
            for k in 0:n-1
                θ = 2π * k / n
                xpos = CX + r_phys * N_RESOL * cos(θ)
                ypos = CY + r_phys * N_RESOL * sin(θ)
                p_local = p_at(Pmean, xpos, ypos)
                cp = 2 * p_local / U2
                @printf io "%.4f,%.6f\n" rad2deg(θ) cp
            end
        end
    end

    # u_x(y) profiles
    ys = -3.0:0.05:3.0
    for x_phys in (2.0, 5.0, 10.0)
        path = joinpath(OUTDIR, @sprintf("ux_profile_x%dD.csv", round(Int, x_phys)))
        dump_line_y(path, x_phys, ys)
        @info "wrote $path"
    end

    # wake-centerline u_x(x). Capped at 11D so we stay inside the WL
    # domain (cylinder at 4D from inlet, domain extends to 12D downstream
    # in the n_diam_x=16 case).
    xs = 0.6:0.05:11.0
    path = joinpath(OUTDIR, "ux_centerline.csv")
    dump_line_x(path, xs, 0.0)
    @info "wrote $path"

    # Cp(θ) on cylinder surface
    path = joinpath(OUTDIR, "Cp_surface.csv")
    dump_cp(path; n=144)
    @info "wrote $path"

    # Summary
    @info "Time-averaging window: $(mean.t[1]) → $(mean.t[end])"
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
