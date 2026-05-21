#!/usr/bin/env julia
#
# 2D cylinder Re=100 in WaterLily — matches the OpenFOAM
# `incompressibleFluid/cylinder` tutorial (with `Uinlet=1.5`, `D=0.001`,
# `nu=1.5e-5`) at the *non-dimensional* level. The output is a time
# series of drag and lift coefficients written to
#   <out_dir>/cylinder_waterlily.csv
# with columns `t  Cd  Cl`, plus a final summary printed to stdout.
#
# Run:
#   julia --project=. scripts/cylinder_Re100_waterlily.jl
#
# Parameters (non-dimensional):
#   Re      = 100
#   U       = 1
#   L (=D)  = N_resol = 64    (cells per diameter)
#   nu      = U L / Re        (chosen to give Re=100)
#   domain  = (16D, 8D)       (downstream-biased)
#   end     ≈ 250 tU/D        (~40 shedding cycles at St ≈ 0.165)
#
# Quantities of interest:
#   Cd   = 2 * F_x / (ρ U² D)   (with ρ=1, U=1, D=N_resol)
#   Cl   = 2 * F_y / (ρ U² D)
#   St   = f_shed * D / U       (extracted from Cl time series)

using WaterLily
using Printf

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_waterlily"))
mkpath(OUTDIR)

# ------------------------------------------------------------------------------
# Simulation set-up
# ------------------------------------------------------------------------------

function cylinder_simulation(; N_resol::Int=64, Re::Real=100,
                              n_diam_x::Int=16, n_diam_y::Int=8)
    # Domain dimensions in cells (powers of 2 for the multigrid).
    nx = n_diam_x * N_resol     # 16 * 64 = 1024
    ny = n_diam_y * N_resol     # 8  * 64 =  512

    # Cylinder placed 4 D downstream of the inlet, mid-height.
    cx, cy = 4 * N_resol, ny / 2

    radius = N_resol / 2
    sdf(x, t) = √sum(abs2, x .- (cx, cy)) - radius

    U = 1.0
    Simulation(
        (nx, ny),               # dims
        (U, 0),                 # uBC
        N_resol;                # length scale L = D
        ν = U * N_resol / Re,
        body = AutoBody(sdf),
        T = Float32,
    )
end

# ------------------------------------------------------------------------------
# Drag / lift coefficients
# ------------------------------------------------------------------------------

# WaterLily's `pressure_force` returns ∫ p · n dS with *outward* normal,
# i.e. the force ON THE FLUID from the body. The drag on the body is the
# negative of the x-component (Newton's third law). OpenFOAM's forceCoeffs
# reports +Cd for a body in a +x stream, so we negate to match conventions.
function coefficients(sim)
    fp = WaterLily.pressure_force(sim)
    fv = WaterLily.viscous_force(sim)
    F  = .-(fp .+ fv)                # force on body, +x stream → +drag
    s  = 0.5 * sim.L * sim.U^2
    return F[1] / s, F[2] / s
end

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

function main(; N_resol::Int=64, Re::Real=100, t_end::Real=250.0)
    @info "Setting up WaterLily cylinder simulation" N_resol Re

    sim = cylinder_simulation(; N_resol, Re)

    @info "Stepping to t = $t_end tU/D"

    # Save a CSV of (t, Cd, Cl) at every step.
    csv_path = joinpath(OUTDIR, "cylinder_waterlily.csv")
    io = open(csv_path, "w")
    println(io, "t,Cd,Cl")

    # Take one initial step before sampling, so that the BCs settle.
    sim_step!(sim)

    t_samples = Float64[]; Cd_samples = Float64[]; Cl_samples = Float64[]

    step_count = 0
    while sim_time(sim) < t_end
        sim_step!(sim; remeasure=false)
        step_count += 1
        if step_count % 5 == 0    # subsample to keep the file small
            Cd, Cl = coefficients(sim)
            t = sim_time(sim)
            push!(t_samples, t); push!(Cd_samples, Cd); push!(Cl_samples, Cl)
            @printf(io, "%.6f,%.6f,%.6f\n", t, Cd, Cl)
            flush(io)
            if step_count % 200 == 0
                @info @sprintf("t=%6.2f  Cd=%6.3f  Cl=%6.3f  (steps=%d)",
                               t, Cd, Cl, step_count)
            end
        end
    end
    close(io)

    @info "Wrote $csv_path with $(length(t_samples)) samples"

    # Quick scalar summary: average Cd over the second half (steady-state),
    # and the peak-to-peak / 2 amplitude of Cl as a Strouhal proxy via FFT.
    half = length(t_samples) ÷ 2
    Cd_mean = sum(@view Cd_samples[half:end]) / (length(Cd_samples) - half + 1)
    Cl_pp   = maximum(@view Cl_samples[half:end]) -
              minimum(@view Cl_samples[half:end])
    @printf("\n=== WaterLily cylinder Re=%d summary ===\n", Re)
    @printf("  Cd (mean, t>%.0f)   = %.4f\n", t_samples[half], Cd_mean)
    @printf("  Cl peak-to-peak     = %.4f\n", Cl_pp)
    @printf("  samples             = %d\n", length(t_samples))
    @printf("=========================================\n")

    return (t = t_samples, Cd = Cd_samples, Cl = Cl_samples)
end

# Allow `include`-able use (returns the named tuple) and `run-as-script`.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
