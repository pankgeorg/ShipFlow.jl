#!/usr/bin/env julia
#
# OpenFOAM-side spatial validation: time-average U across snapshots
# in [t_start, t_end] and sample on the same vertical lines / centerline
# / surface ring as the WaterLily profile script.
#
# Output: runs/cylinder_of_profiles/{ux_profile_x*D.csv, ux_centerline.csv, Cp_surface.csv}

using ShipFlow.Harness: read_volVectorField, read_volScalarField
using Printf
using Statistics: mean

const CASE   = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_fresh"))
const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_of_profiles"))
mkpath(OUTDIR)

const D    = 0.001          # cylinder diameter [m]
const Uinf = 1.5            # inlet speed       [m/s]

# Time-averaging window, in physical seconds. The shedding is well-developed
# by t ≈ 0.06 s in this case.
const T_AVG_START = 0.060
const T_AVG_END   = Inf

# Helpers — sort the time directories numerically.
function time_dirs(case)
    times = String[]
    for d in readdir(case)
        isdir(joinpath(case, d)) || continue
        try
            v = parse(Float64, d)
            v ≥ 0 && push!(times, d)
        catch
        end
    end
    return sort(times; by = s -> parse(Float64, s))
end

# Read U at every snapshot in [t_start, t_end], return cell-wise time-average.
function average_field(case::AbstractString, fieldname::AbstractString;
                       t_start::Real, t_end::Real, parser)
    times = time_dirs(case)
    selected = filter(t -> t_start ≤ parse(Float64, t) ≤ t_end, times)
    isempty(selected) && error("no time directories in window [$t_start, $t_end]")
    @info "Averaging $fieldname over $(length(selected)) snapshots: " *
          "$(first(selected)) → $(last(selected))"
    accum = nothing
    n_cells = nothing
    for t in selected
        path = joinpath(case, t, fieldname)
        isfile(path) || continue
        n, vals = parser(path)
        if accum === nothing
            accum = deepcopy(vals)
            n_cells = n
        else
            n == n_cells || error("inconsistent cell count at t=$t")
            for i in eachindex(accum)
                accum[i] = accum[i] .+ vals[i]
            end
        end
    end
    for i in eachindex(accum)
        accum[i] = accum[i] ./ length(selected)
    end
    return n_cells, accum
end

# Read cell-center coordinates (from 0/C, written by writeCellCentres).
function read_cell_centers(case)
    path = joinpath(case, "0", "C")
    isfile(path) || error("0/C not found — run `foamPostProcess -func writeCellCentres` first")
    n, vecs = read_volVectorField(path)
    xs = [v[1] for v in vecs]
    ys = [v[2] for v in vecs]
    zs = [v[3] for v in vecs]
    return (n=n, x=xs, y=ys, z=zs)
end

# Nearest-cell sampling. We work in 2D (z is uniform at z=0).
function nearest_cell_index(cc, xp, yp)
    best_d2 = Inf
    best_i  = 0
    @inbounds for i in 1:cc.n
        d2 = (cc.x[i] - xp)^2 + (cc.y[i] - yp)^2
        if d2 < best_d2
            best_d2 = d2
            best_i  = i
        end
    end
    return best_i
end

function main()
    cc = read_cell_centers(CASE)
    @info "Cell centers loaded" cc.n

    _, U_mean = average_field(CASE, "U"; t_start=T_AVG_START, t_end=T_AVG_END,
                              parser=read_volVectorField)
    _, p_mean = average_field(CASE, "p"; t_start=T_AVG_START, t_end=T_AVG_END,
                              parser=read_volScalarField)
    @info "Time-averaging done"

    # ─────────────────────────────────────────────────────────────────────
    # u_x(y) profiles at x = 2D, 5D, 10D downstream of the cylinder centre
    # ─────────────────────────────────────────────────────────────────────
    for xD in (2.0, 5.0, 10.0)
        xp = xD * D
        path = joinpath(OUTDIR, @sprintf("ux_profile_x%dD.csv", round(Int, xD)))
        open(path, "w") do io
            println(io, "y_over_D,u_over_U")
            for yD in -3.0:0.05:3.0
                yp = yD * D
                k = nearest_cell_index(cc, xp, yp)
                ux = U_mean[k][1] / Uinf
                @printf io "%.6f,%.6f\n" yD ux
            end
        end
        @info "wrote $path"
    end

    # ─────────────────────────────────────────────────────────────────────
    # wake centerline u_x(x)
    # ─────────────────────────────────────────────────────────────────────
    path = joinpath(OUTDIR, "ux_centerline.csv")
    open(path, "w") do io
        println(io, "x_over_D,u_over_U")
        for xD in 0.6:0.05:14.0
            xp = xD * D
            k = nearest_cell_index(cc, xp, 0.0)
            ux = U_mean[k][1] / Uinf
            @printf io "%.6f,%.6f\n" xD ux
        end
    end
    @info "wrote $path"

    # ─────────────────────────────────────────────────────────────────────
    # Cp(θ) on a ring just outside the cylinder surface (r = 0.55 D)
    # ─────────────────────────────────────────────────────────────────────
    path = joinpath(OUTDIR, "Cp_surface.csv")
    open(path, "w") do io
        println(io, "theta_deg,Cp")
        rp = 0.55 * D
        for k in 0:143
            θ = 2π * k / 144
            xp = rp * cos(θ)
            yp = rp * sin(θ)
            kc = nearest_cell_index(cc, xp, yp)
            cp = 2 * p_mean[kc] / Uinf^2
            @printf io "%.4f,%.6f\n" rad2deg(θ) cp
        end
    end
    @info "wrote $path"

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
