#!/usr/bin/env julia
#
# Extract the water-front position vs time from OpenFOAM damBreak,
# for cross-validation against VoF.jl (when density coupling lands)
# and the Martin–Moyce 1952 experimental data.
#
# "Front position" = rightmost x where α_water > 0.5 along the bottom
# row of cells (y < 0.05).  Output: runs/damBreak_front/front_vs_t.csv

using ShipFlow.Harness: read_volVectorField, read_volScalarField
using Printf

const CASE   = abspath(joinpath(@__DIR__, "..", "runs", "damBreak"))
const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "damBreak_front"))
mkpath(OUTDIR)

# Initial water column from system/setFieldsDict
const COLUMN_W = 0.1461   # m
const COLUMN_H = 0.292    # m

# Martin–Moyce 1952 non-dimensionalisation
# τ = t · sqrt(2g/L) where L = column width;  X = x_front / L
g_grav = 9.81

function time_dirs(case)
    times = String[]
    for d in readdir(case)
        isdir(joinpath(case, d)) || continue
        try
            v = parse(Float64, d); v ≥ 0 && push!(times, d)
        catch; end
    end
    return sort(times; by = s -> parse(Float64, s))
end

function read_cell_centers(case)
    _, vecs = read_volVectorField(joinpath(case, "0", "C"))
    return [v[1] for v in vecs], [v[2] for v in vecs], [v[3] for v in vecs]
end

function front_at(case, t, xs, ys; y_band=0.05)
    # Front = rightmost x at which the water column extends above y_band.
    # We define "water present" as α > 0.5 in ANY cell with y < y_band.
    path = joinpath(case, t, "alpha.water")
    isfile(path) || return missing
    _, α = read_volScalarField(path)
    best_x = 0.0
    for k in eachindex(α)
        ys[k] < y_band || continue
        α[k] > 0.5 || continue
        xs[k] > best_x && (best_x = xs[k])
    end
    return best_x
end

function main()
    @info "Reading damBreak…"
    xs, ys, zs = read_cell_centers(CASE)
    @info "Cells: $(length(xs)); x ∈ [$(minimum(xs)), $(maximum(xs))]"

    times = time_dirs(CASE)
    @info "Time directories: $(length(times)); range $(first(times)) → $(last(times))"

    csv = joinpath(OUTDIR, "front_vs_t.csv")
    open(csv, "w") do io
        println(io, "t,x_front,tau,X")
        for t in times
            tval = parse(Float64, t)
            x_front = front_at(CASE, t, xs, ys)
            x_front === missing && continue
            τ = tval * sqrt(2 * g_grav / COLUMN_W)
            X = x_front / COLUMN_W
            @printf io "%.4f,%.6f,%.4f,%.4f\n" tval x_front τ X
        end
    end
    @info "wrote $csv"

    # Print a brief summary
    println("\n--- damBreak water front advance ---")
    println("t [s]    x_front [m]   τ=t√(2g/L)   X=x/L")
    for t in times
        tval = parse(Float64, t)
        x_front = front_at(CASE, t, xs, ys)
        x_front === missing && continue
        τ = tval * sqrt(2 * g_grav / COLUMN_W)
        X = x_front / COLUMN_W
        @printf "%-8.3f %-13.4f %-12.3f %-6.3f\n" tval x_front τ X
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
